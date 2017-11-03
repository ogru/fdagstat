# Kriging functions:
#' @useDynLib fdagstat
#' @importFrom Rcpp sourceCpp

#' @title Creates an \code{fstat} structure.
#'
#' @param g           - An existing fstat structure (optional). Defaults to \code{NULL}.
#' @param vName       - Variable name. Note this will be the name used in all analyses.
#' @param Coordinates - All coordinates/covariates associated with the functional data (a \code{data.frame}).
#' @param Functions   - A data frame with functional data stacked in columns.
#' @param scalar      - Scalars or Functions? (default: \code{scalar = FALSE}). Note: this will determine forecasting, variograms, etc.
#'
#' @details Returns an \code{fstat} structure which is a simple \code{R} list consisting of the following slots:
#' a) \code{data} - a list with all data, each element has two slots: coordinates and functions;
#' b) \code{variograms} - an empty slot reserved for empirical variograms computed with \code{fvariogram}
#' c) \code{models}     - an empty slot reserved for variogram models computed with \code{fit.variogram} or \code{fit.lmc}
#'
#' When adding functional data, the code will check whether covariates passed in coordinates data frame match the covariates
#' in the coordinates data frame that is already present at slot \code{$data[[1]]$coordinates}. If they don't an error will pop up.
#' Extra covariates are removed without prompt. You've been warned!
#'
#' Note that if parameter scaling (i.e. unit cube) is required, that should be performed as a part of pre-processing.
#' see \code{help(scaleInput)} function of this package.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu}), Alessandra Menafoglio (\email{alessandra.menafoglio@polimi.it})
#'
#' @export
#'
fstat <- function(g = NULL, vName, Coordinates, Functions, scalar = FALSE){

  if(is.null(vName)) stop('Variable name has to be provided. Exiting!')
  if(is.null(Coordinates)) stop('The coordinates were not provided. Exiting!')
  if(is.null(Functions)) stop('The functions were not provided. Exiting!')

  # Checks:
  if(!is.data.frame(Coordinates)) stop('The coordinates should be provided in a data frame. Exiting!')
  if(!is.data.frame(Functions)) stop('The functions should be provided in a data frame. Exiting!')

  if(nrow(Coordinates)!=ncol(Functions) & !scalar) stop('The number of coordinates is not equal to the number of functions. Exiting!')

  if(!is.null(g)){
    if(attr(g$data[[1]], 'isScalar') != scalar){
      stop("Mixed data is not allowed! Either all variables are functional or they are all scalar.")
    }
  }

  if(is.null(g)){ # create
    g <- vector('list', 3)
    names(g)      <- c('data', 'variogram', 'model')
    g$data        <- vector('list', 1)
    g$data[[1]]   <- vector('list', 2)
    names(g$data) <- vName

  } else { # append

    g$data[[length(g$data) + 1]]  <- vector('list', 2)
    names(g$data)[length(g$data)] <- vName

  }

  # need to check if added coordinates are the same as already present coordinates.
  # we cannot mix apples and oranges.

  names(g$data[[vName]])      <- c('coordinates', 'functions')
  g$data[[vName]]$coordinates <- Coordinates
  g$data[[vName]]$functions   <- Functions

  # Attribute setting for various purposes.
  attr(g$data[[vName]], 'hasDrift') <- FALSE


  attr(g$data[[vName]], 'isScalar')   <- scalar

  class(g) <- 'fstat'

  return(g)

}

#' @title Computes empirical trace-variograms and trace-cross-variograms
#'
#' @param formula - Formula specifying which coordinates to use in variogram computation.
#' @param g - Object of class \code{fstat} as produced by \code{fstat} function
#' @param Nlags - Number of variogram lags (required)
#' @param LagMax - Maximum lag distance (also required)
#' @param LagTolerance - Optional. If not provided the code will set it to LagStep/2
#' @param ArgStep - Discretization step of the functional data. Used for numerical integration.
#'                  Note that the code assumes the same sampling for all functional data in contained within g.
#' @param directions - It can be: 'omni' when it will compute omni-directional variogram (default),
#'                                'all' when it computes variograms along each input dimension,
#'                                an \code{nPxd} matrix where \code{nP} = the number of coordinates, \code{d} = the number of directions.
#'                                The code will normalize all vectors!
#' @param angTolerance - Tolerance for computing directional variograms
#' @param crossMode    - Either \code{conventional} or \code{pseudo}, for conventional or pseudo trace-cross-variogram computation,
#'                       set this value to "NA" if you want to skip trace-cross-variogram computation all together.
#' @param useResidual  - This is applicable only when the drift term was already estimated with \code{estimateDrift} function.
#' @param comments     - Whether or not to provide runtime comments!
#'
#' @details Returns an augmented \code{fstat} version of \code{g} by inserting a \code{gstatVariogram} slot.
#'          If there are NA's in functional data, then the code will not compute regular trace variograms. Instead, it will use an
#'          approximation proposed by Gromenko and Kokoszka.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu}) and Alessandra Menafoglio (\email{alessandra.menafoglio@poilimi.it}).
#'
#' @export
#'
fvariogram <- function(formula      = NULL,
                       g            = NULL,
                       Nlags        = NULL,
                       LagMax       = NULL,
                       LagTolerance = NULL,
                       ArgStep      = NULL,
                       directions   = 'omni',
                       angTolerance = 0,
                       crossMode    = 'conventional',
                       useResidual  = FALSE,
                       comments     = TRUE){

  # Checkups: ################################################################
  if(class(g) != 'fstat' | is.null(g)) stop('The function requires objects of class "fstat". Exiting!')

  if((is.null(Nlags)))    stop('The number of variogram lags was not provided! Exiting')

  if(is.null(LagMax))     stop('Maximum lag distance was not provided! Exiting!')

  if(!attr(g$data[[1]], 'isScalar')){
      if(is.null(ArgStep))    stop('Argument step was not provided. Exiting!')
  } else {
    ArgStep = 1
  }

  crossMode <- match.arg(crossMode, c("conventional", "pseudo", "NA"))

  if(comments){
    if(useResidual && attr(g$data[[1]], 'hasDrift')){
      print('Drift was found! Variograms will be evaluated on the residuals.')
    } else {
      print('Variograms will be evaluated on the raw data')
    }
  }

  # This will throw an error if formula is wrong.
  .ModelFrame <- model.frame(formula, data = g$data[[1]]$coordinates)

  # Figure out directions:  ##################################################
  dirNames = NULL

  if(is.matrix(directions)) {

    if(any(apply(directions, 2, function(x) sqrt(sum(x^2))) == 0)){
      stop(paste('One of the provided directions is a zero vector. This is not allowed!',
                 'fvariogram will now exit. Please fix that and try again!', sep = " ")
      )
    }

    if(nrow(directions) != ncol(.ModelFrame)) {
      stop('The number of rows in directions matrix is not equal to the number of selected coordinates/covariates. Exiting!')
    }

    # normalize (default behavior):
    directions <- apply(directions, 2, function(x) x/sqrt(sum(x^2)))

    if(is.null(colnames(directions))){
      dirNames   = c('omni', paste('dir', 1:ncol(directions), sep = ""))
    } else {
      dirNames   = c("omni", colnames(directions))
    }

    directions <- cbind(0, directions)

  } else {

    if(directions == 'omni') { # dummy variable to fool Cpp code.
      directions = matrix(0, ncol = 1, nrow = ncol(.ModelFrame))
      dirNames   = 'omni'

    } else {
      if(directions == 'all') { # directions matrix, the first column is for omni directional.
        directions = cbind(0, diag(ncol(.ModelFrame)))
        dirNames   = c('omni', names(.ModelFrame))
      }
    }

  }

  # Work on variogram lags: ##################################################
  Lags    <- seq(0, LagMax, length.out = Nlags)
  lagStep <- (Lags[3] - Lags[2])

  if(is.null(LagTolerance)){
    LagTolerance <- lagStep/2
  } else {
    if(LagTolerance > lagStep/2){
      LagTolerance = lagStep/2
      if(comments){
        print(
          paste('The provided lag tolerance is higher than lag half. Adjusting to: ',
                round(LagTolerance, 6), sep = ""))
      }
    }
  }

  #Avoid zero lag! gstat doesn't like that!
  Lags[1] <- lagStep/2

  # This is the main part that calls Cpp routines: ###########################
  # First it computes trace variograms for each variable, then
  # it computes cross variograms where it distinguishes between two types:
  # conventional and pseudo. The former requires isotopic coordinates.
  # Code will throw an error if the number of isotopic coordinates is small,
  # and it will recommend to switch to pseudo cross variogram.
  ############################################################################

  for(i in 1:length(g$data)){

    vName         <- names(g$data)[i]

    # Apply formula to coordinates:
    .coordinates <- model.frame(formula, data = g$data[[vName]]$coordinates)
    .coordinates <- as.matrix(.coordinates)

    # Check whether to use residuals or raw data:
    if(useResidual && attr(g$data[[vName]], 'hasDrift')){
      .functions <- as.matrix(g$drift[[vName]]$Residuals)
    } else {
      .functions <- as.matrix(g$data[[vName]]$functions)
    }

    # Check if data is incomplete:
    if(!attr(g$data[[1]], 'isScalar')){
      if(any(is.na(.functions))){
        partialData = TRUE
        if(comments){
          print("Partial data was detected! Trace variograms will be computed with Gromenko's formula!")
        }
      } else {
        partialData = FALSE
        if(comments){
          print("Functions are complete! Trapezoidal method will be used for trace variograms!")
        }
      }
    } else {
      partialData = FALSE
      if(comments){
        print("Scalar data detected. Conventional (scalar) variograms will be computed")
      }
    }

    .CppVariogram <- empTraceVariogram(.coordinates,
                                       .functions,
                                       as.matrix(Lags),
                                       as.matrix(directions),
                                       LagTolerance,
                                       ArgStep,
                                       angTolerance,
                                       partialData)

    singleVariogram <- gstatify(.CppVariogram, label = vName, dir.label = dirNames)

    if(i == 1) {
      OutputStructure <- singleVariogram
    } else {
      OutputStructure <- rbind.data.frame(singleVariogram, OutputStructure)
    }

    if(crossMode != "NA" && i > 1) { # NA is to skip cross variogram computation all together.

      for(j in 1:(i-1)) {

        aux_vName <- names(g$data)[j] # name of the auxiliary data.

        if(crossMode == "conventional"){  # conventional cross-variogram

          # Apply formula to coordinates:
          .auxCoordinates  <- model.frame(formula, data = g$data[[aux_vName]]$coordinates)
          .auxCoordinates  <- as.matrix(.auxCoordinates)

          isoCoord <- isoCoordinates(.auxCoordinates,
                                     .coordinates)

          isoCoord <- isoCoord[!is.na(isoCoord[,2]), ]

          if(nrow(isoCoord) == length(sum(is.na(isoCoord[,2])))){
            stop(paste(
              'Error encountered while computing conventional cross-variogram ',
              'for variables: ', aux_vName, ' and ', vName, '.',
              'The data is heterotopic! Please consider using pseudo-cross-variogram! ',
              'Exiting!', sep('')
            )
            )
          }

          # Check whether to use residuals or raw data:
          if(useResidual && attr(g$data[[aux_vName]], 'hasDrift')){
            .auxFunctions <- as.matrix(g$drift[[aux_vName]]$Residuals)
          } else {
            .auxFunctions <- as.matrix(g$data[[aux_vName]]$functions)
          }

          # This segement is adjusting the input to Cpp function depending on they type of data
          # we are forecasting (scalars or functions)
          .isCoord <- .auxCoordinates[isoCoord[, 1], ]

          if(ncol(.functions) > 1) { # functions
            .auxF  <- .auxFunctions[, isoCoord[, 1]]
            .primF <- .functions[, isoCoord[, 2]]
          } else { # scalars!
            .auxF  <- as.matrix(.auxFunctions[isoCoord[, 1],])
            .primF <- as.matrix(.functions[isoCoord[, 2], ])
          }

          .CppxVario <- empTraceCrossVariogram(.isCoord,
                                               .auxF,
                                               .primF,
                                               as.matrix(Lags),
                                               as.matrix(directions),
                                               LagTolerance,
                                               ArgStep,
                                               angTolerance)

        } else { # pseudo cross-variogram

          # Apply formula to coordinates:
          .auxCoordinates  <- model.frame(formula, data = g$data[[aux_vName]]$coordinates)
          .auxCoordinates  <- as.matrix(.auxCoordinates)

          # Check whether to use residuals or raw data:
          if(useResidual && attr(g$data[[aux_vName]], 'hasDrift')){
            .auxFunctions <- as.matrix(g$drift[[aux_vName]]$Residuals)
          } else {
            .auxFunctions <- as.matrix(g$data[[aux_vName]]$functions)
          }

          if(ncol(.functions) > 1) { # functions
            .auxF  <- .auxFunctions
            .primF <- .functions
          } else { # scalars!
            .auxF  <- as.matrix(.auxFunctions[,1])
            .primF <- as.matrix(.functions[,1])
          }

          .CppxVario <- empTracePseudoCrossVariogram(.coordinates,
                                                     .primF,
                                                     .auxCoordinates,
                                                     .auxF,
                                                     as.matrix(Lags),
                                                     as.matrix(directions),
                                                     LagTolerance,
                                                     ArgStep,
                                                     angTolerance)


        }

        xVario <- gstatify(.CppxVario, label = paste(aux_vName, vName, sep = "."), dir.label = dirNames)

        # Append to the output structure:
        OutputStructure <- rbind.data.frame(xVario, OutputStructure)

      }

    }

  }

  OutputStructure <- OutputStructure[!apply(OutputStructure, 1, function(x) any(is.na(x))), ]

  g$variogram     <- OutputStructure
  g$KrigFormula   <- formula

  return(g)

}





#' @title A function that fits variograms to empirical estimates
#'
#' @param g               - An \code{fstat} structure.
#' @param model           - Variogram model of class \code{'vgm'}.
#' @param dir             - Direction you wish to fit. Useful for fine tunning.
#' @param correctDiagonal - Default = 1. Corrects the diagonal of the LMC sill matrix.
#' @param fitSills        - This is a convenience access parameter for fit.variogram
#' @param fitRanges       - This is a convenience access parameter for fit.variogram
#' @param forceNugget     - Whether to force the nugget to user specified value or to fit it (Boolean).
#' @param posNugget       - Whether to enforce the positivity of the nugget (Boolean). Not to be used
#'                          with conventional cross variogram computation!
#'
#' @details It returns the \code{fstat} structure with an addtional slot \code{models}.
#'          This function will cycle through all available directions and fit either a single-variate
#'          variogram or an LMC structure (depending on the number of id's).
#'          variable \code{dir} is available for the purpose of fine tunning, i.e. different variogram for
#'          a different direction.
#'          In the backend, this code relies on the \code{fit.variogram} routine from the \code{gstat} package.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu})
#'
#' @export
#'
fitVariograms <- function(g = NULL,
                          model = NULL,
                          dir = NULL,
                          correctDiagonal = 1,
                          fitSills = TRUE,
                          fitRanges = TRUE,
                          forceNugget = FALSE,
                          posNugget = FALSE){

  if(class(g) != 'fstat' | is.null(g)) stop('This function requires objects of class "fstat". Exiting!')
  if(is.null(g$variogram)) stop('Empirical variogram structure was not found. Exiting!')

  if(is.null(dir)){
    .directions <- unique(g$variogram$dir.hor)
  } else {
    if(is.na(match(.dir,unique(g$variogram$dir.hor)))) {
      stop('The provided direction was not found in the fstat structure. Exiting!')
    }else {
      .directions <- dir
    }
  }

  # it fits one variogram for each direction.
  for(.dir in .directions){

    # subset with the current direction:
    .currentSet <- subset(g$variogram, dir.hor == .dir)

    if(length(unique(.currentSet$id)) == 1) {

      if(forceNugget){
        tempModel <- fit.variogram(.currentSet, model, warn.if.neg = TRUE, fit.sills = fitSills, fit.ranges = fitRanges)
        tempModel$psill <- c(model$psill[1], sum(tempModel$psill))

        g$model[[as.character(.dir)]][[as.character(unique(.currentSet$id))]] <- fit.variogram(.currentSet, tempModel, warn.if.neg = TRUE,
                                                                                               fit.sills = FALSE, fit.ranges = TRUE)

      } else{
      g$model[[as.character(.dir)]][[as.character(unique(.currentSet$id))]] <- fit.variogram(.currentSet, model, warn.if.neg = TRUE,
                                                                                             fit.sills = fitSills, fit.ranges = fitRanges)
      }

    } else {

      if(forceNugget){ # Individual fitting. This is here because it helps us work with compact scripts.

        # loop through every id and fit one variogram:
        for(.id in unique(.currentSet$id)){

          tempModel <- fit.variogram(subset(.currentSet, id == .id), model, warn.if.neg = TRUE, fit.sills = fitSills, fit.ranges = fitRanges)
          tempModel$psill <- c(model$psill[1], sum(tempModel$psill))

          g$model[[as.character(.dir)]][[as.character(.id)]] <- fit.variogram(subset(.currentSet, id == .id), tempModel, warn.if.neg = TRUE,
                                                                                    fit.sills = FALSE, fit.ranges = fitRanges)
        }

      } else { # lmc:
        g$model[[as.character(.dir)]] <- finterpfitLMC(.currentSet, names(g$data), model, fit.ranges = fitRanges,
                                                       correct.diagonal = correctDiagonal, fNugget = posNugget)
      }

    }

  }
  return(g)
}

#' @title Plots available variogram structure (model and empirical)
#'
#' @param g     - An \code{fstat} structure
#' @param theme - \code{ggplot} theme (optional)
#' @param ggReturn - Boolean, specifying whether to return a ggplot structure or to plot.
#' @param dir      - plot only one specific direction
#' @param npSize   - change the size of empirical values based on the number of points (boolean).
#'
#' @details The function will plot empirical varigram by default. If a variogram model is available
#'          then it will plot the model alongside the empirical values.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu}) and Alessandra Menafoglio (\email{alessandra.menafoglio@polimi.it})
#'
#' @export
#'
plotVariogram <- function(g, theme = NULL, ggReturn = FALSE, dir = NULL, npSize = FALSE){

  if(class(g) != 'fstat' | is.null(g)) stop('This function requires objects of class "fstat". Exiting!')
  if(is.null(g$variogram)) stop('Empirical variogram structure was not found. Exiting!')

  # This part determines whetehr to plot variogram line or not. If one exists in the structure it will plot!
    if(!is.null(g$model)){

      .directions <- unique(g$variogram$dir.hor)

      for(i in 1:length(.directions)){
          .dir <-  .directions[i]
          .currentSet <- subset(g$variogram, dir.hor == .dir)
          .vars <- unique(.currentSet$id)
          .currentSet$fit = NA

          for(.var in .vars){
            .currentSet$fit[.currentSet$id == .var] =
              variogramLine( g$model[[as.character(.dir)]][[.var]],
                             dist_vector = .currentSet$dist[.currentSet$id == .var]
              )$gamma
          }

          if(i==1) {
            .ggData = .currentSet
          } else {
            .ggData = rbind(.ggData, .currentSet)
          }

      }

      .ggLine <- geom_line(aes(x = dist,y = fit, group = dir.hor, colour = dir.hor))

    } else {
      .ggData <- g$variogram
      .ggLine <- NULL
    }


  .ggData$dir.hor <- factor(.ggData$dir.hor)

  # Sets flags for plotting points of variable size depending on the number of pairs.
  if(npSize){
    npSize = "np"
  } else {
    npSize = NULL
  }

  # Subsets the data to plot only one direction.
  if(!is.null(dir)){
   .ggData <- subset(.ggData, dir.hor == dir)
  }

  .p <- ggplot(.ggData) + geom_point(aes_string(x = "dist", y = "gamma", colour = "dir.hor", group = "dir.hor", size = npSize)) +
    .ggLine + facet_wrap(~id, scales = "free_y") + theme

  # Return ggPlot structure or not? This is useful for external visual tunning (axis width, font etc)
  if(ggReturn){
    return(.p)
  } else {
    print(.p)
  }

}

computeCovariance <- function(g = NULL, type = NULL, bigData = FALSE){

  if(class(g) != 'fstat' | is.null(g)) stop('This function requires objects of class "fstat". Exiting!')

  type = match.arg(type, c('sum', 'product', 'productSum', 'omnidirectional'))

}

#' @title A function that converts Cpp output to gstat variogram type class:
#'
#' @param M - A matrix produced by one of the Cpp variogram computing routines (i.e. \code{empTraceVariogram});
#' @param label - A label you wish to assign as \code{id} in variogram structure (The Cpp code does not return any variable id);
#' @param dir.label - A vector of labels for directional variogram. The Cpp code returns [0,1,2,...].
#'
#' @return A data frame of class: \code{gstatVariogram}.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu})
#'
#' @details Not to be called directly! Serves as a suport function only!
#'          We expose it anyways because it might be useful in research.
#'
#' @export
#'
gstatify <- function(M = NULL, label = NULL, dir.label = NULL){

  if(is.null(M))     stop('Matrix not provided! Exiting.')
  if(is.null(label)) stop('Label not provided! Exiting.')

  if(is.null(dir.label)){

    if(length(unique(M[, 4])) == 1){

      dir.hor = unique(M[, 4])

    } else {

      dir.hor = M[, 4]

    }

  } else {

    if(length(unique(M[, 4])) == length(dir.label)) {

      uniqueList <- sort(unique(unique(M[, 4])))

      dir.hor = M[, 4]

      for(i in length(uniqueList):1){

        dir.hor = gsub(uniqueList[i], dir.label[i], dir.hor, fixed = TRUE)

      }

    } else {

      print('Error in gstatify, the number of labels is not maching the number of directions. Exiting!')

    }

  }

  Out <- data.frame(np      = M[, 2],
                    dist    = M[, 1],
                    gamma   = M[, 3],
                    dir.hor = dir.hor,
                    dir.ver = 0,
                    id      = label,
                    dist.av = M[, 5])

  class(Out) <- c("gstatVariogram", "data.frame")

  return(Out)
}


#' @title A modified version of the 'fit.lmc' function from the gstat package
#'
#' @details Not to be called directly. This is an internal function!
#'
#' @author Modified by: Ogy Grujic (\email{ogyg@stanford.edu})
#'
#'
finterpfitLMC <- function (.v, .n, .model, fit.ranges = FALSE, fit.lmc = !fit.ranges,
                           correct.diagonal = 1, fNugget = FALSE, ...)
{
  posdef = function(X) {
    q = eigen(X)
    d = q$values
    d[d < 0] = 0
    q$vectors %*% diag(d, nrow = length(d)) %*% t(q$vectors)
  }

  .g = NULL

  for (i in 1:length(.n)) {
    for (j in i:length(.n)) {
      name = ifelse(i == j, .n[i], cross.name(.n[i], .n[j]))
      .x = .v[.v$id == name, ]
      if (nrow(.x) == 0)
        stop(paste("gstatVariogram", name, "not present"))

      if(fNugget){
        .temp <- fit.variogram(.x, .model, fit.ranges = fit.ranges, ...)
        if(any(.temp$psill < 0)){
          .totalSill = sum(.temp$psill)
          .temp$psill[.temp$psill<0] = 0
          .temp$psill[.temp$psill > 0] = .totalSill
          .g$model[[name]] = fit.variogram(.x, .temp, fit.ranges = fit.ranges, fit.sills=FALSE, ...)
        } else {
          .g$model[[name]] = .temp
        }
      } else {
      .g$model[[name]] = fit.variogram(.x, .model, fit.ranges = fit.ranges,
                                      ...)
      }
    }
  }
  if (fit.lmc) {
    .m = .g$model[[.n[1]]]
    for (k in 1:nrow(.m)) {
      psill = matrix(NA, nrow = length(.n), ncol = length(.n))
      for (i in 1:length(.n)) {
        for (j in i:length(.n)) {
          name = ifelse(i == j, .n[i], cross.name(.n[i],
                                                 .n[j]))
          psill[i, j] = psill[j, i] = .g$model[[name]][k, "psill"]
        }
      }
      psill =  posdef(psill) #as.matrix(nearPD(psill, keepDiag = TRUE)$mat)
      diag(psill) = diag(psill) * correct.diagonal
      for (i in 1:length(.n)) {
        for (j in i:length(.n)) {
          name = ifelse(i == j, .n[i], cross.name(.n[i],
                                                 .n[j]))
          .g$model[[name]][k, "psill"] = psill[i, j]
        }
      }
    }
  }
  return(.g$model)
}


#' @title Computes covariance matrices from available variogram models within the fstat structure
#'
#' @param .g   -  An 'fstat' structure
#' @param type -  The type of covariance to compute ('omni', 'separable', 'sum', 'product', 'productSum')
#'
#' @details It will compute a covariance matrix for each variable in the fstat structure.
#'          Note: 'separable, 'product', 'sum', productSum' are experimental in the current version.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu})
#'
#' @export
addCovariance <- function(.g = NULL, type = 'omni'){

  if(class(.g) != 'fstat' | is.null(.g)) stop('This function requires objects of class "fstat". Exiting!')

  type <- match.arg(type, c("omni", "sum", "product", "productSum"))

  if(is.null(.g$model)) stop('Variogram model was not found in the fstat structure. Exiting!')

  # initialize covariance structure.
   .g$covariance <- NULL

  # HELPER FUNCTIONS: ##########################################################
  .covCompute <- function(.direction, .formula){

    .names <- names(.g$model[[.direction]])

    .covReturn <- vector('list', length(.names))

    names(.covReturn) <- .names

    for(.name in .names) {

      .vars <- unlist(strsplit(.name, "[.]"))

      if(length(.vars) == 1){ # its trace variogram structure:

        # Compute Distance:
        .C <- model.frame(.formula, data = .g$data[[.name]]$coordinates)
        .Dist <- as.matrix(dist(.C))
        # .varioLine <- variogramLine(.g$model[[.direction]][[.name]], dist_vector = .Dist)
        .sill <- sum(.g$model[[.direction]][[.name]]$psill)

        # .covReturn[[.name]] <- .sill - .varioLine
        .covReturn[[.name]] <- variogramLine(.g$model[[.direction]][[.name]],
                                             dist_vector = .Dist, covariance = TRUE)
        attr(.covReturn[[.name]], 'sill') <- .sill

      } else { # Its a cross variogram structure:

        .C1 <- model.frame(.formula, data = .g$data[[.vars[1]]]$coordinates)
        .C2 <- model.frame(.formula, data = .g$data[[.vars[2]]]$coordinates)
        .C  <- rbind(.C1, .C2)

        .n1   <- nrow(.g$data[[.vars[1]]]$coordinates)
        .n2   <- nrow(.g$data[[.vars[2]]]$coordinates)

        .Dist <- as.matrix(dist(.C))[(1:.n1), -(1:.n1)]

        # .varioLine <- variogramLine(.g$model[[.direction]][[.name]], dist_vector = .Dist)
        .sill <- sum(.g$model[[.direction]][[.name]]$psill)

        # .covReturn[[.name]] <- .sill - .varioLine
        .covReturn[[.name]] <- variogramLine(.g$model[[.direction]][[.name]],
                                             dist_vector = .Dist, covariance = TRUE)
        attr(.covReturn[[.name]], 'sill') <- .sill
      }

    }

    return(.covReturn)
  }

  .covAggregate <- function(.operator){

    .covReturn <- .covStruct[[1]]

    .nElements <- length(.covStruct)
    .nSubElem  <- length(.covReturn)

    for(i in 1:.nElements){
      for(j in 1:.nSubElem){
        .command <- paste(".covReturn[[j]] <- .covReturn[[j]]", .operator, ".covStruct[[i]][[j]]", sep = " ")
        .command <- parse(text = .command)
        eval(.command)
      }
    }

    return(.covReturn)

  }

  ##############################################################################

  if(type == "omni"){ #omni is always present by design.

    .g$covariance[['omni']] <- .covCompute('omni', .g$KrigFormula)

    attr(.g$covariance, 'covType') <- 'omni'

  } else { # all other directions:

    .directions <- names(.g$model)
    .directions <- .directions[-which(.directions == 'omni')]

    .covStruct <- vector('list', length(.directions))
    names(.covStruct) <- .directions

    for(.dir in .directions){
      .covStruct[[.dir]] <- .covCompute(.dir, paste("~", .dir, sep = ""))
    }

    .operation <- switch(type,
                        sum     = "+",
                        product = "*")

    .aggregate    <- .covAggregate(.operation)

    .g$covariance[[type]] <- .aggregate

    attr(.g$covariance, 'covType') <- type

  }

  return(.g)

}

#' @title The function estimates functional or scalar drift.
#'
#' @param .formula - formula specifying which variables to use for drift;
#' @param .g       - an fstat structure containing all data.
#' @param .type    - drift type (OLS, GLS, MEAN)
#' @param .regCoef - Regularization coefficient (just like in ridge regression)
#'
#' @details It adds a "drift" slot to the fstat structure (.g) with the following components:
#'          Betas         - Coefficients assigned to each variable specified by the formula;
#'          Design matrix - A design matrix for each variable (useful for UK);
#'          Residuals     - Functional or scalar residuals;
#'
#'          This code does not implement functional regression (Ramsay and Silverman, 2005), instead it fits one regression
#'          for each time step (piece-wise). The actual functional regression will be implemented in future version(s).
#'
#'          for .type=MEAN it just subtracts the mean from the functional data.
#'
#'@author Ogy Grujic (\email{ogyg@stanford.edu}) and Alessandra Menafoglio (\email{alessandra.menafoglio@polimi.it})
#'
#'@export
#'
estimateDrift <- function(.formula, .g, .type = "OLS", Intercept = TRUE, regCoef = 0){

  if(class(.g) != 'fstat' | is.null(.g)) stop('This function requires objects of class "fstat". Exiting!')

  if(is.null(.formula)) stop('Formula was not provided. Exiting!')

  if(is.character(.formula)){
    .formula.list <- vector("list", length(.g$data))
    names(.formula.list) <- names(.g$data)
    .formula.list <- lapply(.formula.list, function(x) return(.formula))
    .formula <- .formula.list
  } else if(is.list(.formula)){
    if(length(.formula) != length(.g$data)){
      stop("The length of formula list is not the same as the number of outputs. Exiting!")
    } else {
      names(.formula) <- names(.g$data)
    }
  }

  .type <- match.arg(.type, c("OLS", "GLS", "MEAN")) # only three methods have been implemented so far.
  # Plan is to implement the method by Ramsay and Silverman with constrains on drift parameters (betas)
  # Namely betas should sum to one for every time step.

  # If its GLS check if covariance structure exists!

  ### Helper Functions: ########
  .olsDrift <- function(.X, .Ys){

    .X  <- as.matrix(.X)
    .Ys <- as.matrix(.Ys)

    .XtX       <- t(.X)%*%.X

    .XtX_inv   <- solve(regCoef*diag(nrow(.XtX)) + .XtX)   # Potential problem (consider pseudo inverse!)

    .Hat       <- .X %*% .XtX_inv %*% t(.X)

    if(ncol(.Ys) == 1) .Ys <- t(.Ys)

    .Ys_hat    <- t(.Hat %*% t(.Ys))

    .Residuals <- .Ys - .Ys_hat

    if(nrow(.Ys) == 1)  .Residuals <- matrix(.Residuals, ncol = 1)

    .Betas     <- t(.XtX_inv %*% t(.X) %*% t(.Ys))

    return(list(Betas       = .Betas,
                fHat        = .Ys_hat,
                Residuals   = .Residuals
    )
    )
  }
  .glsDrift <- function(.X, .Ys, .sigma){

    .X         <- as.matrix(.X)
    .Ys        <- as.matrix(.Ys)
    .sigma     <- as.matrix(.sigma)

    .L         <- t(chol(.sigma))

    Linv_X     <- solve(.L) %*% .X

    if(ncol(.Ys) == 1){
      # Linv_Y     <- solve(.L) %*% .Ys
      .Betas <- solve(t(.X)%*%solve(.sigma)%*%.X)%*%t(.X)%*%solve(.sigma)%*%.Ys
      .Ys_hat    <- .X %*% .Betas
    } else {
      Linv_Y     <- solve(.L) %*% t(.Ys)
      Xt_Sigma_X <- t(Linv_X) %*% Linv_X

      .Betas     <- solve(Xt_Sigma_X + regCoef*diag(nrow(Xt_Sigma_X))) %*% t(Linv_X) %*% Linv_Y
      .Ys_hat    <- t(.X %*% .Betas)
    }

    .Residuals <- .Ys - .Ys_hat

    return(list(Betas       = t(.Betas),
                fHat        = .Ys_hat,
                Residuals   = .Residuals
    )
    )
  }
  .mean     <- function(.Ys){

    .Ys_hat    <- apply(.Ys, 1, mean)

    .Residuals <- .Ys - .Ys_hat

    return(list(Betas       = NA,
                fHat        = .Ys_hat,
                Residuals   = .Residuals
    )
    )
  }

  ###### MAIN #######
  .g$drift        <- vector('list', length(names(.g$data)))
  names(.g$drift) <- names(.g$data)

  for(.name in names(.g$data)){

    # Compute design matrix:
    if(.type != "MEAN") {
      .ModelFrame <- model.frame(.formula[[.name]], data = .g$data[[.name]]$coordinates)

      # add intercept: This needs to be expanded so it works with the formula!
      if(Intercept){
        .ModelFrame <- cbind(1L, .ModelFrame)
      }
    } else {
        .ModelFrame <- rep(1, nrow(.g$data[[.name]]$coordinates))
    }

    # Compute Drift: Needs modification for covariance in the case when its something else but "omni"
    .driftEstimate <- switch(.type,
                             OLS  = .olsDrift(.ModelFrame, .g$data[[.name]]$functions),
                             GLS  = .glsDrift(.ModelFrame, .g$data[[.name]]$functions, .g$covariance[[1]][[.name]]),
                             MEAN = .mean(.g$data[[.name]]$functions)
                             )

    .g$drift[[.name]]                  <- .driftEstimate
    .g$drift[[.name]]$DesignMatrix     <- .ModelFrame

    attr(.g$drift, 'type')             <- .type
    attr(.g$drift, 'formula')          <- .formula

    # Adjust global attribute:
    attr(.g$data[[.name]], 'hasDrift') <- TRUE
    attr(.g$drift[[.name]], 'hasIntercept') <- Intercept
  }
  return(.g)
}

#' @title A general predict function for fstat structure
#'
#' @param .g         - an fstat structure
#' @param .what      - variable you wish to predict (must match one of the already present names)
#' @param .type      - kriging type any of the following: 'OK', 'UK', 'OcoK', and 'UcoK'. Applies to both scalar and functional cases.
#' @param ErrorCheck - whether to check for errors on each segment and print info.
#' @param algIndependent - Boolean, specifying whether drifts in cokriging are independent (default = FALSE)
#'
#' @details In all methods we use covariances (no variograms) by subtracting the variogram from the sill. The type of covariance structure is
#'          extracted from the fstat structure.
#'
#'          We don't provide functionality for simple kriging and simple co-kriging since doing them is just silly.
#'
#'          UcoK can have three types of drifts. Algebraically independent, algebraically dependent and mixed (see Chiles and Delfiner, 1998)
#'          This code implements the first two cases (algIndependent flag). Mixed drift will be implemented in one of the future versions of this code.
#'
#' @author Ogy Grujic \email{ogyg@stanford.edu}
#'
#' @export
#'
predictFstat <- function(.g, .newCoordinates = NULL, .what = NULL, .type = 'OK', ErrorCheck = FALSE, algIndependent = FALSE){

  if(class(.g) != 'fstat' | is.null(.g)) stop('This function requires objects of class "fstat". Exiting!')

  if(is.null(.newCoordinates)) stop('New coordinates have to be specifed! Exiting!')

  if(is.null(.what)) stop('Please specify which variable to predict! (.what = ?)')

  .what <- match.arg(.what, names(.g$data))

  .type <- match.arg(.type, c("OK", "UK", "OcoK", "UcoK"))

  if((.type == "UcoK" | .type == "UK") && !attr(.g$data[[.what]], 'hasDrift')) {
    stop('You need to pre-estimate the drift in order to use Universal (co)Kriging. Either estimate the drift or try OK or OcoK')
  }

  # HELPER FUNCTIONS: ##########################################################
  # computes covariance for a given direction, variable index, and covariance type (cross/auto).
  .covCompute <- function(.DIR, .INDEX, .COVNAME, .FORMULA) {

    .d          <- NULL
    .newCoord   <- model.frame(.FORMULA, data = .newCoordinates)
    .trainCoord <- model.frame(.FORMULA, data = .g$data[[.INDEX]]$coordinates)

    .d          <- as.matrix(dist(rbind(.trainCoord, .newCoord)))
    .d          <- .d[1:nrow(.trainCoord),-(1:nrow(.trainCoord))]

    .d          <- as.matrix(.d)

    # .cSill      <- sum(.g$model[[.DIR]][[.COVNAME]]$psill)

    # return(.cSill - variogramLine(.g$model[[.DIR]][[.COVNAME]], dist_vector = .d))
    return(variogramLine(.g$model[[.DIR]][[.COVNAME]], dist_vector = .d, covariance=TRUE))
  }

  # Aggregates unidirectional covariances:
  .covAggregate <- function(.covStruct, .operator){

    .covReturn <- .covStruct[[1]]

    .nElements <- length(.covStruct)
    .nSubElem  <- length(.covReturn)

    for(i in 1:.nElements){
      for(j in 1:.nSubElem){
        .command <- paste(".covReturn[[j]] <- .covReturn[[j]]", .operator, ".covStruct[[i]][[j]]", sep = " ")
        .command <- parse(text = .command)
        eval(.command)
      }
    }

    return(.covReturn)

  }

  # Boolean operator that checks if matrix is positive definite.
  is.posdef = function(M) {
    !any(eigen(M)$values < 0)
  }

  # This function creates an F matrix for universal co-kriging!
  .createF <- function(){

    nD <- ncol(.g$drift[[1]]$DesignMatrix)
    nP <- length(.g$data)

    .Fout <-  NULL

    Pts <- unlist(lapply(.g$drift, function(x) nrow(x$DesignMatrix)))

    nPts <- sum(Pts)

    # algebraically dependent drifts:
    if(!algIndependent){
      for(i in 1:nP){
        fMiddle <- .g$drift[[i]]$DesignMatrix # Always the same!
        .Fout <- rbind(.Fout, fMiddle)
      }
    } else { # algebraic independence

      .leadMatrix <- NULL
      .tailMatrix <- matrix(0, ncol = nD,
                            nrow = nPts)

      for(i in 1:nP){

        .designMatrix <- as.matrix(.g$drift[[i]]$DesignMatrix)
        if(i == 1){
          .leadMatrix = NULL
        } else {
          .leadMatrix   <- matrix(0, nrow = sum(Pts[1:(i-1)]), ncol=ncol(.designMatrix))
        }

        if(i == nP){
          .tailMatrix = NULL
        } else {
          .tailMatrix   <- matrix(0, nrow = sum(Pts[(i+1):length(Pts)]), ncol=ncol(.designMatrix))
        }

        #
        # # augment tail:
        # .tailMatrix <- .tailMatrix[-(1:nrow(.designMatrix)),]
        #
        .currentMatrix <- rbind(.leadMatrix, .designMatrix, .tailMatrix)
        #
        # # augment lead:
        # .leadMatrix <- rbind(.leadMatrix, .designMatrix * 0)

        # append:
        .Fout <- cbind(.Fout, .currentMatrix)

      }
    }

    return(as.matrix(.Fout))

  }

  if(.type == "OcoK" | .type == "UcoK"){     # COKRIGING

    k     = which(names(.g$data) == .what)   # Get the index of the forecasted variable.

    # PREPARE SIMPLE CO-KRIGING SYSTEM OF EQUAITONS: START *********************
    lhs   = NULL  # initialize covariance matrix:
    rhs   = NULL  # initialize right hand side:

    for(i in 1:length(names(.g$data))){

      # LEFT HAND SIDE COMPONENT: START ----------------------------------------
      H = NULL    # initialize horizontal matrix:

      for(j in 1:length(names(.g$data))){
        if(i==j){
          .covName = names(.g$data)[j]
          h = .g$covariance$omni[[.covName]]
        } else if(i<j) {
          .covName = paste(names(.g$data)[i], names(.g$data)[j], sep = ".")
          h = .g$covariance$omni[[.covName]]
        } else {
          .covName = paste(names(.g$data)[j], names(.g$data)[i], sep = ".")
          h = t(.g$covariance$omni[[.covName]])
        }
        H = cbind(H,h)  # add to horizontal matrix
      }

      lhs <- rbind(lhs, H) # Add to simple co-Kriging matrix (left hand side (lhs))

      # LEFT HAND SIDE COMPONENT: END ------------------------------------------


      # RIGHT HAND SIDE COMPONENT: START ---------------------------------------

      # Determine variogram type (based on the output)
      if(i == k){
        .covName = names(.g$data)[k]
      } else if(i < k) {
        .covName = paste(names(.g$data)[i], names(.g$data)[k], sep = ".")
      } else {
        .covName = paste(names(.g$data)[k], names(.g$data)[i], sep = ".")
      }

      .c <- NULL   # initialize rhs component.

      # Compute covariance:
      if(attr(.g$covariance, 'covType') == "omni"){     # omni is always present by design.

        .c     <-  .covCompute("omni", i, .covName, .g$KrigFormula)

      } else { # all other types (SUM, PRODUCT)

        .directions <- names(.g$model)  # this matches kriging formula!
        .directions <- .directions[-which(.directions == 'omni')]

        # Initialize components of the sum/product covariance:
        .covStruct <- vector('list', length(.directions))
        names(.covStruct) <- .directions

        for(.dir in .directions){
          .covStruct[[.dir]] <- .covCompute(.dir, i, .covName, paste("~", .dir, sep = ""))
        }

        # Assemble covariance:
        .aggregate    <- .covAggregate(.covStruct, switch(attr(.g$covariance, 'covType'),
                                                          sum     = "+",
                                                          product = "*"
        )
        )
        .c <- .aggregate
      }

      rhs <- rbind(rhs, .c) # Append to the right hand side.

      # RIGHT HAND SIDE COMPONENT: END -----------------------------------------

    } # PREPARE  A SIMPLE CO-KRIGING SYSTEM OF EQUAITONS: END ******************

    # check positive definiteness:
    if(ErrorCheck) {
      print(paste("Simple co-kriging matrix ", ifelse(is.posdef(lhs), "is", "is NOT"), " positive definite!"))
    }

    # Matrix augmentation: adds unbiasedness constraints to covariance matrix
    # hence creating an ordinary co-kriging matrix.
    .repNames  <- unlist(lapply(1:length(.g$data),
                                function(i)
                                  rep(names(.g$data)[i], nrow(.g$data[[i]]$coordinates))
    )
    )
    .augVector1                     <- rep(0, length(.repNames))
    .augVector2                     <- .augVector1
    .augVector1[.repNames == .what] <- 1
    .augVector2[.repNames != .what] <- 1
    .augMatrix <- cbind(.augVector1, .augVector2)


    # Ordinary co-Kriging matrix:

    # UPGRADE**
    if(.type == "UcoK"){ # More augmentation for Universal co-Kriging. -----------------

      # NOTE: this segment assumes that drift exists, function header checks that it does!

      # It should create either algebraically dependent or independent drift.

      .Fmatrix <- .createF()

      .X <- model.frame(attr(.g$drift, 'formula')[.what][[1]], data = .newCoordinates)

      if(attr(.g$drift[[.what]], 'hasIntercept')) {
        .X <- cbind(1, .X)
      }

      # Augmentation for algebraically INDEPENDENT drifts.
      if(algIndependent){

        .FrhsAug <- NULL
        Pts <- unlist(lapply(.g$drift, function(x) nrow(x$DesignMatrix)))

        for(k in 1:length(.g$data)){
          if(.what == names(.g$data)[[k]]){
            if(k == 1){
              .FrhsAug <- .X
            } else {
              .FrhsAug <- cbind(.FrhsAug, .X)
            }
          } else {
            if(is.null(.FrhsAug)){
              .FrhsAug <- matrix(0, nrow = nrow(.X), ncol = ncol(.g$drift[[k]]$DesignMatrix))
            } else {
              .FrhsAug <- cbind(.FrhsAug, matrix(0, nrow = nrow(.X), ncol = ncol(.g$drift[[k]]$DesignMatrix)))
            }
          }
        }
        .X <- .FrhsAug
      }

      .solution <- krigSolver(lhs, rhs, .Fmatrix, .X)

      # Compute forecast: note we are using original curves for this! Modify for scalars!
      if(attr(.g$data[[1]], 'isScalar')){
        .fs <- as.matrix(unlist(lapply(.g$data, function(x) x$functions)))
      } else {
        .fs <- as.matrix(Reduce(cbind, lapply(.g$data, function(x) x$functions)))
      }

      # Check for scalars and adjust appropriatelly:
      if(ncol(.fs) > 1){
        .forecast <-  .fs %*% .solution$weights
      } else {
        .forecast <-  t(t(.fs) %*% .solution$weights)
      }


      # .forecast <- .fs %*% .solution$weights


    } else {             # ORDINARY CO-KRIGING SYSTEM OF EQUATIONS --------------------

      .augMatrix2 <- matrix(1, nrow = nrow(.newCoordinates), ncol = 2)
      .augMatrix2[,2] <- 0

      .solution <- krigSolver(lhs, rhs, .augMatrix, .augMatrix2)

      if(attr(.g$data[[1]], 'isScalar')){
        .fs <- as.matrix(unlist(lapply(.g$data, function(x) x$functions)))
      } else {
        .fs <- as.matrix(Reduce(cbind, lapply(.g$data, function(x) x$functions)))
      }

      # Check for scalars and adjust appropriatelly:
      if(ncol(.fs) > 1){
        .forecast <-  .fs %*% .solution$weights
      } else {
        .forecast <-  t(t(.fs) %*% .solution$weights)
      }

    }

    # Return Co-kriging structure (C and c are there for debugging purposes!):
    return(list(C = lhs, c = rhs, Variance = .solution$variance, Forecast = .forecast, Weights = .solution$weights))

  } else {            # KRIGING methods: ***************************************

    # FORM SIMPLE KRIGING SYSTEM OF EQUATIONS: START ---------------------------
    # get appropriate matrix:
    lhs <- .g$covariance[[attr(.g$covariance, 'covType')]][[.what]]

    k   <- which(names(.g$data) == .what)   # Get the index of the forecasted variable.

    # Get the right hand side:
    if(attr(.g$covariance, 'covType') == "omni"){     # omni is always present by design.

      rhs     <-  .covCompute("omni", k, .what, .g$KrigFormula)

    } else { # all other types (SUM, PRODUCT)

      .directions <- names(.g$model)  # this matches kriging formula!
      .directions <- .directions[-which(.directions == 'omni')]

      # Initialize components of the sum/product covariance:
      .covStruct <- vector('list', length(.directions))
      names(.covStruct) <- .directions

      for(.dir in .directions){
        .covStruct[[.dir]] <- .covCompute(.dir, k, .what, paste("~", .dir, sep = ""))
      }
      # Assemble covariance:
      .aggregate    <- .covAggregate(.covStruct, switch(attr(.g$covariance, 'covType'),
                                                        sum     = "+",
                                                        product = "*"
      ))
      rhs <- .aggregate
    }

    # FORM SIMPLE KRIGING SYSTEM OF EQUATIONS: END -----------------------------

    if(.type == "OK"){ # Ordinary kriging --------------------------------------

      .DM <- matrix(1, ncol = 1, nrow = nrow(lhs))
      .X  <- matrix(1, ncol = 1, nrow = ncol(rhs))

      .solution <- krigSolver(lhs, rhs, .DM, .X)

      .fs <- as.matrix(.g$data[[.what]]$functions)
      .forecast <-  .fs %*% .solution$weights

      return(list(C = lhs, c = rhs, Variance = .solution$variance, Forecast = .forecast, Weights = .solution$weights))

    } else {            # Universal Kriging --------------------------------------

      .X <- model.frame(attr(.g$drift, 'formula')[[.what]], data = .newCoordinates)

      if(attr(.g$drift[[.what]], 'hasIntercept')) {
        .X <- cbind(1, .X)
      }

      .DM <- as.matrix(.g$drift[[.what]]$DesignMatrix)

      .solution <- krigSolver(lhs, rhs, .DM, .X)

      .fs <- as.matrix(.g$data[[.what]]$functions)

      if(ncol(.fs) > 1){
        .forecast <-  .fs %*% .solution$weights
      } else {
        .forecast <-  t(t(.fs) %*% .solution$weights)
      }

      return(list(C = lhs, c = rhs, Variance = .solution$variance, Forecast = .forecast, Weights = .solution$weights))

    }
  }
}


#' @title A function that computes kriging weights and variance
#'
#' @param .LHS  -  Left hand side of the simple kriging matrix
#' @param .RHS  -  Right hand side of the simple kriging matrix
#' @param .F    -  Drift design matrix of the training data (ALL of it!)
#' @param .f    -  Drift design matrix of the target points.
#' @param .d    -  A matrix of distances between targets and training data
#'                 This is here to enable upgrate for variable neighborhood size (not implemented)
#' @param .k    -  Number of closest points to consider.
#'
#' @details  This function is using the method outlined in Chilles and Delfiner (1999) page 169 that overcomes
#' the problem with indefiniteness of ordinary or universal kriging matrices.
#' This function is also a syntactic sugar for predictFinterp function.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu})
#'
#' @export
#'
krigSolver <- function(.LHS, .RHS, .F = NULL, .f = NULL, .d = NULL, .k = NULL){

  # Sanity checks:
  if(is.null(.LHS) | is.null(.RHS)) stop('Left hand side or right hand side was not provided. Exiting!')
  if(nrow(.LHS) != nrow(.RHS)) stop('Left hand side is incompatible with the right hand side. Exiting!')

  if(is.null(.F) | is.null(.f)){ # simple kriging

    # Solve for weights: Cholesky should be implemented!

    .weights <- solve(.LHS) %*% .RHS

    # Solve for variance:
    .variance <- .LHS[1,1] - apply(.weights * .RHS, 2, sum)

    return(weights = .weights, variance = .variance)

  } else { # Universal or ordinary kriging.

    # This is per Chiles and Delfiner (1999), page 169.
    .SigmaInv <- solve(.LHS)
    .X1       <- .SigmaInv %*% .F
    .lambda   <- .SigmaInv %*% .RHS

    .Q  <- t(.F) %*% .X1
    .R  <- t(.F) %*% .lambda - t(.f)

    .mu <- solve(.Q) %*% .R

    .weights  <- .lambda - .X1 %*% .mu

    .variance <- .LHS[1,1] - apply(.weights * .RHS, 2, sum) - apply(.mu * t(.f), 2, sum)

    return(list(weights = .weights, variance = .variance))
  }

}


#' @title Performs k-fold cross validation on the kriging model:
#'
#' @param Coordinates     - Coordinates of your functional data
#' @param Functions       - Functions stacked in columns
#' @param Variogram Model - Class vgm, produced as a call to gstat::vgm
#' @param .covType        - Covariance type: "omni", "product", "sum".
#' @param  F               - Universal kriging matrix. For ordinary kriging a vector of 1's.
#'                          For simple kriging \code{NULL} which is default!
#'
#' @details This function performs k-fold cross validation for a given variogram model
#'          If F=NULL it assumes simple kriging in which case it must work with covariances.
#'          Covariances are produced from the variogram by substracting the variogram values
#'          from the fitted sill. Of course, this assumes second order stationarity! It's users
#'          responsibility to verify that this is truly appropriate!
#'
#'          In all other cases of universal kriging (OK as a special case), the code will
#'          work with variograms.
#'
#'          Note: This is an experimental function! Use at your own responsibility!
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu})
#'
#' @export
#'
traceKrigeCV <- function(.coordinates, .functions, .model, .F = NULL){

  if(nrow(.coordinates) != ncol(.functions)) {
    stop('The number of coordinates is not equal to the number of functions. Exiting!')
  }

  # Create a kriging matrix: Omni or Separable/Sum
  .d <- as.matrix(dist(as.matrix(.coordinates)))
  .G <- variogramLine(.model, dist_vector = .d)

  if(!is.null(.F)) .F <- as.matrix(.F)

  # Perform LOO x-Validation:
  .forecast <- matrix(NA, ncol = ncol(.functions), nrow = nrow(.functions))

  for(i in 1:ncol(.functions)){

    # Simple Kriging:
    if(is.null(.F)){

      .M <- .G[-i, ]
      .c <- .G[ , i]
      .M <- .M[ ,-i]
      .weights <- solve(.M) %*% .c

    } else {

      # this is per Chiles and Delfiner, page 169.
      .SigmaInv <- solve(.G[-i,-i])
      .X1       <- .SigmaInv %*% .F[-i,]
      .lambda   <- .SigmaInv %*% .G[-i, i]

      .Q  <- t(.F[-i,]) %*% .X1
      .R  <- t(.F[-i,]) %*% .lambda - .F[i, ]

      .mu <- solve(.Q) %*% .R

      .weights <- .lambda - .X1 %*% .mu

    }

    .forecast[,i] <- as.matrix(.functions[, -i]) %*% .weights

  }

  return(.forecast)

}










