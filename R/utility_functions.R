#' @title Numerical integration via trapezoid rule
#'
#' @details This is a support function for \code{\link{gammaTraceEmp}}. As such it should not be called directly!
#'
#' @param step - Step along the abscissa
#' @param y    - Function evaluated on an equispaced grid with \code{step = step}
#'
#' @return  Integral of \code{y}
#'
#' @author Alessandra Menafoglio (\email{alessandra.menafoglio@polimi.it})
#'
#' @references \url{https://en.wikipedia.org/wiki/Trapezoidal_rule}
#'
#' @export
#'
trapzc <- function(step, y)
{

  if(is.null(step)) stop('Step was not provided. Exiting!')

  if(is.null(y))    stop('Function was not provided! Exiting!')

  int <- step * (0.5 * y[1] + sum(y[2:(length(y) - 1)]) + 0.5 * y[length(y)])

  return (int)

}

#' @title  A simple scaling of a vector of parameters (i.e. to (0,1) scale)
#'
#' @param  .Input - Input vector you wish to scale.
#' @param  .lB    - lower bound you wish to use for scaling (default = min(.Input))
#' @param  .uB    - upper bound you wish to use for scaling (default = max(.Input))
#'
#' @author Ogy Grujic
#'
#' @export
#'
scaleInput <- function(.Input, .lB = NULL, .uB = NULL){

  if(is.null(.lB)) {
    .lB   <- min(.Input)
  }

  if(is.null(.uB)) {
    .uB   <- max(.Input)
  }

  .output <- (.Input - .lB) / (.uB - .lB)

  return(.output)

}

#'@title Plots summary statistics for regression model(s)
#'
#'@param .g - An fstat structure
#'
#'@author Ogy Grujic (ogyg@stanford.edu)
#'
#'@export
#'
summary.fstat <- function(.g=NULL){

  if(is.null(.g$drift)) stop("Drift not present in the .g structure! Exiting")
  if(ncol(.g$data[[1]]$functions) > 1) stop("Summary statistics are available only for scalar outputs")

  for(.name in names(.g$drift)){
    .Betas <- t(.g$drift[[.name]]$Betas)
    .X <- as.matrix(.g$drift[[.name]]$DesignMatrix)

    if(is.null(.g$covariance)){
      .Sigma <- diag(nrow(.X))
    } else {
      .Sigma <- .g$covariance$omni[[.name]]
    }
    .rss <- sum(.g$drift[[.name]]$Residuals^2)
    .rdf <- nrow(.X) - length(.Betas)
    .var <- (.rss/.rdf) * solve(t(.X) %*% solve(.Sigma) %*% .X)
    .sd <- sqrt(diag(.var))
    .summary <- cbind(.Betas, .sd, .Betas/.sd, 2*pt(abs(.Betas/.sd), df=.rdf, lower.tail = FALSE))
    colnames(.summary) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    cat(paste('Drift summary for variable ', .name, "\n", sep=""))
    printCoefmat(.summary, P.values=TRUE)
  }
  # 2*pt(abs(-2.105), 188, lower.tail = FALSE)
}


