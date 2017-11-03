// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

/* Computes if two vectors are parallel within given tolerance (degrees) */
int isPar (arma::vec x,
           arma::vec y, double tol){

  int out ;
  double const pi = 3.14 ;

  double dotProd = dot(x,y) ;
  double denom   = norm(x) * norm(y) ;

  double angle = acos(dotProd/denom) * (180/pi) ;

  if(angle <= tol){
    out = 1 ;
  } else {
    out = 0 ;
  }

  return(out) ;

}

/* Returns a vector with 1's and 0's indicating whether the x is parallel (or not)
* within given tolerance (in degrees) to any of the row vectors provided in matrix Y.
* It works with any direction in any dimmension!
*/

NumericVector isParallel (NumericVector x,
                          NumericMatrix Y,
                          double tol = 1e-10){

  int nDir = Y.ncol() ;

  NumericVector outVect(nDir) ;

  for(int k = 0; k < nDir; ++k){

    NumericMatrix::Column y = Y.column(k);

    outVect(k) = isPar(x, y, tol) ;

  }

  return(outVect) ;
}

/* Computes Euclidean distance between vectors a and b. */
double Distance(NumericVector a, NumericVector b){

  double d = 0;

  for (int i = 0; i < a.size(); ++i){

    d += (a[i] - b[i])*(a[i] - b[i]) ;

  }

  d = sqrt(d);
  return d;
}


/* Computes partial functional variogram as per Gromenko and Kokoszka */
double gammaPartial(NumericVector f1, NumericVector f2){

  double gammaOut = 0 ;
  double counter  = 0 ;
  double tempSum  = 0 ;

  /* It assumes that both vectors are of the same length */
  for(int j = 0; j < f1.size(); ++j){

    if(f1(j) == f1(j) && f2(j) == f2(j)){

      tempSum += (f1(j) - f2(j)) * (f1(j) - f2(j)) ;
      counter += 1 ;

    }

  }

  gammaOut = tempSum / counter ;

  return gammaOut ;

}

/* Computes Integrated squared difference between two functions. Function
 * uses trapezoidal method for integration.
 */

double gammaTrace(NumericVector f1, NumericVector f2, const double argStep){

  double gammaTrace = 0;

  /* Compute squared difference */
    NumericVector y(f1.size()) ;

  for(int j=0; j < f1.size(); ++j){

    y(j) = (f1(j) - f2(j))*(f1(j) - f2(j)) ;

  }

  /*multiply the first and the last with 0.5 */
    gammaTrace = 0.5*(y(0) + y(y.size()-1)) ;


    for (int i=1; i< y.size() - 1; ++i){

      gammaTrace += y(i) ;

    }

    gammaTrace *= argStep;

    return gammaTrace;

}

double gammaCrossTrace(NumericVector pf1, NumericVector pf2,
                       NumericVector sf1, NumericVector sf2,
                       const double argStep){

  double gammaTrace = 0;

  /* Compute squared difference */
  NumericVector y(pf1.size());

  for(int j=0; j < pf1.size(); ++j){
    y(j) = (pf1(j) - pf2(j))*(sf1(j) - sf2(j));
  }

  /*multiply the first and the last with 0.5 */
  gammaTrace = 0.5*(y(0) + y(y.size()-1));


  for (int i = 1; i < y.size() - 1; ++i){

    gammaTrace += y(i);

  }

  gammaTrace *= argStep;

  return gammaTrace;

}

double gammaPseudoCrossTrace(NumericVector pf,
                             NumericVector sf,
                             const double argStep){

  double gammaTrace = 0;

  /* Compute squared difference */
  NumericVector y(pf.size());

  for(int j=0; j < pf.size(); ++j){
    y(j) = (pf(j) - sf(j)) * (pf(j) - sf(j));
  }

  /*multiply the first and the last with 0.5 */
  gammaTrace = 0.5*(y(0) + y(y.size()-1));


  for (int i = 1; i < y.size() - 1; ++i){

    gammaTrace += y(i);

  }

  gammaTrace *= argStep;

  return gammaTrace;

}

// [[Rcpp::export]]
NumericMatrix empTraceVariogram(NumericMatrix Coords,
                                NumericMatrix Fs,
                                NumericVector Lags,
                                NumericMatrix Directions,
                                double distTolerance,
                                double Step = 0,
                                double angTolerance = 0,
                                bool   isPartial    = false){

  double d ;
  float gamma ;
  float sf1 ;
  float sf2 ;
  int    tempK ;

  int nCoords = Coords.nrow() ;
  int nLags   = Lags.size() ;
  int nDir    = Directions.ncol() ;
  int nArgVal = Fs.ncol() ;

  NumericVector parVect(nDir) ;
  NumericMatrix outMatrix(nLags * nDir, 5);

  /* Checkups: */

  for (int i = 0; i < nCoords; ++i){

    for(int j = 0; j < i; ++j){

      NumericMatrix::Row a = Coords.row(i);
      NumericMatrix::Row b = Coords.row(j);

      d       = Distance(a, b) ;
      parVect = isParallel(a - b, Directions, angTolerance) ;

      /* Two modes single variate and functional */
      if(nArgVal > 1){

        /* Compute gammaTrace for current pair */
        NumericMatrix::Column f1 = Fs.column(i);
        NumericMatrix::Column f2 = Fs.column(j);


        if(isPartial == true){

          /* Compute average gammaTrace for current pair on all avaialbe time steps! */
          gamma = gammaPartial(f1, f2);

        } else {

          /* Compute complete gamma with trapezoidal method */
          gamma = gammaTrace(f1, f2, Step);
        }

      } else { // Scalar case

        sf1 = Fs(i, 0);
        sf2 = Fs(j, 0);

        gamma = (sf1 - sf2)*(sf1 - sf2) ;

      }

      /* Determine lag */
      for(int k = 0; k < nLags; ++k){

        if(d > (Lags(k) - distTolerance) && d < (Lags(k) + distTolerance)){

          /* Omni Directional Variogram */
          outMatrix(k, 0)     = Lags(k) ;
          outMatrix(k, 1)     = outMatrix(k, 1) + 1 ;
          outMatrix(k, 2)     = outMatrix(k, 2) + gamma ;
          outMatrix(k, 3)     = 0 ;
          outMatrix(k, 4)     = outMatrix(k, 4) + d ;

          /* See if this gamma fits into any other direction besides omni */
          for(int u = 1 ; u < nDir ; ++u){

            if(parVect(u) == 1){

              tempK = u * nLags + k ;
              outMatrix(tempK, 0) = Lags(k) ;
              outMatrix(tempK, 1) = outMatrix(tempK, 1) + 1 ;
              outMatrix(tempK, 2) = outMatrix(tempK, 2) + gamma ;
              outMatrix(tempK, 3) = u ;
              outMatrix(tempK, 4) = outMatrix(tempK, 4) + d ;

            }

          }

        }
      }
    }
  }

  outMatrix(_, 2) = outMatrix(_, 2) / (2*outMatrix(_, 1));
  outMatrix(_, 4) = outMatrix(_, 4) / (outMatrix(_, 1));

  return outMatrix;
}


// [[Rcpp::export]]
NumericMatrix empTraceCrossVariogram(NumericMatrix Coords,
                                     NumericMatrix Fs,
                                     NumericMatrix SFs,
                                     NumericVector Lags,
                                     NumericMatrix Directions,
                                     double distTolerance,
                                     double Step = 0,
                                     double angTolerance = 0
){

  double d ;
  double gamma ;
  double pf1 ;
  double pf2 ;
  double sf1 ;
  double sf2 ;
  int    tempK ;

  int nCoords = Coords.nrow() ;
  int nLags   = Lags.size() ;
  int nDir    = Directions.ncol() ;
  int nArgVal = Fs.ncol() ;

  NumericVector parVect(nDir) ;
  NumericMatrix outMatrix(nLags * nDir, 5);

  /* Checkups: */

  for (int i = 0; i < nCoords; ++i){

    for(int j = 0; j < i; ++j){

      NumericMatrix::Row a = Coords.row(i);
      NumericMatrix::Row b = Coords.row(j);

      d       = Distance(a, b) ;
      parVect = isParallel(a - b, Directions, angTolerance) ;

      /* Two modes single variate and functional */
      if(nArgVal > 1){

        /* Compute gammaTrace for current pair */
        NumericMatrix::Column pf1  = Fs.column(i);
        NumericMatrix::Column pf2  = Fs.column(j);
        NumericMatrix::Column sf1 = SFs.column(i);
        NumericMatrix::Column sf2 = SFs.column(j);

        gamma = gammaCrossTrace(pf1, pf2, sf1, sf2, Step);

      } else {

        pf1 = Fs(i, 0);
        pf2 = Fs(j, 0);
        sf1 = SFs(i, 0);
        sf2 = SFs(j, 0);

        gamma = (pf1 - pf2)*(sf1 - sf2) ;

      }



      /* Determine lag */
      for(int k = 0; k < nLags; ++k){

        if(d > (Lags(k) - distTolerance) && d < (Lags(k) + distTolerance)){

          /* Omni Directional Variogram */
          outMatrix(k, 0)     = Lags(k) ;
          outMatrix(k, 1)     = outMatrix(k, 1) + 1 ;
          outMatrix(k, 2)     = outMatrix(k, 2) + gamma ;
          outMatrix(k, 3)     = 0 ;
          outMatrix(k, 4) = outMatrix(k, 4) + d ;

          /* See if this gamma fits into any other direction besides omni */
          for(int u = 1 ; u < nDir ; ++u){

            if(parVect(u) == 1){

              tempK = u * nLags + k ;
              outMatrix(tempK, 0) = Lags(k) ;
              outMatrix(tempK, 1) = outMatrix(tempK, 1) + 1 ;
              outMatrix(tempK, 2) = outMatrix(tempK, 2) + gamma ;
              outMatrix(tempK, 3) = u ;
              outMatrix(tempK, 4) = outMatrix(tempK, 4) + d ;



            }

          }

        }
      }
    }
  }

  /* Computes cross-gamma */
  outMatrix(_, 2) = outMatrix(_, 2) / (2 * outMatrix(_, 1));
  outMatrix(_, 4) = outMatrix(_, 4) / (outMatrix(_, 1));

  return outMatrix;
}


// [[Rcpp::export]]
NumericMatrix isoCoordinates(NumericMatrix CoordsF1,
                             NumericMatrix CoordsF2){

  int nCoordsF1   = CoordsF1.nrow();
  int nCoordsF2   = CoordsF2.nrow();
  int isoElements = 0;

  NumericMatrix isoIndices(nCoordsF1, 2);

  for(int i = 0; i < nCoordsF1; i++){

    isoIndices(i, 0) = i + 1;
    isoIndices(i, 1) = NAN;

    NumericMatrix::Row r1 = CoordsF1.row(i);

    for(int j=0; j<nCoordsF2; j++){
      NumericMatrix::Row r2 = CoordsF2.row(j);



      if(Distance(r1,r2) == 0){
        isoIndices(i, 1) = j + 1;
        isoElements += 1;
        break;
      }
    }
  }

  return(isoIndices);
}


// [[Rcpp::export]]
NumericMatrix empTracePseudoCrossVariogram(NumericMatrix primaryCoords,
                                           NumericMatrix Fs,
                                           NumericMatrix secondaryCoords,
                                           NumericMatrix SFs,
                                           NumericVector Lags,
                                           NumericMatrix Directions,
                                           double distTolerance,
                                           double Step = 0,
                                           double angTolerance = 0
){

  double d ;
  double gamma ;
  double pf ;
  double sf ;
  int    tempK ;

  int nPrimaryCoords   = primaryCoords.nrow() ;
  int nSecondaryCords  = secondaryCoords.nrow() ;
  int nLags            = Lags.size() ;
  int nDir             = Directions.ncol() ;
  int nArgVal          = Fs.ncol() ;

  NumericVector parVect(nDir) ;
  NumericMatrix outMatrix(nLags * nDir, 5);

  /* Checkups: */

  for (int i = 0; i < nPrimaryCoords; ++i){ // Loop through the primary points

    for(int j = 0; j < nSecondaryCords; ++j){ //Loop through the secondary points

      NumericMatrix::Row a = primaryCoords.row(i);
      NumericMatrix::Row b = secondaryCoords.row(j);

      d       = Distance(a, b) ;
      parVect = isParallel(a - b, Directions, angTolerance) ;

      /* Two modes single variate and functional */
      if(nArgVal > 1){

        /* Compute gammaTrace for current pair */
        NumericMatrix::Column pf  = Fs.column(i);
        NumericMatrix::Column sf  = SFs.column(j);

        gamma = gammaPseudoCrossTrace(pf, sf, Step);

      } else {

        pf = Fs(i, 0);
        sf = SFs(j, 0);

        gamma = (pf - sf)*(pf - sf) ;

      }



      /* Determine lag */
      for(int k = 0; k < nLags; ++k){

        if(d > (Lags(k) - distTolerance) && d < (Lags(k) + distTolerance)){

          /* Omni Directional Variogram */
          outMatrix(k, 0)     = Lags(k) ;
          outMatrix(k, 1)     = outMatrix(k, 1) + 1 ;
          outMatrix(k, 2)     = outMatrix(k, 2) + gamma ;
          outMatrix(k, 3)     = 0 ;
          outMatrix(k, 4)     = outMatrix(k, 4) + d ;

          /* See if this gamma fits into any other direction besides omni */
          for(int u = 1 ; u < nDir ; ++u){

            if(parVect(u) == 1){

              tempK = u * nLags + k ;
              outMatrix(tempK, 0) = Lags(k) ;
              outMatrix(tempK, 1) = outMatrix(tempK, 1) + 1 ;
              outMatrix(tempK, 2) = outMatrix(tempK, 2) + gamma ;
              outMatrix(tempK, 3) = u ;
              outMatrix(tempK, 4) = outMatrix(tempK, 4) + d ;



            }

          }

        }
      }
    }
  }

  /* Computes cross-gamma */
  outMatrix(_, 2) = outMatrix(_, 2) / (2 * outMatrix(_, 1));
  outMatrix(_, 4) = outMatrix(_, 4) / (outMatrix(_, 1));

  return outMatrix;
}



