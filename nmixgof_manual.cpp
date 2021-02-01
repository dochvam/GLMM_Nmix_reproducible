// THIS CODE REPURPOSED FROM KNAPE ET AL. 2018 (DOI: 10.1111/2041-210X.13062)

#include <Rcpp.h>
using namespace Rcpp;

// These functions compute the CDF of the sum of a sequence of
// independent binomial trials with the same index N. 


// [[Rcpp::export]]
NumericVector pbinsumRow(NumericVector y, double N, NumericVector p) {
  NumericVector res(3);
  NumericVector ySum(1);
  ySum(0) = -1;
  // Handle NAs by setting p to 0 and y  to 1 => dbinom(y,N,p) = 0.
  for (int i=0; i<y.length(); i++) {
    if(NumericVector::is_na(y[i]) || NumericVector::is_na(p[i])) {
      y[i] = 1;
      p[i] = 0.0;
    } else {
      if (ySum[0] == -1) {
        ySum[0] = y[i];
      } else {
        ySum[0] += y[i];
      }
    }
  }
  if (ySum[0] == -1) {
    res[0] = res[1] = res[2] = NA_REAL;
    return res;
  }
  if (!Rcpp::traits::is_finite<REALSXP>(ySum(0))) {
    res[0] = ySum[0];
    res[1] = res[2] = R_NaN;
    return res;
  }
  NumericMatrix pMat(ySum[0] + 1, p.length());
  NumericVector pT(ySum[0] + 1);
  NumericVector pF(ySum[0] + 1);
  for (int col=0; col<pMat.ncol(); col++) {
    for (int k = 0; k <= ySum[0]; k++) {
      pMat(k, col) = Rf_dbinom(k, N, p[col], false);
    }
  }
  for (int k = 0; k <= ySum[0]; k++) {
    pT[k] = pMat(k, 0);
    pF[k] = pT[k];
  }
  for (int col = 1; col<pMat.ncol(); col++) {
    for (int j = 0; j <= ySum[0]; j++) {
      pF[j] = 0;
      for (int i = 0; i<=j; i++) {
        pF[j] += pT[i]*pMat(j-i, col);
      }
    }
    for (int k = 0; k <= ySum[0]; k++) {
      pT[k] = pF[k];
    }
  }
  for (int k = 0; k < ySum[0]; k++)
    res[1] += pF[k];
  res[2] = res[1] + pF[ySum[0]];
  res[0] = ySum[0];
  return res; 
}



// [[Rcpp::export]]
NumericMatrix pbinsum(NumericMatrix y, NumericVector N, NumericMatrix p) {
  NumericMatrix cumProb(y.nrow(),3);
  if(y.ncol() != p.ncol() || y.nrow() != p.nrow()) {
    stop("Dimensions of y do not match those of p.");
  }
  if(y.nrow() != N.length()) {
    stop("Length of N does not match the number of rows of y or p.");
  }
  /*
   for (int i = 0; i < y.nrow(); i++) {
   for (int j = 0; j < y.ncol(); j++) {
   if (NumericVector::is_na(y(i,j)) || NumericVector::is_na(p(i,j))) {
   y(i,j) = NA_REAL;
   p(i,j) = NA_REAL;
   }
   }
   }
   */
  for (int i = 0; i < y.nrow(); i++) {
    cumProb.row(i) = pbinsumRow(y.row(i), N[i], p.row(i));
  }
  return cumProb;
}


#include <Rcpp.h>
using namespace Rcpp;

// These functions compute the CDF of the sum of a sequence of
// independent beta-binomial trials with the same index N.


// [[Rcpp::export]]
NumericVector pbbinsumRow(NumericVector y, double N, NumericVector p, NumericVector theta) {
  NumericVector res(3);
  NumericVector ySum(1);
  double vif = 1/theta[0];
  ySum(0) = -1;
  // Handle NAs by setting p to 0 and y  to 1 => dbinom(y,N,p) = 0.
  for (int i=0; i<y.length(); i++) {
    if(NumericVector::is_na(y[i]) || NumericVector::is_na(p[i])) {
      y[i] = 1;
      p[i] = -1.0;
    } else {
      if (ySum[0] == -1) {
        ySum[0] = y[i];
      } else {
        ySum[0] += y[i];
      }
    }
  }
  if (ySum[0] == -1) {
    res[0] = res[1] = res[2] = NA_REAL;
    return res;
  }
  if (!Rcpp::traits::is_finite<REALSXP>(ySum(0))) {
    res[0] = ySum[0];
    res[1] = res[2] = R_NaN;
    return res;
  }
  NumericMatrix pMat(ySum[0] + 1, p.length());
  NumericVector pT(ySum[0] + 1);
  NumericVector pF(ySum[0] + 1);
  for (int col=0; col<pMat.ncol(); col++) {
    for (int k = 0; k <= ySum[0]; k++) {
      if (p[col] < 0) {
        if (k==0)
          pMat(k, col) = 1.0;
        else
          pMat(k, col) = 0.0;
      } else {
        pMat(k, col) = exp(Rf_lchoose(N, k) + Rf_lbeta(k+p[col]*vif, N-k+(1-p[col])*vif) - Rf_lbeta(p[col]*vif, (1-p[col])*vif));
        if (NumericVector::is_na(pMat(k,col))) { // Rough trick to get around numerical underflow
          pMat(k, col) = 0.0;
        }
      }
    }
  }
  for (int k = 0; k <= ySum[0]; k++) {
    pT[k] = pMat(k, 0);
    pF[k] = pT[k];
  }
  for (int col = 1; col<pMat.ncol(); col++) {
    for (int j = 0; j <= ySum[0]; j++) {
      pF[j] = 0;
      for (int i = 0; i<=j; i++) {
        pF[j] += pT[i]*pMat(j-i, col);
      }
    }
    for (int k = 0; k <= ySum[0]; k++) {
      pT[k] = pF[k];
    }
  }
  for (int k = 0; k < ySum[0]; k++)
    res[1] += pF[k];
  res[2] = res[1] + pF[ySum[0]];
  res[0] = ySum[0];
  return res; 
}  

// [[Rcpp::export]]
NumericMatrix pbbinsum(NumericMatrix y, NumericVector N, NumericMatrix p, NumericVector theta) {
  NumericMatrix cumProb(y.nrow(),3);
  if(y.ncol() != p.ncol() || y.nrow() != p.nrow()) {
    stop("Dimensions of y do not match those of p.");
  }
  if(y.nrow() != N.length()) {
    stop("Length of N does not match the number of rows of y or p.");
  }
  if(theta.length() != 1) {
    stop("delta should have length 1.");
  }
  /*
   for (int i = 0; i < y.nrow(); i++) {
   for (int j = 0; j < y.ncol(); j++) {
   if (NumericVector::is_na(y(i,j)) || NumericVector::is_na(p(i,j))) {
   y(i,j) = NA_REAL;
   p(i,j) = NA_REAL;
   }
   }
   }
   */
  for (int i = 0; i < y.nrow(); i++) {
    cumProb.row(i) = pbbinsumRow(y.row(i), N[i], p.row(i), theta);
  }
  return cumProb;
}
