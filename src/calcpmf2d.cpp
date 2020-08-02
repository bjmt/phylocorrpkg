#include <Rcpp.h>

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector calcPMF2d(const Rcpp::IntegerVector &xInd,
    const Rcpp::IntegerVector &yInd, const Rcpp::NumericMatrix &PMFmat) {
  Rcpp::NumericVector out(xInd.size());
  for (R_xlen_t i = 0; i < xInd.size(); ++i) {
    out(i) = PMFmat(xInd[i], yInd[i]);
  }
  return out;
}
