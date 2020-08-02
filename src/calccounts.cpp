#include <Rcpp.h>

// [[Rcpp::export(rng = false)]]
Rcpp::List calcCounts(const Rcpp::IntegerVector &xInd,
    const Rcpp::IntegerVector &yInd, const Rcpp::NumericMatrix &kernRatio,
    const Rcpp::NumericVector &kernRatioVec,
    const Rcpp::NumericVector &kernTrueNormVec,
    const Rcpp::NumericVector &kernFalseNormVec) {

  R_xlen_t n = kernRatioVec.size();

  Rcpp::NumericVector Tcounts(n), Fcounts(n);
  double r;

  for (R_xlen_t i = 0; i < n; ++i) {
    r = kernRatio(xInd[i], yInd[i]);
    Rcpp::LogicalVector rInd = kernRatioVec >= r;
    Rcpp::NumericVector Tsub = kernTrueNormVec[rInd];
    Rcpp::NumericVector Fsub = kernFalseNormVec[rInd];
    Tcounts[i] = Rcpp::sum(Tsub);
    Fcounts[i] = Rcpp::sum(Fsub);
  }

  return Rcpp::List::create(Rcpp::_["Tcount"] = Tcounts,
      Rcpp::_["Fcount"] = Fcounts);

}
