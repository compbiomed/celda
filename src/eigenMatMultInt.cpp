// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMatMultInt(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map< Eigen::MatrixXi> B){
  Eigen::MatrixXd C = A.transpose() * B.cast<double>();
  
  return Rcpp::wrap(C);
}
