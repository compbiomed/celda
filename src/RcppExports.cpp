// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// eigenMatMultInt
SEXP eigenMatMultInt(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map< Eigen::MatrixXi> B);
RcppExport SEXP _celda_eigenMatMultInt(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map< Eigen::MatrixXi> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMultInt(A, B));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _colSumByGroup(SEXP, SEXP);
RcppExport SEXP _colSumByGroup_real(SEXP, SEXP);
RcppExport SEXP _colSumByGroupChange(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP _mvAdd(SEXP, SEXP, SEXP);
RcppExport SEXP _mvMult(SEXP, SEXP, SEXP);
RcppExport SEXP _rowSumByGroup(SEXP, SEXP);
RcppExport SEXP _rowSumByGroup_real(SEXP, SEXP);
RcppExport SEXP _rowSumByGroupChange(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_celda_eigenMatMultInt", (DL_FUNC) &_celda_eigenMatMultInt, 2},
    {"_colSumByGroup",         (DL_FUNC) &_colSumByGroup,         2},
    {"_colSumByGroup_real",    (DL_FUNC) &_colSumByGroup_real,    2},
    {"_colSumByGroupChange",   (DL_FUNC) &_colSumByGroupChange,   4},
    {"_mvAdd",                 (DL_FUNC) &_mvAdd,                 3},
    {"_mvMult",                (DL_FUNC) &_mvMult,                3},
    {"_rowSumByGroup",         (DL_FUNC) &_rowSumByGroup,         2},
    {"_rowSumByGroup_real",    (DL_FUNC) &_rowSumByGroup_real,    2},
    {"_rowSumByGroupChange",   (DL_FUNC) &_rowSumByGroupChange,   4},
    {NULL, NULL, 0}
};

RcppExport void R_init_celda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
