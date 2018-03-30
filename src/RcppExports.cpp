// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// prox
arma::sp_mat prox(arma::vec u, arma::mat w);
RcppExport SEXP _hierbasis2_prox(SEXP uSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(prox(u, w));
    return rcpp_result_gen;
END_RCPP
}
// solvehierbasis
List solvehierbasis(arma::mat design_mat, arma::vec y, arma::vec ak, arma::mat weights, int n, double lam_min_ratio, int nlam, double max_lambda);
RcppExport SEXP _hierbasis2_solvehierbasis(SEXP design_matSEXP, SEXP ySEXP, SEXP akSEXP, SEXP weightsSEXP, SEXP nSEXP, SEXP lam_min_ratioSEXP, SEXP nlamSEXP, SEXP max_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type design_mat(design_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ak(akSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lam_min_ratio(lam_min_ratioSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< double >::type max_lambda(max_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(solvehierbasis(design_mat, y, ak, weights, n, lam_min_ratio, nlam, max_lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hierbasis2_prox", (DL_FUNC) &_hierbasis2_prox, 2},
    {"_hierbasis2_solvehierbasis", (DL_FUNC) &_hierbasis2_solvehierbasis, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_hierbasis2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
