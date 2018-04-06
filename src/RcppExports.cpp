// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GetProxOne
arma::vec GetProxOne(arma::vec y, arma::vec weights);
RcppExport SEXP _hierbasis2_GetProxOne(SEXP ySEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(GetProxOne(y, weights));
    return rcpp_result_gen;
END_RCPP
}
// FitAdditive
List FitAdditive(arma::vec y, arma::mat weights, arma::vec ak, NumericVector x, arma::mat beta, double max_lambda, double lam_min_ratio, double alpha, double tol, int p, int J, int n, int nlam, double max_iter, bool beta_is_zero, arma::vec active_set);
RcppExport SEXP _hierbasis2_FitAdditive(SEXP ySEXP, SEXP weightsSEXP, SEXP akSEXP, SEXP xSEXP, SEXP betaSEXP, SEXP max_lambdaSEXP, SEXP lam_min_ratioSEXP, SEXP alphaSEXP, SEXP tolSEXP, SEXP pSEXP, SEXP JSEXP, SEXP nSEXP, SEXP nlamSEXP, SEXP max_iterSEXP, SEXP beta_is_zeroSEXP, SEXP active_setSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ak(akSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type max_lambda(max_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lam_min_ratio(lam_min_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< double >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type beta_is_zero(beta_is_zeroSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type active_set(active_setSEXP);
    rcpp_result_gen = Rcpp::wrap(FitAdditive(y, weights, ak, x, beta, max_lambda, lam_min_ratio, alpha, tol, p, J, n, nlam, max_iter, beta_is_zero, active_set));
    return rcpp_result_gen;
END_RCPP
}
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
    {"_hierbasis2_GetProxOne", (DL_FUNC) &_hierbasis2_GetProxOne, 2},
    {"_hierbasis2_FitAdditive", (DL_FUNC) &_hierbasis2_FitAdditive, 16},
    {"_hierbasis2_prox", (DL_FUNC) &_hierbasis2_prox, 2},
    {"_hierbasis2_solvehierbasis", (DL_FUNC) &_hierbasis2_solvehierbasis, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_hierbasis2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
