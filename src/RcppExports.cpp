// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calThruNodes_cpp
Rcpp::List calThruNodes_cpp(const NumericVector& ances, const NumericMatrix& states_R, const NumericMatrix& forTime_R, const NumericVector& lambdas, const NumericVector& mus, const NumericMatrix& Q, int num_threads);
RcppExport SEXP _secsseCPP_calThruNodes_cpp(SEXP ancesSEXP, SEXP states_RSEXP, SEXP forTime_RSEXP, SEXP lambdasSEXP, SEXP musSEXP, SEXP QSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type ances(ancesSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type states_R(states_RSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type forTime_R(forTime_RSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mus(musSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(calThruNodes_cpp(ances, states_R, forTime_R, lambdas, mus, Q, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_secsseCPP_calThruNodes_cpp", (DL_FUNC) &_secsseCPP_calThruNodes_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_secsseCPP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
