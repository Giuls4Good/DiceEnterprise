// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// findIntervalSingle
int findIntervalSingle(double U, NumericVector breaks);
RcppExport SEXP _DiceEnterprise_findIntervalSingle(SEXP USEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type U(USEXP);
    Rcpp::traits::input_parameter< NumericVector >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(findIntervalSingle(U, breaks));
    return rcpp_result_gen;
END_RCPP
}
// updateFunCpp
int updateFunCpp(int currentState, arma::uvec B, arma::vec U, bool connected, bool fine, List P_cumsum, List P_moves_list);
RcppExport SEXP _DiceEnterprise_updateFunCpp(SEXP currentStateSEXP, SEXP BSEXP, SEXP USEXP, SEXP connectedSEXP, SEXP fineSEXP, SEXP P_cumsumSEXP, SEXP P_moves_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type currentState(currentStateSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type U(USEXP);
    Rcpp::traits::input_parameter< bool >::type connected(connectedSEXP);
    Rcpp::traits::input_parameter< bool >::type fine(fineSEXP);
    Rcpp::traits::input_parameter< List >::type P_cumsum(P_cumsumSEXP);
    Rcpp::traits::input_parameter< List >::type P_moves_list(P_moves_listSEXP);
    rcpp_result_gen = Rcpp::wrap(updateFunCpp(currentState, B, U, connected, fine, P_cumsum, P_moves_list));
    return rcpp_result_gen;
END_RCPP
}
// CFTPCpp
List CFTPCpp(int k, NumericVector probs, bool connected, bool fine, List P_cumsum, List P_moves_list, bool monotonic, int min, int max, bool verbose);
RcppExport SEXP _DiceEnterprise_CFTPCpp(SEXP kSEXP, SEXP probsSEXP, SEXP connectedSEXP, SEXP fineSEXP, SEXP P_cumsumSEXP, SEXP P_moves_listSEXP, SEXP monotonicSEXP, SEXP minSEXP, SEXP maxSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< bool >::type connected(connectedSEXP);
    Rcpp::traits::input_parameter< bool >::type fine(fineSEXP);
    Rcpp::traits::input_parameter< List >::type P_cumsum(P_cumsumSEXP);
    Rcpp::traits::input_parameter< List >::type P_moves_list(P_moves_listSEXP);
    Rcpp::traits::input_parameter< bool >::type monotonic(monotonicSEXP);
    Rcpp::traits::input_parameter< int >::type min(minSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(CFTPCpp(k, probs, connected, fine, P_cumsum, P_moves_list, monotonic, min, max, verbose));
    return rcpp_result_gen;
END_RCPP
}
// construct_discrete_simplex
arma::mat construct_discrete_simplex(int d, int m);
RcppExport SEXP _DiceEnterprise_construct_discrete_simplex(SEXP dSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_discrete_simplex(d, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DiceEnterprise_findIntervalSingle", (DL_FUNC) &_DiceEnterprise_findIntervalSingle, 2},
    {"_DiceEnterprise_updateFunCpp", (DL_FUNC) &_DiceEnterprise_updateFunCpp, 7},
    {"_DiceEnterprise_CFTPCpp", (DL_FUNC) &_DiceEnterprise_CFTPCpp, 10},
    {"_DiceEnterprise_construct_discrete_simplex", (DL_FUNC) &_DiceEnterprise_construct_discrete_simplex, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DiceEnterprise(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
