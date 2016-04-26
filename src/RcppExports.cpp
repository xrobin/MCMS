// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getLikelihood_MCMC_Cpp
double getLikelihood_MCMC_Cpp(const S4& aModel, const List& aVarianceModelAsList, const NumericVector& aC, const List& anO, const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2, const bool verbose);
RcppExport SEXP MCMS_getLikelihood_MCMC_Cpp(SEXP aModelSEXP, SEXP aVarianceModelAsListSEXP, SEXP aCSEXP, SEXP anOSEXP, SEXP c_prior_sdSEXP, SEXP o_prior_shape1SEXP, SEXP o_prior_shape2SEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const S4& >::type aModel(aModelSEXP);
    Rcpp::traits::input_parameter< const List& >::type aVarianceModelAsList(aVarianceModelAsListSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type aC(aCSEXP);
    Rcpp::traits::input_parameter< const List& >::type anO(anOSEXP);
    Rcpp::traits::input_parameter< const double >::type c_prior_sd(c_prior_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type o_prior_shape1(o_prior_shape1SEXP);
    Rcpp::traits::input_parameter< const double >::type o_prior_shape2(o_prior_shape2SEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(getLikelihood_MCMC_Cpp(aModel, aVarianceModelAsList, aC, anO, c_prior_sd, o_prior_shape1, o_prior_shape2, verbose));
    return __result;
END_RCPP
}
// getPrior_MCMC_Cpp
double getPrior_MCMC_Cpp(const S4& aModel, const List& aVarianceModelAsList, const NumericVector& aC, const List& anO, const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2, const bool verbose);
RcppExport SEXP MCMS_getPrior_MCMC_Cpp(SEXP aModelSEXP, SEXP aVarianceModelAsListSEXP, SEXP aCSEXP, SEXP anOSEXP, SEXP c_prior_sdSEXP, SEXP o_prior_shape1SEXP, SEXP o_prior_shape2SEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const S4& >::type aModel(aModelSEXP);
    Rcpp::traits::input_parameter< const List& >::type aVarianceModelAsList(aVarianceModelAsListSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type aC(aCSEXP);
    Rcpp::traits::input_parameter< const List& >::type anO(anOSEXP);
    Rcpp::traits::input_parameter< const double >::type c_prior_sd(c_prior_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type o_prior_shape1(o_prior_shape1SEXP);
    Rcpp::traits::input_parameter< const double >::type o_prior_shape2(o_prior_shape2SEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(getPrior_MCMC_Cpp(aModel, aVarianceModelAsList, aC, anO, c_prior_sd, o_prior_shape1, o_prior_shape2, verbose));
    return __result;
END_RCPP
}
// run_MCMC_Cpp
NumericVector run_MCMC_Cpp(const S4& aModel, const List& aVarianceModelAsList, const unsigned long n, const unsigned long n_out, const unsigned long burn_in, const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2, const double o_restrict, const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale, const bool verbose, const bool cooling, const NumericVector& seed);
RcppExport SEXP MCMS_run_MCMC_Cpp(SEXP aModelSEXP, SEXP aVarianceModelAsListSEXP, SEXP nSEXP, SEXP n_outSEXP, SEXP burn_inSEXP, SEXP c_prior_sdSEXP, SEXP o_prior_shape1SEXP, SEXP o_prior_shape2SEXP, SEXP o_restrictSEXP, SEXP prior_move_proportionSEXP, SEXP c_sdSEXP, SEXP o_sdSEXP, SEXP o_k_scaleSEXP, SEXP verboseSEXP, SEXP coolingSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const S4& >::type aModel(aModelSEXP);
    Rcpp::traits::input_parameter< const List& >::type aVarianceModelAsList(aVarianceModelAsListSEXP);
    Rcpp::traits::input_parameter< const unsigned long >::type n(nSEXP);
    Rcpp::traits::input_parameter< const unsigned long >::type n_out(n_outSEXP);
    Rcpp::traits::input_parameter< const unsigned long >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const double >::type c_prior_sd(c_prior_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type o_prior_shape1(o_prior_shape1SEXP);
    Rcpp::traits::input_parameter< const double >::type o_prior_shape2(o_prior_shape2SEXP);
    Rcpp::traits::input_parameter< const double >::type o_restrict(o_restrictSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_move_proportion(prior_move_proportionSEXP);
    Rcpp::traits::input_parameter< const double >::type c_sd(c_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type o_sd(o_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type o_k_scale(o_k_scaleSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type cooling(coolingSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type seed(seedSEXP);
    __result = Rcpp::wrap(run_MCMC_Cpp(aModel, aVarianceModelAsList, n, n_out, burn_in, c_prior_sd, o_prior_shape1, o_prior_shape2, o_restrict, prior_move_proportion, c_sd, o_sd, o_k_scale, verbose, cooling, seed));
    return __result;
END_RCPP
}
