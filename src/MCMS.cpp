#include "conversions.hpp"
#include "Parameters.hpp"
#include <Rcpp.h>
#include "RcppConversions.hpp"
#include "RcppHelpers.hpp"
#include <string>
#include <vector>
#include "VarianceModel.hpp"

using Rcpp::as;
using Rcpp::CharacterVector;
using Rcpp::S4;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using std::string;

// [[Rcpp::export]]
NumericVector run_MCMC_Cpp(const S4& aModel, const List& aVarianceModelAsList,
	const unsigned long n, const unsigned long n_out, const unsigned long burn_in,
	const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2, const double o_restrict,
	const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale, const bool verbose = false) {
	VarianceModel aVarianceModel = as<VarianceModel>(aVarianceModelAsList);

	MonteCarlo m = convertS4ToMonteCarlo(as<PeptidesModel>(aModel), aVarianceModel, c_prior_sd, o_prior_shape1, o_prior_shape2, o_restrict,
		prior_move_proportion, c_sd, o_sd, o_k_scale);

	if (verbose) {
		Rcpp::Rcout << m;
	}

	NumericMatrix res = m.iterate(n, n_out, burn_in);
	Rcpp::colnames(res, as<CharacterVector>(Rcpp::wrap(m.getIterateNames())));

	//Rcpp::Rcout << m.getCByReference() << std::endl;
	//Rcpp::Rcout << m.getOByReference() << std::endl;

	return res;
}

