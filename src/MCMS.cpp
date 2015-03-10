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
NumericVector run_MCMC_Cpp(const S4& aModel, const List& aVarianceModelAsList, const unsigned long n, const unsigned long n_out,
	const double scale, const double shape1, const double shape2, const double o_restrict,
	const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale) {
	VarianceModel aVarianceModel = as<VarianceModel>(aVarianceModelAsList);

	MonteCarlo m = convertS4ToMonteCarlo(aModel, aVarianceModel, scale, shape1, shape2, o_restrict,
		prior_move_proportion, c_sd, o_sd, o_k_scale);

	Rcpp::Rcout << m;

	NumericMatrix res = m.iterate(n, n_out);
	Rcpp::colnames(res, as<CharacterVector>(Rcpp::wrap(m.getIterateNames())));

	Rcpp::Rcout << m.getCByReference() << std::endl;
	Rcpp::Rcout << m.getOByReference() << std::endl;

	return res;
}
