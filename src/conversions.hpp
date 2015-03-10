#pragma once

//#include "Likelihood.hpp"
#include "MonteCarlo.hpp"
#include "Parameters.hpp"
#include <Rcpp.h>
#include "typedefs.hpp"

/** Various conversion functions outside the Rcpp namespace */

/** Converts a named list of numeric vectors into an o_type object (typically a map)
 * to represent occupancy parameters
 */
oParams::o_type convertListToOMap(const Rcpp::List&);
cParams::c_type convertVectorToCMap(const Rcpp::NumericVector&);
MonteCarlo convertS4ToMonteCarlo(const Rcpp::S4&, const VarianceModel&, const double scale, const double shape1, const double shape2,
	const double o_restrict,
	const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale);


//vector<Peptides> convertListToPeptidesVector(List aDataList, ...) {


/** Initialize a PRNG using the seed from R */
std::mt19937_64 seedFromR();

