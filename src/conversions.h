#pragma once

//#include "Likelihood.h"
#include "MonteCarlo.h"
#include "Parameters.h"
#include "Peptide.h"
#include <Rcpp.h>
#include "S4Aliases.h" // ProteinModel
#include "typedefs.h"

/** Various conversion functions outside the Rcpp namespace */

/** Converts a named list of numeric vectors into an o_type object (typically a map)
 * to represent occupancy parameters
 */
oParams::o_type convertListToOMap(const Rcpp::List&);
cParams::c_type convertVectorToCMap(const Rcpp::NumericVector&);

/** Converts the ProteinModel to a vector of Peptide */
std::vector<Peptide> convertS4ToPeptides(const ProteinModel& aProtein);

/** The actual MC objects */
MonteCarlo convertS4ToMonteCarloWithParams(const PeptidesModel&, const VarianceModel&,
		const cParams::c_type aCMap, const oParams::o_type anOMap,
		const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2,
		const double o_restrict, const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale,
		const Rcpp::NumericVector& seed);

MonteCarlo convertS4ToMonteCarlo(const PeptidesModel&, const VarianceModel&, const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2,
	const double o_restrict,
	const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale,
	const Rcpp::NumericVector& seed);


//vector<Peptides> convertListToPeptidesVector(List aDataList, ...) {


/** Initialize a PRNG using the seed from R */
std::mt19937_64 seedFromR(const Rcpp::NumericVector&);

