#pragma once

#include "Likelihood.hpp"
#include "Parameters.hpp"
#include <Rcpp.h>
#include "typedefs.hpp"

/** Various conversion functions outside the Rcpp namespace */

/** Converts a named list of numeric vectors into an o_type object (typically a map)
 * to represent occupancy parameters
 */
oParams::o_type convertListToOMap(const Rcpp::List&);
cParams::c_type convertVectorToCMap(const Rcpp::NumericVector&);
Likelihood convertS4ToLikelihood(const Rcpp::S4&, const VarianceModel&, const double scale, const double shape1, const double shape2,
	const double prior.move.proportion, const double c.sd, const double o.sd);


//vector<Peptides> convertListToPeptidesVector(List aDataList, ...) {