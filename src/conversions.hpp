#pragma once

#include "Likelihood.hpp"
#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include "typedefs.hpp"

/** Various conversion functions outside the Rcpp namespace */

/** Converts a named list of numeric vectors into an o_type object (typically a map)
 * to represent occupancy parameters
 */
oParams::o_type convertListToOMap(const Rcpp::List&);
cParams::c_type convertVectorToCMap(const Rcpp::NumericVector&);
Likelihood convertS4ToLikelihood(const Rcpp::S4&, const VarianceModel&, const double shape1, const double shape2);


//vector<Peptides> convertListToPeptidesVector(List aDataList, ...) {