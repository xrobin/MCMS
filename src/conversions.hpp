#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include "typedefs.hpp"

/** Various conversion functions outside the Rcpp namespace */

/** Converts a named list of numeric vectors into an o_type object (typically a map)
 * to represent occupancy parameters
 */
oParams::o_type convertListToO(Rcpp::List);


//vector<Peptides> convertListToPeptidesVector(List aDataList, ...) {