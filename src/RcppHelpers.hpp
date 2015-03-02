#pragma once

#include <Rcpp.h>

/** Things that are missing in Rcpp */

namespace Rcpp {
	/** Friendly rownames() and colnames() accessors of a NumericMatrix */
	CharacterVector rownames(const NumericMatrix &x);
	CharacterVector colnames(const NumericMatrix &x);
}
