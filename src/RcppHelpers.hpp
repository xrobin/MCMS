#pragma once

#include <Rcpp.h>

/** Things that are missing in Rcpp */

namespace Rcpp {
	/** Friendly rownames() and colnames() accessors of a NumericMatrix */
	CharacterVector rownames(const NumericMatrix &x);
	CharacterVector colnames(const NumericMatrix &x);
	void rownames(NumericMatrix &x, const CharacterVector& names);
	void colnames(NumericMatrix &x, const CharacterVector& names);
}
