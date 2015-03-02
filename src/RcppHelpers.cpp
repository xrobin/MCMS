#include <Rcpp.h>

/** Things that are missing in Rcpp */

namespace Rcpp {
	/** Friendly rownames() and colnames() accessors of a NumericMatrix */
	CharacterVector rownames(const NumericMatrix &x) {
		List dimnames = x.attr("dimnames");
		return dimnames[0];
	}
	CharacterVector colnames(const NumericMatrix &x) {
		List dimnames = x.attr("dimnames");
		return dimnames[1];
	}
}
