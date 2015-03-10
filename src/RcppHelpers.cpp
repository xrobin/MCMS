#include <Rcpp.h>
#include "RcppHelpers.hpp"

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

	/** Friendly rownames() and colnames() accessors of a NumericMatrix */
	void rownames(NumericMatrix &x, const CharacterVector& names) {
		List dimnames = x.attr("dimnames");
		if (dimnames.size() == 0) {
			dimnames = Rcpp::List::create(names, CharacterVector());
		}
		else {
			dimnames[0] = names;
		}
		x.attr("dimnames") = dimnames;
	}
	void colnames(NumericMatrix &x, const CharacterVector& names) {
		List dimnames = x.attr("dimnames");
		if (dimnames.size() == 0) {
			dimnames = Rcpp::List::create(CharacterVector(), names);
		}
		else {
			dimnames[1] = names;
		}
		x.attr("dimnames") = dimnames;
	}
}
