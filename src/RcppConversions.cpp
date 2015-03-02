#include "conversions.hpp"
#include "LikelihoodParams.hpp"
#include <Rcpp.h>


namespace Rcpp {
	template <> oParams as(SEXP anOList) {
		oParams::o_type anOMap = convertListToOMap(as<List>(anOList));
		oParams anO(anOMap);
		return(anO);
	}

	template <> SEXP wrap(const oParams& anOParams) {
		throw std::runtime_error(std::string("not implemented"));
		//wrap(List)
	}
}