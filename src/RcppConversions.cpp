#include "conversions.hpp"
#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include "VarianceModel.hpp"


namespace Rcpp {
	template <> oParams::o_type as(SEXP anOList) {
		return convertListToOMap(as<List>(anOList));
	}

	template <> SEXP wrap(const oParams& anOParams) {
		throw std::runtime_error(std::string("not implemented"));
		//wrap(List)
	}

	template <> VarianceModel as(SEXP aVarianceModelList) {
		List mdl = as<List>(aVarianceModelList)["mdl"];
		NumericMatrix shape = as<List>(mdl)["shape"];
		NumericMatrix rate = as<List>(mdl)["rate"];
		double shape0 = shape(0, 0);
		double rate0 = rate(1, 0);
		double rate1 = rate(2, 0);
		double rate2 = rate(0, 0);
		return VarianceModel(rate0, rate1, rate2, shape0);
	}
}