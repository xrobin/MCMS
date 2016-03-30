#include "conversions.hpp"
#include "RcppConversions.hpp"

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
		// Get shape
		NumericMatrix shape = as<List>(mdl)["shape"];
		double shape0 = shape(0, 0);

		// Get rate
		double rate0, rate1, rate2;
		VarianceModel::RateType rateType;
		NumericMatrix rate = as<List>(mdl)["rate"];
		if (rate.nrow() == 1) {
			rate0 = rate(0, 0);
			rateType = VarianceModel::RateType::constant;
		}
		else if (rate.nrow() == 3) {
			double rate0 = rate(1, 0);
			double rate1 = rate(2, 0);
			double rate2 = rate(0, 0);
			rateType = VarianceModel::RateType::power;
		}
		else {
			Rcpp::stop("'rate' matrix in variance model object must have 1 or 3 rows");
		}
		return VarianceModel(rateType, rate0, rate1, rate2, shape0);
	}

	template <> PeptidesModel as(SEXP aPeptidesModel) {
		return PeptidesModel(aPeptidesModel);
	}

	template <> ProteinModel as(SEXP aProteinModel) {
		return ProteinModel(aProteinModel);
	}
}