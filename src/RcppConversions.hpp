#pragma once

#include "Parameters.hpp"
#include <Rcpp.h>
#include "VarianceModel.hpp"

namespace Rcpp {
	// oParams
	template <> oParams::o_type as(SEXP);
	template <> SEXP wrap(const oParams&);

	// VarianceModel
	template <> VarianceModel as(SEXP);
	//template <> SEXP wrap(const VarianceModel<>&);

//	// LikelihoodParams
//	template <> LikelihoodParams as(Rcpp::S4 someParams);
//	template <> Rcpp::List wrap(const LikelihoodParams<> &someParams);
}
