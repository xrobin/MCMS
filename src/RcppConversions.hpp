#pragma once

#include "LikelihoodParams.hpp"
#include <Rcpp.h>

namespace Rcpp {
	// oParams
	template <> oParams as(SEXP&);
	template <> SEXP wrap(const oParams<>&);

//	// LikelihoodParams
//	template <> LikelihoodParams as(Rcpp::S4 someParams);
//	template <> Rcpp::List wrap(const LikelihoodParams<> &someParams);
}
