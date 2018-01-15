#pragma once

#include "Parameters.h"
#include <Rcpp.h>
#include "S4Aliases.h"
#include "VarianceModel.h"

namespace Rcpp {
	// oParams
	template <> oParams::o_type as(SEXP);
	template <> SEXP wrap(const oParams&);

	// VarianceModel
	template <> VarianceModel as(SEXP);
	//template <> SEXP wrap(const VarianceModel<>&);

	// PeptidesModel
	template <> PeptidesModel as(SEXP);
	template <> ProteinModel as(SEXP);
}
