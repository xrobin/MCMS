#include "conversions.hpp"
#include "Parameters.hpp"
#include <Rcpp.h>
#include "RcppConversions.hpp"
#include "S4Aliases.hpp"
#include <vector>
#include "VarianceModel.hpp"

using Rcpp::as;
using Rcpp::S4;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// [[Rcpp::export]]
double getLikelihood_MCMC_Cpp(const S4& aModel, const List& aVarianceModelAsList,
	const NumericVector& aC, const List& anO,
	const double c_prior_sd, const double o_prior_shape1, const double o_prior_shape2,
	const bool verbose = false) {

	VarianceModel aVarianceModel = as<VarianceModel>(aVarianceModelAsList);
	PeptidesModel aPeptidesModel = as<PeptidesModel>(aModel);
	ProteinModel aProteinModel = aPeptidesModel.getProteinModel();
	const NumericMatrix sampleDependency = aProteinModel.slot("sample.dependency");

	// The following variables aren't used but still needed to build the MonteCarlo object.
	// We set them to signaling NaNs so we're sure we don't use them
	double o_restrict, prior_move_proportion, c_sd, o_sd, o_k_scale;
	o_restrict = prior_move_proportion = c_sd = o_sd = o_k_scale = std::numeric_limits<double>::signaling_NaN();
	Constants constants(aVarianceModel, sampleDependency, c_prior_sd, o_prior_shape1, o_prior_shape2,
		o_restrict, prior_move_proportion, c_sd, o_sd, o_k_scale);

	// Get the parameters
	const cParams::c_type aCMap = convertVectorToCMap(aC);
	const oParams::o_type anOMap = convertListToOMap(anO);

	cParams c(aCMap, constants.sampleDependenceMatrix);
	oParams o(anOMap);

	const std::vector<Peptide> peptides = convertS4ToPeptides(aProteinModel);

	// Make the likelihood object
	Likelihood l(peptides, c, o, constants);

	if (verbose) {
		Rcpp::Rcout << l;
	}

	return l.getLikelihoodValue();
}
