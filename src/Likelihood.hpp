#pragma once
#include <map>
#include <Rcpp.h>
#include <vector>

#include "LikelihoodParams.hpp"
#include "Peptide.hpp"

class PeptideLikelihood {
	public:
	Peptide peptide;
	double likelihood;
};

class Likelihood {
	std::vector<PeptideLikelihood> peptideLikelihoods;
	//std::vector<Peptide> peptides;
	//std::vector<double> likelihoods;
	double likelihood;

	// Store likelihood parameters
	cParams c;
	oParams o;


	// mapping from param to peptides that must be updated
	// Maps an index on the c to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<oPair*>> onupdate_c;
	//   c           peptide indices
	// Maps an index on the o to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<std::vector<oPair*>>> onupdate_o;
	//  sample     site       peptide indices

	/** Computes the whole likelihood */
	void computeLikelihood();
	//void bindPeptidesAndParameters();

	public:
//	Likelihood(
//		std::vector<Peptide> somePeptides,
//		LikelihoodParams someLikelihoodParams
//		):
//		peptides(somePeptides),
//		likelihoods(somePeptides.size()),
//		likelihood(0), // We'll update it just after
//		params(someLikelihoodParams),
//		onupdate_c(), onupdate_o() // We'll update those too...
//	{
//		bindPeptidesAndParameters();
//		computeLikelihood();
//	};

	Likelihood(const Rcpp::S4& aProtein, const Rcpp::S4&);
};
