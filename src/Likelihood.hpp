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
	/** Constructors */
	PeptideLikelihood(const Peptide& aPeptide):
		peptide(aPeptide), likelihood(0) {};
	PeptideLikelihood(const Peptide& aPeptide, const double aLikelihood):
		peptide(aPeptide), likelihood(aLikelihood) {};
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
	std::vector<std::vector<PeptideLikelihood*>> onupdate_c;
	//   c           peptide indices
	// Maps an index on the o to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<std::vector<PeptideLikelihood*>>> onupdate_o;
	//  sample     site       peptide indices

	/** Link the peptides and parameters */
	void linkParamsAndPeptides();
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

	Likelihood(const std::vector<Peptide>&, const cParams&, const oParams&);
};
