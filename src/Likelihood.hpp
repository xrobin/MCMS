#pragma once

#include <iostream> // cout
#include "LikelihoodParams.hpp"
#include <map>
#include "Peptide.hpp"
#include <Rcpp.h>
#include <vector>


class PeptideLikelihood {
	public:
	Peptide peptide;
	double likelihood;
	/** Constructors */
	PeptideLikelihood(const Peptide& aPeptide):
		peptide(aPeptide), likelihood(0) {};
	PeptideLikelihood(const Peptide& aPeptide, const double aLikelihood):
		peptide(aPeptide), likelihood(aLikelihood) {};
	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const PeptideLikelihood &aPeptideLikelihood);
};


class Likelihood {
	std::vector<PeptideLikelihood> peptideLikelihoods;
	//std::vector<Peptide> peptides;
	//std::vector<double> likelihoods;
	double likelihood;

	// Store likelihood parameters
	cParams c;
	oParams o;
	const LikelihoodConstants constants;


	// mapping from param to peptides that must be updated
	// Maps an index on the c to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<PeptideLikelihood*>> onupdate_c;
	//   c           peptide indices
	// Maps an index on the o to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<std::vector<PeptideLikelihood*>>> onupdate_o;
	//  sample     site       peptide indices

	/** Link the peptides and parameters */
	void linkParamsAndPeptides();
	/** Computes the whole likelihood and stores it */
	void updateLikelihood();
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

	/** Constructor */
	Likelihood(const std::vector<Peptide>&, const cParams&, const oParams&, const LikelihoodConstants&);
	/** Delete copy/assign constructors */
	Likelihood(const Likelihood&) = delete;
	Likelihood& operator=(const Likelihood&) = delete;
	/** Use a move constructor instead */
	Likelihood(Likelihood&& old) :
		peptideLikelihoods(std::move(old.peptideLikelihoods)),
		likelihood(std::move(old.likelihood)),
		c(std::move(old.c)),
		o(std::move(old.o)),
		constants(std::move(old.constants)),
		onupdate_c(std::move(old.onupdate_c)),
		onupdate_o(std::move(old.onupdate_o)) {
			// No need to re-link after a move (swap)
			// linkParamsAndPeptides(); // Just in case. TODO: check it is not needed and remove
		}

	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const Likelihood &Likelihood);
};
