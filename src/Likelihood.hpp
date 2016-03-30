#pragma once

#include <iostream> // cout
#include "Parameters.hpp"
#include <map>
#include "Peptide.hpp"
#include <Rcpp.h>
#include <vector>


class PeptideLikelihood {
	public:
	Peptide peptide;
	private:
	double likelihoodValue, temptativeLikelihoodValue;
	const VarianceModel& varianceModel;
	public:
	/** Constructors */
	/** With a given likelihood (avoids computing it before we've linked the parameters) */
	PeptideLikelihood(const Peptide& aPeptide, const VarianceModel& aVarianceModel, const double aLikelihood):
		peptide(aPeptide), likelihoodValue(aLikelihood),
		temptativeLikelihoodValue(0), varianceModel(aVarianceModel) {};
	PeptideLikelihood(const Peptide& aPeptide, const VarianceModel& aVarianceModel):
		PeptideLikelihood(aPeptide, aVarianceModel, aPeptide.computeLikelihood(aVarianceModel)) {};
	//PeptideLikelihood(const Peptide& aPeptide, const double aLikelihood):
	//	peptide(aPeptide), likelihood(aLikelihood) {};

	/** MC functions */
	double changed() { /** Called just after a parameter was changed */
		double previousLikelihoodValue = likelihoodValue;
		likelihoodValue = peptide.computeLikelihood(varianceModel);
		return likelihoodValue - previousLikelihoodValue;
	}
	/** The parameter just changed temptatively, update temptativePriorValue */
	double temptativeChanged() {
		temptativeLikelihoodValue = peptide.computeLikelihood(varianceModel);
		return temptativeLikelihoodValue - likelihoodValue;
	}
	double accept() {
		likelihoodValue = temptativeLikelihoodValue;
		return likelihoodValue;
	}

	double update() { // just update and return
		likelihoodValue = peptide.computeLikelihood(varianceModel);
		return likelihoodValue;
	}

	double getValue() const {
		return likelihoodValue;
	}

	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const PeptideLikelihood &aPeptideLikelihood);
};


class Likelihood {
	std::vector<PeptideLikelihood> peptideLikelihoods;
	//std::vector<Peptide> peptides;
	//std::vector<double> likelihoods;
	double likelihoodValue, temptativeLikelihoodValue;

	// Store likelihood parameters
	//cParams c;
	//oParams o;
	//const Constants constants;


	// mapping from param to peptides that must be updated
	// Maps an index on the c to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<PeptideLikelihood*>> onupdate_c;
	//   c           peptide indices
	// Maps an index on the o to a vector of indices on Peptides and likelihoods to be updated
	std::vector<std::vector<std::vector<PeptideLikelihood*>>> onupdate_o;
	//  sample     site       peptide indices

	/** Link the peptides and parameters */
	void linkParamsAndPeptides(cParams&, oParams&);
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
	Likelihood(const std::vector<Peptide>&, cParams&, oParams&, const Constants&);
	/** Delete copy/assign constructors */
	Likelihood(const Likelihood&) = delete;
	Likelihood& operator=(const Likelihood&) = delete;
	/** Use a move constructor instead */
	Likelihood(Likelihood&& old) :
		peptideLikelihoods(std::move(old.peptideLikelihoods)),
		likelihoodValue(std::move(old.likelihoodValue)),
		temptativeLikelihoodValue(std::move(old.temptativeLikelihoodValue)),
		onupdate_c(std::move(old.onupdate_c)),
		onupdate_o(std::move(old.onupdate_o)) {
			// No need to re-link after a move (swap)
			// linkParamsAndPeptides(); // Just in case. TODO: check it is not needed and remove
		}

	double getLikelihoodValue() {return likelihoodValue;}

	/** MC functions */

	/** Updates the likelihood (potentially only some peptides linked with the given parameter)
	 * and returns the new total likelihood */
	double updateAll();
	double updateC(const size_t);
	double updateO(const size_t, const size_t);

	/** After a parameter was changed, updates the likelihood and returns the change */
	double changedC(const size_t);
	double changedO(const size_t, const size_t);

	/** After a parameter was changed temptatively, updates the temptative likelihood values and returns the change
	 * Does not change the likelihood itself. */
	double temptativeChangedC(const size_t);
	double temptativeChangedO(const size_t, const size_t);

	/** After a tempative change, makes the change effective.
	 * Warning: no checking is done, make sure to pass the same indices
	 */
	double acceptC(const size_t);
	double acceptO(const size_t, const size_t);

//	double getC(const size_t i) {
//		return c.getC(i);
//	}
//
//
//	double getO(const size_t sample, const size_t site) {
//		return o.getO(sample, site);
//	}

//	/** Computes the whole likelihood and stores it */
//	double getLikelihood() {
//		return likelihoodValue;
//	}

	/** Returns the number of likelihoods that will be affected by the change in c or o */
	size_t onUpdateCSize(const size_t i) {
		return onupdate_c.at(i).size();
	}
	size_t onUpdateOSize(const size_t sample, const size_t site) {
		return onupdate_o.at(sample).at(site).size();
	}

	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const Likelihood &Likelihood);
};
