#include "Likelihood.hpp"
//#include "prettyprint.hpp"
#include <Rcpp.h>

using namespace Rcpp;
using std::string;
using std::vector;

void Likelihood::linkParamsAndPeptides(cParams &aC, oParams &anO) {
	// Clear onupdate_c and onupdate_o
	onupdate_c.clear();
	onupdate_o.clear();
	// Set the right sizes
	onupdate_c.resize(aC.size());
	onupdate_o.resize(anO.size());
	for (size_t i = 0; i < onupdate_o.size(); ++i) {
		onupdate_o[i].resize(anO.size(i));
	}

	// Loop over the peptides
	for (PeptideLikelihood &peptideLikelihood: peptideLikelihoods) {
		// Get the relevant names
		string sampleName = peptideLikelihood.peptide.getSampleName();
		string refName = peptideLikelihood.peptide.getRefName();
		string pairName = peptideLikelihood.peptide.getPairName();

		// Link peptide to c
		size_t redundantCIdx = aC.getIndexOnRedundantC(pairName);
		peptideLikelihood.peptide.setC(aC.getPointerToRedundantC(redundantCIdx));
		// Get the non redundant c(s) for this peptide
		vector<size_t> cIndices = aC.getNonRedundantCFromRedundantC(redundantCIdx);
		for (size_t idx : cIndices) {
			// Link parameter to peptides to update
			onupdate_c.at(idx).push_back(&peptideLikelihood);
		}

		// Iterate over sites
		for (SiteSpecs& site: peptideLikelihood.peptide.getSiteSpecs()) {
			// Link peptide to o
			site.params.sample = anO.getPointerToO(sampleName, site.siteName);
			site.params.reference = anO.getPointerToO(refName, site.siteName);

			// Link o to peptide
			const size_t sampleIdx = anO.getFirstIndexOnO(sampleName);
			const size_t refIdx = anO.getFirstIndexOnO(refName);
			onupdate_o.at(sampleIdx).at(anO.getSecondIndexOnO(sampleIdx, site.siteName)).push_back(&peptideLikelihood);
			onupdate_o.at(refIdx).at(anO.getSecondIndexOnO(refIdx, site.siteName)).push_back(&peptideLikelihood);
		}
	}
}


Likelihood::Likelihood(const std::vector<Peptide> &peptides, cParams &aC, oParams &anO,
                       const Constants &someConstants) {

	peptideLikelihoods.reserve(peptides.size());
	for (auto &peptide: peptides) {
		peptideLikelihoods.push_back(PeptideLikelihood(peptide, someConstants.varianceModel, 0));
	}

	linkParamsAndPeptides(aC, anO);
	updateAll();
}


double Likelihood::updateAll() {
	// Reset likelihood
	likelihoodValue = 0;
	// Iterate over all peptides
	for (auto& pl: peptideLikelihoods) {
		// Ask the PeptideLikelihood to update itself and use it to increment the final value
		likelihoodValue += pl.update();
	}
	return likelihoodValue;
}

double Likelihood::updateC(const size_t i) {
	// Get the PeptideLikelihoods that are affected
	std::vector<PeptideLikelihood*> affected = onupdate_c.at(i);
	// Update them all
	for (PeptideLikelihood* pl: affected) {
		likelihoodValue -= pl->getValue(); // subtract previous likelihood
		likelihoodValue += pl->update(); // add new likelihood
	}
	// Returns the new total likelihood value
	return likelihoodValue;
}

double Likelihood::updateO(const size_t sample, const size_t site) {
	// Get the PeptideLikelihoods that are affected
	std::vector<PeptideLikelihood*> affected = onupdate_o.at(sample).at(site);
	// Update them all
	for (PeptideLikelihood* pl: affected) {
		likelihoodValue -= pl->getValue(); // subtract previous likelihood
		likelihoodValue += pl->update(); // add new likelihood
	}
	// Returns the new total likelihood value
	return likelihoodValue;
}


double Likelihood::changedC(const size_t i) {
	std::vector<PeptideLikelihood*> affected = onupdate_c.at(i);
	// Keep previous value
	double previousLikelihoodValue = likelihoodValue;
	// execute changed() on all peptides and update value as we go()
	for (PeptideLikelihood* pl: affected) {
		likelihoodValue += pl->changed();
	}
	// Return the change
	return likelihoodValue - previousLikelihoodValue;
}

double Likelihood::changedO(const size_t sample, const size_t site) {
	std::vector<PeptideLikelihood*> affected = onupdate_o.at(sample).at(site);
	// Keep previous value
	double previousLikelihoodValue = likelihoodValue;
	// execute changed() on all peptides and update value as we go()
	for (PeptideLikelihood* pl: affected) {
		likelihoodValue += pl->changed();
	}
	// Return the change
	return likelihoodValue - previousLikelihoodValue;
}

double Likelihood::temptativeChangedC(const size_t i) {
	std::vector<PeptideLikelihood*> affected = onupdate_c.at(i);
	// Start from current accepted likelihood
	temptativeLikelihoodValue = likelihoodValue;
	for (PeptideLikelihood* pl: affected) {
		// Update value with all affected peptideLikelihoods
		temptativeLikelihoodValue += pl->temptativeChanged();
	}
	// Return the total change
	return temptativeLikelihoodValue - likelihoodValue;


}
double Likelihood::temptativeChangedO(const size_t sample, const size_t site) {
	std::vector<PeptideLikelihood*> affected = onupdate_o.at(sample).at(site);
	// Start from current accepted likelihood
	temptativeLikelihoodValue = likelihoodValue;
	for (PeptideLikelihood* pl: affected) {
		// Update value with all affected peptideLikelihoods
		temptativeLikelihoodValue += pl->temptativeChanged();
	}
	// Return the total change
	return temptativeLikelihoodValue - likelihoodValue;
}

double Likelihood::acceptC(const size_t i) {
	// Accept the local temptative likelihood value
	likelihoodValue = temptativeLikelihoodValue;
	std::vector<PeptideLikelihood*> affected = onupdate_c.at(i);
	// Also do it for all affected peptideLikelihoods
	for (PeptideLikelihood* pl: affected) {
		pl->accept();
	}
	return likelihoodValue;
}
double Likelihood::acceptO(const size_t sample, const size_t site) {
	// Accept the local temptative likelihood value
	likelihoodValue = temptativeLikelihoodValue;
	std::vector<PeptideLikelihood*> affected = onupdate_o.at(sample).at(site);
	// Also do it for all affected peptideLikelihoods
	for (PeptideLikelihood* pl: affected) {
		pl->accept();
	}
	return likelihoodValue;
}
