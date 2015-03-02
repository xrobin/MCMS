#include "Likelihood.hpp"
#include <Rcpp.h>

using namespace Rcpp;
using std::string;
using std::vector;

void Likelihood::computeLikelihood() {
	// Reset likelihood
	likelihood = 0;
	// Iterate over all peptides
	for (auto& pl: peptideLikelihoods) {
		// store individual likelihood in the peptideLikelihood
		pl.likelihood = pl.peptide.computeLikelihood();
		// and increment the global likelihood
		likelihood += pl.likelihood;
	}
}

void Likelihood::linkParamsAndPeptides() {
	// Clear onupdate_c and onupdate_o
	onupdate_c.clear();
	onupdate_o.clear();
	// Set the right sizes
	onupdate_c.resize(c.size());
	onupdate_o.resize(o.size());
	for (size_t i = 0; i < onupdate_o.size(); ++i) {
		onupdate_o[i].resize(o.size(i));
	}

	// Loop over the peptides
	for (PeptideLikelihood &peptideLikelihood: peptideLikelihoods) {
		// Get the relevant names
		string sampleName = peptideLikelihood.peptide.getSampleName();
		string refName = peptideLikelihood.peptide.getRefName();
		string pairName = peptideLikelihood.peptide.getPairName();

		// Link peptide to c
		size_t redundantCIdx = c.getIndexOnRedundantC(pairName);
		peptideLikelihood.peptide.setC(c.getPointerToRedundantC(redundantCIdx));
		// Get the non redundant c(s) for this peptide
		vector<size_t> cIndices = c.getNonRedundantCFromRedundantC(redundantCIdx);
		for (size_t idx : cIndices) {
			// Link parameter to peptides to update
			onupdate_c.at(idx).push_back(&peptideLikelihood);
		}

		// Iterate over sites
		for (siteSpec& site: peptideLikelihood.peptide.getSiteSpecs()) {
			// Link peptide to o
			site.params.sample = o.getPointerToO(sampleName, site.siteName);
			site.params.reference = o.getPointerToO(refName, site.siteName);

			// Link o to peptide
			const size_t sampleIdx = o.getFirstIndexOnO(sampleName);
			const size_t refIdx = o.getFirstIndexOnO(refName);
			onupdate_o.at(sampleIdx).at(o.getSecondIndexOnO(sampleIdx, site.siteName)).push_back(&peptideLikelihood);
			onupdate_o.at(refIdx).at(o.getSecondIndexOnO(refIdx, site.siteName)).push_back(&peptideLikelihood);
		}
	}
	c.prettyprint();
	std::cout << "onupdate_c = " << onupdate_c << std::endl;
	o.prettyprint();
	std::cout << "onupdate_o = " << onupdate_o << std::endl;
}


Likelihood::Likelihood(const std::vector<Peptide> &peptides, const cParams &aC, const oParams &anO):
	c(aC), o(anO) {
	peptideLikelihoods.reserve(peptides.size());
	for (auto &peptide: peptides) {
		peptideLikelihoods.push_back(PeptideLikelihood(peptide));
	}

	linkParamsAndPeptides();
	computeLikelihood();
}