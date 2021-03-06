#include <cmath>
#include "Peptide.h"
#include <stdexcept>
#include <string>
#include "VarianceModel.h"


double Peptide::calcRatio() const {
	double sitesSum = 0;
	for (const SiteSpecs& site: siteSpecs) {
		double o_sample = *(site.params.sample);
		double o_ref = *(site.params.reference);
		sitesSum += log((1 - o_sample) / (1 - o_ref));
		if (site.siteActivity) {
			sitesSum += log(o_sample / (1 - o_sample)) - log(o_ref / (1 - o_ref));
		}
	}
	return *c + sitesSum;
}

double Peptide::computeLikelihood(const VarianceModel &aVarianceModel) const {
	double predictedRatio = calcRatio();
	double gamma = aVarianceModel.calcGamma(predictedRatio, nEff, q);
	double nuTilde = aVarianceModel.calcNuTilde(nEff);
	return - nuTilde * log1p(gamma * std::pow(ratio - predictedRatio, 2)) + std::log(gamma) / 2;
}