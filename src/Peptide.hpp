#pragma once

#include <Rcpp.h>
//#include <vector>
#include <cstddef>
#include <boost/container/vector.hpp>

class Peptide {

	const size_t nSites;
	const double ratio, nEff, q;
	const lv_t siteActivity;
	const lv_ptrdiff_t offsetC, offsetOSample, offsetORef; // how far are the relevant parameters from their start?

	public:
	Peptide(const size_t aNSites,
			const double aRatio,
			const double aNEff,
			const double aQ,
			const lv_t& someSiteActivity,
			const lv_ptrdiff_t anOffsetC,
			const lv_ptrdiff_t anOffsetOSample,
			const lv_ptrdiff_t anOffsetORef):
			nSites(aNSites), ratio(aRatio), nEff(aNEff), q(aQ),
			siteActivity(someSiteActivity),
			offsetC(anOffsetC),
			offsetOSample(anOffsetOSample),
			offsetORef(anOffsetORef)
			{};

	double computeLikelihood(double *c, double[] *oSample, double[] *oRef);
};