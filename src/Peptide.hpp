#pragma once

//#include <Rcpp.h>
//#include <vector>
//#include <cstddef>
#include <string>
#include "typedefs.hpp"
#include <vector>

// Describes a pair of pointers to occupancy ratios for the sample and reference class.
class oPair {
	public:
	double *reference, *sample;
};

class siteSpec {
	public:
	const bool siteActivity;
	const std::string siteName;
	oPair params;
	siteSpec(const bool aSiteActivity, const std::string &aSiteName):
		siteActivity(aSiteActivity), siteName(aSiteName) {};
	siteSpec(const bool aSiteActivity, const std::string &aSiteName, const oPair &aParams):
		siteActivity(aSiteActivity), siteName(aSiteName), params(aParams) {};
};

class Peptide {
	//const size_t nSites;
	const double ratio, nEff, q;
	const std::string sampleName, refName, pairName;


	// Pointers to the parameters
	double *c;
	std::vector<siteSpec> siteSpecs;
	//o_value_type *oRef, *oSample;

	//double previousLikelihood;

	public:
	Peptide(const double aRatio,
			const double aNEff,
			const double aQ,
			const std::string &theSampleName,
			const std::string &theRefName,
			const std::string &thePairName,
			std::vector<siteSpec> &theSiteSpecs
			):
			ratio(aRatio), nEff(aNEff), q(aQ),
			sampleName(theSampleName), refName(theRefName), pairName(thePairName),
			c(nullptr), siteSpecs(theSiteSpecs)
			{};

	//double computeLikelihood(double *c, double[] *oSample, double[] *oRef);
	double computeLikelihood();

	/** Getters */
	std::string getSampleName() const {
		return(sampleName);
	}
	std::string getRefName() const {
		return(refName);
	}
	std::string getPairName() const {
		return(pairName);
	}
	/** Give access to the site specs by reference */
	std::vector<siteSpec>& getSiteSpecs() {
		return siteSpecs;
	}


	/** Functions to set the pointers */
	void setC(double * targetC) {
		c = targetC;
	}
//
//	void setO(std::vector<oPair>& newOPairs) {
//		oPairs = newOPairs;
//	}

//	void setORef(c_type *targetORef) {
//		oRef = targetORef;
//	}
//	void setOSample(c_type *targetOSample) {
//		oSample = targetOSample;
//	}

};