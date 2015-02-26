#pragma once

//#include <Rcpp.h>
//#include <vector>
//#include <cstddef>
#include "typedefs.hpp"
#include <vector>

// Describes a pair of pointers to occupancy ratios for the sample and reference class.
class oPair {
	public:
	double *reference, *sample;
};

class Peptide {
	//const size_t nSites;
	const double ratio, nEff, q;
	const lv_type siteActivity;
	const std::string sampleName, refName, pairName;


	// Pointers to the parameters
	double *c;
	std::vector<oPair> oPairs;
	//o_value_type *oRef, *oSample;

	//double previousLikelihood;

	public:
	Peptide(const double aRatio,
			const double aNEff,
			const double aQ,
			const lv_type& someSiteActivity,
			const std::string &theSampleName,
			const std::string &theRefName,
			const std::string &thePairName
			):
			ratio(aRatio), nEff(aNEff), q(aQ),
			siteActivity(someSiteActivity),
			sampleName(theSampleName), refName(theRefName), pairName(thePairName)
			{};

	//double computeLikelihood(double *c, double[] *oSample, double[] *oRef);
	double computeLikelihood();

	/** Functions to set the pointers */
	void setC(double *targetC) {
		c = targetC;
	}

	void setO(std::vector<oPair>& newOPairs) {
		oPairs = newOPairs;
	}

//	void setORef(c_type *targetORef) {
//		oRef = targetORef;
//	}
//	void setOSample(c_type *targetOSample) {
//		oSample = targetOSample;
//	}

};