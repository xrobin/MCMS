#pragma once

#include <vector>
#include <map>
#include <Rcpp.h>
#include <string>

#include "typedefs.hpp"


/** Constant variables during the Monte Carlo sampling */
class LikelihoodConstants {
	public:
	const Rcpp::NumericMatrix& sampleDependence;
	LikelihoodConstants(const Rcpp::NumericMatrix& aSampleDependenceMatrix):
		sampleDependence(aSampleDependenceMatrix) {};
};

/** Parameters for the Monte Carlo sampling */
class LikelihoodParams {
	const LikelihoodConstants& constants;
	c_type c, redundantC;
	o_type o;

	void updateRedundantC();

	public:
	LikelihoodParams(const LikelihoodConstants& someLikelihoodConstants, const c_type& aC, const o_type& anO):
		constants(someLikelihoodConstants), c(aC), redundantC(someLikelihoodConstants.sampleDependence.nrow()),
		o(anO) {
		updateRedundantC();
	}

	c_type::iterator getC() {
		return(redundantC.begin());
	}
	c_type::iterator getO(std::string key) {
		return(o.at(key).begin());
	}
};