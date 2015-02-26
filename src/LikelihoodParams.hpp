#pragma once

#include <vector>
#include <map>
#include <Rcpp.h>
#include <string>

//#include "typedefs.hpp"


///** Constant variables during the Monte Carlo sampling */
//class LikelihoodConstants {
//	public:
//	const Rcpp::NumericMatrix& sampleDependence;
//	explicit LikelihoodConstants(const Rcpp::NumericMatrix& aSampleDependenceMatrix):
//		sampleDependence(aSampleDependenceMatrix) {};
//};

///** Parameters for the Monte Carlo sampling */
//class LikelihoodParams {
//	//const LikelihoodConstants constants;
//	//cParams c;
//	//oParams o;
//	// Indicative names for c, redundantC and o
//	//CharacterVector cNames, redundantCNames, oNames;
//
//	public:
//	LikelihoodParams(const LikelihoodConstants& someLikelihoodConstants, const c_type& aC, const o_type& anO):
//		constants(someLikelihoodConstants), c(aC), redundantC(someLikelihoodConstants.sampleDependence.nrow()),
//		o(anO) {
//		updateRedundantC();
//	}
//
//	LikelihoodParams(const LikelihoodConstants& someLikelihoodConstants): constants(someLikelihoodConstants) {}
//
////	c_type getC() {
////		return(redundantC);
////	}
////	double getC(c_type::size_type i) {
////		return(redundantC[c]);
////	}
////
////	c_type::iterator getO(std::string key) {
////		return(o.at(key).begin());
////	}
//};


/** Storing the c parameters */
class cParams {
	typedef std::vector<double>::size_type size_t;
	const Rcpp::NumericMatrix sampleDependence;
	std::vector<double> c;
	std::vector<double> redundantC;
	std::unordered_map<std::string, size_t> cNames; // so from a name we know which c to update

	public:
	cParams() {} // Default constructor. TODO: remove?
	cParams(std::vector<double> someC, std::vector<std::string> someNames): c(someC) {
		if (someC.size() != someNames.size()) {
			throw std::invalid_argument( "someC and someNames sizes don't match" );
		}
		for (size_t i = 0; i < someNames.size(); ++i) {
			cNames.emplace(someNames[i], i);
		}
		updateRedundantC();
	}
	/** Updates the value of c at element i or key. */
	void updateC(size_t i, double newC) {
		c.at(i) = newC;
		updateRedundantC();
	}
	void updateC(std::string key, double newC) {
		updateC(cNames.at(key), newC);
	}

	/** Get a pointer to the redundant C at position i or key*/
//	*double getPointerToC() {
//		&c;
//	}
	double * getPointerToC(size_t i) {
		return &(c.at(i));
	}
	double * getPointerToC(std::string key) {
		return &(c.at(cNames.at(key)));
	}

	private:
	void updateRedundantC() {
		for (int row = 0; row < sampleDependence.nrow(); ++row) {
			redundantC[row] = 0;
			for (int col = 0; col < sampleDependence.ncol(); ++col) {
				int elem = sampleDependence(row, col);
				if (elem != 0) {
					redundantC[row] += elem * c[col];
				}
			}
		}
	}
};


/** Storing the c parameters */
class oParams {
	typedef std::vector<double>::size_type size_t;
	std::vector<std::vector<double>> o;
	std::unordered_map<std::string, size_t> sampleNames; // so from a name we know which sample number it is
	std::vector<std::unordered_map<std::string, size_t>> siteNames; // so from sample number and name we know which site it is

	public:
	typedef std::map<std::string, std::map<std::string, double>> o_type;
	/** Updates the value of o */
	void updateO(size_t sample, size_t site, double newO) {
		o.at(sample).at(site) = newO;
	}
	void updateO(size_t sample, std::string site, double newO) {
		updateO(sample, siteNames.at(sample).at(site), newO);
	}
	void updateO(std::string sample, std::string site, double newO) {
		updateO(sampleNames.at(sample), site, newO);
	}

	/** Get a pointer to the relevant o */
	double * getPointerToO(size_t sample, size_t site) {
		return &(o.at(sample).at(site));
	}
	double * getPointerToO(size_t sample, std::string site) {
		return getPointerToO(sample, siteNames.at(sample).at(site));
	}
	double * getPointerToO(std::string sample, std::string site) {
		return getPointerToO(sampleNames.at(sample), site);
	}

	oParams() {} // Default constructor. TODO: remove?
	oParams(o_type anOMap) {
		size_t sampleNumber = 0;
		for (auto& samplePair: anOMap) {
		//for(auto sample_iterator = anOMap.begin(); sample_iterator != anOMap.end(); ++sample_iterator) {
			std::string sampleName = samplePair.first;
			std::map<std::string, double> sampleSites = samplePair.second;
			std::vector<double> sampleOs;
			sampleOs.reserve(sampleSites.size());
			std::unordered_map<std::string, size_t> sampleSiteNames;
			sampleSiteNames.reserve(sampleSites.size());
			// Loop over sites
			size_t siteNumber = 0;
			for (auto& sitesPair: sampleSites) {
			//for(auto site_iterator = sampleSites.begin() ; site_iterator != sampleSites.end(); ++site_iterator) {
				std::string siteName = sitesPair.first;
				double siteO = sitesPair.second;
				sampleSiteNames.emplace(siteName, siteNumber);
				sampleOs.push_back(siteO);
				++siteNumber;
			}
			siteNames.push_back(sampleSiteNames);
			o.push_back(sampleOs);
			sampleNames[sampleName] = sampleNumber;
			++sampleNumber;
		}
    // iterator->first = key
    // iterator->second = value
    // Repeat if you also want to iterate through the second map.
	}
};
