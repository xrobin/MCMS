#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include "RcppHelpers.hpp" // for colnames
#include <string>
#include <vector>

using Rcpp::as;
using std::vector;
using std::string;

void cParams::updateRedundantC(){
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

oParams::oParams(const o_type &anOMap) {
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
}

cParams::cParams(const c_type &aCMap, const Rcpp::NumericMatrix &aSampleDependenceMatrix):
		sampleDependence(aSampleDependenceMatrix),
		c(), redundantC(aSampleDependenceMatrix.nrow()),
		cNames(), redundantCNames(),
		redundantCToC(aSampleDependenceMatrix.nrow()) {
	// Populate the main c and cNames
	size_t cNumber = 0;
	for (auto& samplePairPar: aCMap) {
		c.push_back(samplePairPar.second);
		cNames[samplePairPar.first] = cNumber;
		++cNumber;
	}
	// populate redundantCNames from the matrix
	vector<string> redundantNames = as<vector<string>>(Rcpp::rownames(aSampleDependenceMatrix));
	for (size_t redundantCNumber = 0; redundantCNumber < redundantNames.size(); ++redundantCNumber) {
		redundantCNames[redundantNames[redundantCNumber]] = redundantCNumber;
	}
	// Populate the redundantCToC indices

	vector<string> nonRedundantNames = as<vector<string>>(Rcpp::colnames(aSampleDependenceMatrix));
	int ncol = aSampleDependenceMatrix.ncol();
	int nrow = aSampleDependenceMatrix.nrow();
	// Iterate over the columns first
	for (int j = 0; j < ncol; ++j) {
		// Get the index on c for the column
		size_t colIdx = getIndexOnC(nonRedundantNames[j]);
		// Then go through the rows
		for (int i = 0; i < nrow; ++i) {
			if (aSampleDependenceMatrix(i, j) != 0) {
				// and then find the rows that are != 0
				redundantCToC.at(i).push_back(static_cast<size_t>(j));
			}
		}
	}

	updateRedundantC();
}