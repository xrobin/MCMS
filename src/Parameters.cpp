#include "Parameters.h"
//#include "prettyprint.h"
#include <Rcpp.h>
#include "RcppHelpers.h" // for colnames
#include <string>
#include <vector>

using Rcpp::as;
using std::vector;
using std::string;

void cParams::updateRedundantC() {
	for (size_t i = 0; i < dependencyPairs.size(); ++i) { // over redundantC
		updateRedundantC(i);
	}
}

inline void cParams::updateRedundantC(const size_t i) {
	const auto& cVec = dependencyPairs[i];
	double newRedundantC = 0;
	for (size_t j = 0; j < cVec.size(); ++j) { // over c
		const DependencyPair &pair = cVec[j];
		newRedundantC += pair.mult * c[pair.i];
	}
	redundantC[i] = newRedundantC;
}

void cParams::updateRedundantC(const std::vector<size_t> &is) {
	for (size_t i = 0; i < is.size(); ++i) { // over redundantC
		updateRedundantC(is[i]);
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
		my_universal_unordered_map sampleSiteNames;
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
		dependencyPairs(),
		c(), redundantC(aSampleDependenceMatrix.nrow()),
		cNames(), redundantCNames(),
		redundantCToC(aSampleDependenceMatrix.nrow()), cToRedundantC() {
	// Info from the sample Dependency Matrix
	int ncol = aSampleDependenceMatrix.ncol();
	int nrow = aSampleDependenceMatrix.nrow();

	// Populate the main c and cNames
	size_t cNumber = 0;
	for (auto& samplePairPar: aCMap) {
		c.push_back(samplePairPar.second);
		cNames[samplePairPar.first] = cNumber;
		++cNumber;
	}
	cToRedundantC.resize(c.size());

	// populate redundantCNames from the matrix
	vector<string> redundantNames = as<vector<string>>(Rcpp::rownames(aSampleDependenceMatrix));
	vector<string> nonRedundantNames = as<vector<string>>(Rcpp::colnames(aSampleDependenceMatrix));
	for (size_t redundantCNumber = 0; redundantCNumber < redundantNames.size(); ++redundantCNumber) {
		redundantCNames[redundantNames[redundantCNumber]] = redundantCNumber;
	}

	// Populate dependencyPairs from matrix
	// Loop over redundantC (rows)
	for (int i = 0; i < nrow; ++i) {
		vector<DependencyPair> responsibleCs; // The c's impacting the redundantC i
		size_t rowIdx = getIndexOnRedundantC(redundantNames[i]);
		for (int j = 0; j < ncol; ++j) {
			double aMult = aSampleDependenceMatrix(i, j);
			// Where is column j in c?
			size_t colIdx = getIndexOnC(nonRedundantNames[j]);
			if (aMult != 0) {
				responsibleCs.push_back(DependencyPair(aMult, static_cast<size_t>(colIdx)));
			}
		}
		dependencyPairs.push_back(responsibleCs);
	}

	// Populate the redundantCToC indices
	// Iterate over the columns first
	for (int j = 0; j < ncol; ++j) {
		// Get the index on c for the column
		size_t colIdx = getIndexOnC(nonRedundantNames[j]);
		// Then go through the rows
		for (int i = 0; i < nrow; ++i) {
			if (aSampleDependenceMatrix(i, j) != 0) {
				// and then find the rows that are != 0
				redundantCToC.at(i).push_back(static_cast<size_t>(colIdx));
				cToRedundantC.at(colIdx).push_back(static_cast<size_t>(i));
			}
		}
	}

	updateRedundantC();
}

string cParams::getName(const size_t i) const {
	for (const std::pair<string, size_t>& name: cNames) {
		if (name.second == i) {
			return name.first;
		}
	}
	throw std::out_of_range("Index not found");
}

string cParams::getRedundantName(const size_t i) const {
	for (const std::pair<string, size_t>& name: redundantCNames) {
		if (name.second == i) {
			return name.first;
		}
	}
	throw std::out_of_range("Index not found");
}

string oParams::getSampleName(const size_t sample) const {
	for (const std::pair<string, size_t>& name: sampleNames) {
		if (name.second == sample) {
			return name.first;
		}
	}
	throw std::out_of_range("Index not found");
}


string oParams::getSiteName(const size_t sample, const size_t site) const {
	const my_universal_unordered_map& sampleSiteNames = siteNames[sample];
	for (const std::pair<string, size_t>& name: sampleSiteNames) {
		if (name.second == site) {
			return name.first;
		}
	}
	throw std::out_of_range("Index not found");
}
