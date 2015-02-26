#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include <string>
#include "typedefs.hpp"
#include <vector>

using namespace Rcpp;
using std::map;
using std::string;
using std::vector;

oParams::o_type convertListToO(List anOList) {
	CharacterVector sampleNames = anOList.names();
	oParams::o_type anO;
	//anO.reserve(sampleNames.size());
	for (auto& sampleName: as<vector<string>>(sampleNames)) {
		map<string, double> sampleSitesMap;
		NumericVector sampleSites = anOList[sampleName];
		vector<string> siteNames = sampleSites.names();
		for (int i = 0; i < sampleSites.size(); ++i) {
			string siteName = siteNames[i];
			double siteO = sampleSites[i];
			sampleSitesMap[siteName] = siteO;
		}
		anO[sampleName] = sampleSitesMap;
	}
	return(anO);
}

//std::map<std::string, std::map<std::string, double>>

//vector<Peptides> convertListToPeptidesVector(List aDataList, ) {
//	vector<Peptides> peps;
//	peps.reserve(aDataList.size());
//	vector<double> ratios = aDataList["ratio"],
//	               qs     = aDataList["q"],
//	               nEffs  = aDataList["n"];
//	vector<string> pairs      = aDataList["pair"],
//	               samples    = aDataList["sample"],
//	               references = aDataList["reference"];
//
//	aDataList
//	Rcout << &peps << "\n";
//	Rcout << &peps.size() << "\n";
//	peps.reserve()
//	return(anO);
//}