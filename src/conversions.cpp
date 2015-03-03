#include "conversions.hpp"
#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include "RcppHelpers.hpp" // for colnames
#include <string>
#include "typedefs.hpp"
#include <vector>

using Rcpp::as;
using Rcpp::DataFrame;
using Rcpp::CharacterVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::stop;
using std::map;
using std::string;
using std::vector;

oParams::o_type convertListToOMap(const List &anOList) {
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

cParams::c_type convertVectorToCMap(const NumericVector &aCVector) {
	vector<string> pairNames = as<vector<string>>(aCVector.names());
	cParams::c_type aCMap;
	for (int i = 0; i < aCVector.size(); ++i) {
		aCMap[pairNames[i]] = aCVector[i];
	}
	return(aCMap);
}


Likelihood convertS4ToLikelihood(const Rcpp::S4& aModel, const VarianceModel& aVarianceModel) {
	const string modelClass = string(as<CharacterVector>(aModel.attr("class"))[0]);
	if (modelClass != "Peptides") {
		stop(string("aModel of class Peptides expected, ") + modelClass + " received");
	}
	const Rcpp::S4 aProtein = as<Rcpp::S4>(aModel.slot("protein"));
	const string proteinClass = string(as<CharacterVector>(aProtein.attr("class"))[0]);
	if (proteinClass != "Protein") {
		stop(string("aProtein of class Protein expected, ") + proteinClass + " received");
	}

	// Build the parameters
	// First O
	const oParams::o_type anOMap = convertListToOMap(aModel.slot("o"));
	const oParams anO(anOMap);

	// Then C
	const cParams::c_type aCMap = convertVectorToCMap(aModel.slot("c"));
	const NumericMatrix sampleDependency = aProtein.slot("sample.dependency");
	const cParams aC(aCMap, sampleDependency);
//	NumericVector aC = aModel.slot("c");

	// Create the peptides
	std::vector<Peptide> peptides;
	// Extract them from the data
	const DataFrame data = as<DataFrame>(aProtein.slot("data"));
	const NumericMatrix sitesCoverage = as<NumericMatrix>(aProtein.slot("sites.coverage"));
	const NumericMatrix sitesActivation = as<NumericMatrix>(aProtein.slot("sites.activation"));
	// Extract the columns in data for faster access
	const NumericVector ratios = data["ratio"];
	const NumericVector qs = data["q"];
	const NumericVector nEffs = data["n"];
	const vector<string> sampleNames = data["sample"];
	const vector<string> referenceNames = data["reference"];
	const vector<string> pairNames = data["pair"];
	// Indexed over number of sites
	const vector<string> siteNames = as<vector<string>>(Rcpp::colnames(sitesCoverage));
	// Now loop over all this
	for (int i = 0; i < data.nrows(); ++i) {
		std::vector<SiteSpecs> siteSpecs;
		// Loop over the sites
		for (int j = 0; j < sitesCoverage.ncol(); ++j) {
			if (sitesCoverage(i, j) == 1) {
				bool newSiteActivity = sitesActivation(i, j);
				string newSiteName = siteNames[j];
				siteSpecs.push_back(SiteSpecs(newSiteActivity, newSiteName));
			}
			else {
				string newSiteName = siteNames[j];
			}
		}
		// Create the Peptide object
		Peptide newPeptide(ratios[i], nEffs[i], qs[i],
			sampleNames[i], referenceNames[i], pairNames[i],
			siteSpecs);
		peptides.push_back(newPeptide);
	}

	Likelihood l (peptides, aC, anO, LikelihoodConstants(aVarianceModel));
	//Rcpp::Rcout << l;
	return(l);
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