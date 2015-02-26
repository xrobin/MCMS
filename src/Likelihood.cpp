#include "Likelihood.hpp"
#include <Rcpp.h>

using namespace Rcpp;
using std::string;
using std::vector;

void Likelihood::computeLikelihood() {
	// Reset likelihood
	likelihood = 0;
	// Iterate over all peptides
	for (auto& pl: peptideLikelihoods) {
		// store individual likelihood in the peptideLikelihood
		pl.likelihood = pl.peptide.computeLikelihood();
		// and increment the global likelihood
		likelihood += pl.likelihood;
	}
}

Likelihood::Likelihood(const Rcpp::S4& aProtein, const Rcpp::S4& aModel) {
	string proteinClass = string(as<CharacterVector>(aProtein.attr("class"))[0]);
	if (proteinClass != "Protein") {
		stop(string("aProtein of class Protein expected, ") + proteinClass + " received");
	}
	string modelClass = string(as<CharacterVector>(aModel.attr("class"))[0]);
	if (modelClass != "Peptides") {
		stop(string("aModel of class Peptides expected, ") + modelClass + " received");
	}

	NumericMatrix sampleDependency = aProtein.slot("sample.dependency");
	//params = LikelihoodConstants(sampleDependency);

	// TODO: make it recompile!
	//NumericVector aC = aModel.slot("c");
	//c = cParams(as<vector<double>>(aC), aC.names());
	//oParams::o_type anOMap = convertListToO(aModel.slot("o"));
	//o = oParams(anOMap);
	//LikelihoodParams lp(lc, aC, anO);

	//vector<Peptide> peptides = convertListToPeptidesVector(aProtein.slot("data"), );
	//Rcout << &peptides << "\n";

}