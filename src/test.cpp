#include "conversions.hpp"
#include "LikelihoodParams.hpp"
#include <Rcpp.h>
#include "RcppConversions.hpp"
#include <string>
#include <vector>
#include "VarianceModel.hpp"

using Rcpp::as;
using Rcpp::S4;
using Rcpp::List;
using std::string;

//// For more on using Rcpp click the Help button on the editor toolbar
////' A Function to test the redundantC stuff
////' @param aProtein
////' @param aModel
////' @useDynLib MCMS
////' @export
//// [[Rcpp::export]]
//Rcpp::NumericVector run_MCMC(const Rcpp::S4& aProtein, const Rcpp::S4& aModel) {
//	string proteinClass = string(as<CharacterVector>(aProtein.attr("class"))[0]);
//	if (proteinClass != "Protein") {
//		stop(string("aProtein of class Protein expected, ") + proteinClass + " received");
//	}
//	string modelClass = string(as<CharacterVector>(aModel.attr("class"))[0]);
//	if (modelClass != "Peptides") {
//		stop(string("aModel of class Peptides expected, ") + modelClass + " received");
//	}
//
//	NumericMatrix sampleDependency = aProtein.slot("sample.dependency");
//	//LikelihoodConstants lc(sampleDependency);
//
//	NumericVector aC = aModel.slot("c");
//	oParams::o_type anO = convertListToO(aModel.slot("o"));
//	//LikelihoodParams lp(lc, aC, anO);
//
//	//std::vector<Peptides> peptides = convertListToPeptidesVector(aProtein.slot("data"), );
//	//Rcout << &peptides << "\n";
//
//
//	NumericVector ret(sampleDependency.nrow());
//	//for (int i = 0; i < ret.size(); ++i) {
//	//	ret[i] = *(lp.getC() + i);
//	//}
//	return(ret);
//}
//
//' @useDynLib MCMS
//' @title a test C++ function
//' @examples
//' data(ENSTest)
//' ENSTestProtein <- Protein(ENSTest)
//' ENSTestModel <- PeptidesModel(ENSTestProtein)
//' # testCpp(slot(ENSTestModel, "c"), slot(ENSTestProtein, "sample.dependency"))
//' testCpp(ENSTestProtein, ENSTestModel)
//' @export
// [[Rcpp::export]]
void testCpp(const S4& aModel, const List& aVarianceModelAsList) {
	VarianceModel aVarianceModel = as<VarianceModel>(aVarianceModelAsList);

//	oParams::o_type anOMap = convertListToO(oList);
//	oParams anO(anOMap);
//	double * ptr = anO.getPointerToO("P5", "Y_39");
//	Rcout << ptr << " = " << *ptr << "\n";
//	anO.updateO("P5", "Y_39", 1.01);
//	Rcout << ptr << " = " << *ptr << "\n";

//	Rcout << "converting c to map\n";
//	cParams::c_type aCMap = convertVectorToCMap(oVector);
//	Rcout << "Making real c\n";
//	cParams aC(aCMap, deps);
//	Rcout << "Getting a pointer to c\n";
//	double * ptr = aC.getPointerToC("P1_480");
//	Rcout << ptr << " = " << *ptr << "\n";
//	aC.updateC("P1_480", 1.01);
//	Rcout << ptr << " = " << *ptr << "\n";

	Likelihood l = convertS4ToLikelihood(aModel, aVarianceModel);
	//Rcpp::Rcout << l;
}
