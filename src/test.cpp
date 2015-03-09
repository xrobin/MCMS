#include "conversions.hpp"
#include "Parameters.hpp"
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
//	//Constants lc(sampleDependency);
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
//' testCpp(ENSTestModel, var.model,
//' 	scale = 1, shape1 = .5, shape2 = .5,
//' 	prior_move_proportion = .02, c_sd = 0.05, o_sd = 0.05, o_k_scale = 1/100)
//' @export
// [[Rcpp::export]]
void testCpp(const S4& aModel, const List& aVarianceModelAsList,
	const double scale, const double shape1, const double shape2,
	const double prior_move_proportion, const double c_sd, const double o_sd, const double o_k_scale) {
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

	using Rcpp::Rcout;
	Rcout << "In testCpp" << "\n";


	// Likelihood l = convertS4ToLikelihood(aModel, aVarianceModel, scale, shape1, shape2, prior_move_proportion, c_sd, o_sd, o_k_scale);
	MonteCarlo m = convertS4ToMonteCarlo(aModel, aVarianceModel, scale, shape1, shape2, prior_move_proportion, c_sd, o_sd, o_k_scale);

	Rcout << "Got the MC object" << "\n";

	m.iterate(10);

	Rcout << "Done!" << "\n";
	//Rcpp::Rcout << l;
}
