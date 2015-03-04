#pragma once

#include <map>
#include <vector>
// #include "prettyprint.hpp"
#include "Priors.hpp"
#include <Rcpp.h>
#include <string>
#include "VarianceModel.hpp"

//#include "typedefs.hpp"


/** Constant variables during the Monte Carlo sampling */
class LikelihoodConstants {
	public:
	const VarianceModel varianceModel;
	explicit LikelihoodConstants(const VarianceModel& aVarianceModel):
		varianceModel(aVarianceModel) {};
};

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
	class DependencyPair { // Subclass to store the index and sign of the dependency
		public:
		double mult;
		size_t i; // index on redundantC
		DependencyPair(const double aMult, const size_t anI) :
			mult(aMult), i(anI) {}
	};

	typedef std::vector<double>::size_type size_t;
	//const Rcpp::NumericMatrix sampleDependence;
	std::vector<std::vector<DependencyPair>> dependencyPairs;
	std::vector<double> c, redundantC;
	std::unordered_map<std::string, size_t> cNames, redundantCNames; // so from a name we know which c to update
	std::vector<std::vector<size_t>> redundantCToC; // Maps a redundantC to one or more c on which it depends
	std::vector<std::vector<size_t>> cToRedundantC; // Maps a c to one or more redundantC which depend on it
	NormalPrior normalPrior;

public:
	typedef std::map<std::string, double> c_type;

	cParams(const c_type &aCMap, const Rcpp::NumericMatrix &aSampleDependenceMatrix);

	/** Updates the value of c at element i or key. */
	void updateC(const size_t i, const double newC) {
		c.at(i) = newC;
		updateRedundantC( cToRedundantC.at(i) );
	}
	void updateC(const std::string& key, double newC) {
		updateC(cNames.at(key), newC);
	}

	/** Get a pointer to the C redundant C at position i or key*/
	double * getPointerToC(const size_t i) {
		return &(c.at(i));
	}
	double * getPointerToC(const std::string& key) {
		return &(c.at(cNames.at(key)));
	}
	double * getPointerToRedundantC(const size_t i) {
		return &(redundantC.at(i));
	}
	double * getPointerToRedundantC(const std::string& key) {
		return &(redundantC.at(redundantCNames.at(key)));
	}

	/** Discover index from key */
	size_t getIndexOnC(const std::string& key) const {
		return cNames.at(key);
	}
	size_t getIndexOnRedundantC(const std::string& key) const {
		return redundantCNames.at(key);
	}

	/** When I have a redundant c, which non redundant c(s) are involved? */
	std::vector<size_t> getNonRedundantCFromRedundantC(const size_t i) const {
		return redundantCToC.at(i);
	}

	/** How many parameters do we have? */
	size_t size() const {
		return c.size();
	}
	size_t redundantSize() const {
		return redundantC.size();
	}

	/** Calculate the normal prior we have on all the o parameters */
	double calcPrior();

//	void prettyprint() const {
//		std::cout << "c = " << c << std::endl;
//		std::cout << "cNames = " << cNames << std::endl;
//		std::cout << "redundantCNames = " << redundantCNames << std::endl;
//		std::cout << "redundantCToC = " << redundantCToC << std::endl;
//	}

	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const cParams &aCParams);

	private:
	void updateRedundantC();
	void updateRedundantC(const size_t i);
	void updateRedundantC(const std::vector<size_t> &i);
};


/** Storing the c parameters */
class oParams {
	typedef std::vector<double>::size_type size_t;
	std::vector<std::vector<double>> o;
	std::unordered_map<std::string, size_t> sampleNames; // so from a name we know which sample number it is
	std::vector<std::unordered_map<std::string, size_t>> siteNames; // so from sample number and name we know which site it is
	BetaPrior betaPrior;

public:
	/** Constructor from an o_type map*/
	typedef std::map<std::string, std::map<std::string, double>> o_type;
	oParams(const o_type &anOMap, const BetaPrior& aBetaPrior);
	oParams(const o_type &anOMap, const double aShape1, const double aShape2):
		oParams(anOMap, BetaPrior(aShape1, aShape2)) {};

	/** Updates the value of o */
	void updateO(const size_t sample, const size_t site, const double newO) {
		o.at(sample).at(site) = newO;
	}
	void updateO(const size_t sample, const std::string& site, double newO) {
		updateO(sample, siteNames.at(sample).at(site), newO);
	}
	void updateO(const std::string& sample, const std::string& site, double newO) {
		updateO(sampleNames.at(sample), site, newO);
	}

	/** Get a pointer to the relevant o */
	double * getPointerToO(const size_t sample, const size_t site) {
		return &(o.at(sample).at(site));
	}
	double * getPointerToO(const size_t sample, const std::string& site) {
		return getPointerToO(sample, siteNames.at(sample).at(site));
	}
	double * getPointerToO(const std::string& sample, const std::string& site) {
		return getPointerToO(sampleNames.at(sample), site);
	}

	/** How many parameters do we have? */
	size_t size() const {
		return o.size();
	}
	size_t size(const size_t sample) const {
		return o.at(sample).size();
	}
	size_t size(const std::string& sample) const {
		return size(sampleNames.at(sample));
	}

	/** Calculate the beta prior we have on all the o parameters */
	double calcPrior();

	/** Get indices */
	size_t getFirstIndexOnO(const std::string& sample) const {
		return sampleNames.at(sample);
	}
	size_t getSecondIndexOnO(const std::string& sample, const std::string& site) const {
		return getSecondIndexOnO(sampleNames.at(sample), site);
	}
	size_t getSecondIndexOnO(const size_t sample, const std::string& site) const {
		return siteNames.at(sample).at(site);
	}

//	void prettyprint() const {
//		std::cout << "o = " << o << std::endl;
//		std::cout << "sampleNames = " << sampleNames << std::endl;
//		std::cout << "siteNames = " << siteNames << std::endl;
//	}

	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const oParams &anOParams);
};
