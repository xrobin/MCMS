#pragma once

#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/laplace.hpp>
#include "LikelihoodParams.hpp"
//#include <cmath>
#include <random>
//#include <boost/math/distributions/normal.hpp>
//#include <random>

//class GenericPrior {
//	double pdf(double x) = 0;
//	double sample(std::mt19937_64& rng) = 0;
//};

class BetaPrior {
	typedef boost::math::policies::policy<
		boost::math::policies::overflow_error<
			boost::math::policies::ignore_error>>
		ignore_error_policy;
	boost::math::beta_distribution<double, ignore_error_policy> beta;
	std::uniform_real_distribution<double> unif01;

	public:
	BetaPrior(const double aShape1, const double aShape2):
		beta(aShape1, aShape2), unif01(0.0, 1.0) {}
	double pdf(double x) {
		return boost::math::pdf(beta, x);
	}
	double sample(std::mt19937_64& rng) {
		double uniform = unif01(rng);
		return quantile(beta, uniform);
	}
};

class LaplacePrior {
	boost::math::laplace_distribution<double> laplace;
	std::uniform_real_distribution<double> unif01;

	public:
	LaplacePrior(const double aScale):
		laplace(0, aScale), unif01(0.0, 1.0) {}
	double pdf(double x) {
		return boost::math::pdf(laplace, x);
	}
	double sample(std::mt19937_64& rng) {
		double uniform = unif01(rng);
		return quantile(laplace, uniform);
	}
};

template <typename PriorType>
class ParamPriors {
	double *param; // link directly into cParams/oParams
	PriorType& prior;
	double priorValue;
	double temptativePriorValue; // prior of a proposed changed

	public:
	/** Constructor */
	ParamPriors(double *aParam, PriorType& aPrior):
		param(aParam), prior(aPrior),
		priorValue(aPrior.pdf(*aParam)), temptativePriorValue(0) {};

	/** MC functions */
	double changed() { /** Called just after a parameter was changed */
		double previousPriorValue = priorValue;
		priorValue = prior.pdf(*param);
		return priorValue - previousPriorValue;
	}
	/** The parameter just changed temptatively, update temptativePriorValue */
	double temptativeChanged() {
		temptativePriorValue = prior.pdf(*param);
		return temptativePriorValue - priorValue;
	}
	void accept() {
		priorValue = temptativePriorValue;
	}

	double update() { // just update and return
		priorValue = prior.pdf(*param);
		return priorValue;
	}

	double getValue() const {
		return priorValue;
	}

	/** Output */
	//friend std::ostream& operator<< (std::ostream &out, const PeptideLikelihood &aPeptideLikelihood);
};


/** This class complements the likelihood in the Monte Carlo estimation and computes the prior
 * It takes a reference to c and o, and calculates their priors
 *
 */
class Prior {
	//cParams &c;
	//oParams &o;
	std::vector<ParamPriors<LaplacePrior>> cPriors;
	std::vector<std::vector<ParamPriors<BetaPrior>>> oPriors;
	double priorTotal, temptativePriorTotal;

	static std::vector<ParamPriors<LaplacePrior>> makeCPriors(cParams &aC, LaplacePrior &aLaplacePrior);
	static std::vector<std::vector<ParamPriors<BetaPrior>>> makeOPriors(oParams &anO, BetaPrior &aBetaPrior);

	public:
	double updateAll();
	double updateC(const size_t);
	double updateO(const size_t, const size_t);

	Prior(cParams &aC, oParams &anO,
		LaplacePrior &aLaplacePrior, BetaPrior &aBetaPrior):
			cPriors(makeCPriors(aC, aLaplacePrior)), oPriors(makeOPriors(anO, aBetaPrior)) {
		updateAll();
	}
	double changedO(const size_t sample, const size_t site) {
		return oPriors.at(sample).at(site).changed();
	}
	double temptativeChangedO(const size_t sample, const size_t site) {
		return oPriors.at(sample).at(site).temptativeChanged();
	}
	void acceptO(const size_t sample, const size_t site) {
		oPriors.at(sample).at(site).accept();
	}
	double changedC(const size_t i) {
		return cPriors.at(i).changed();
	}
	double temptativeChangedC(const size_t i) {
		return cPriors.at(i).temptativeChanged();
	}
	void acceptC(const size_t i) {
		cPriors.at(i).accept();
	}
};