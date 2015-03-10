#pragma once

#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <iostream>
#include "Parameters.hpp"
//#include <cmath>
#include <random>
//#include <boost/math/distributions/normal.hpp>
//#include <random>

class GenericPrior {
	public:
	virtual double pdf(double x) const = 0;
	virtual double sample(std::mt19937_64& rng) const = 0;
    virtual void print(std::ostream &out) const = 0;
    virtual bool isValid(double x) const { //valid by default...
    	return true;
    }
};

class BetaPrior: public GenericPrior {
	typedef boost::math::policies::policy<
		boost::math::policies::overflow_error<
			boost::math::policies::ignore_error>>
		ignore_error_policy;
	const boost::math::beta_distribution<double, ignore_error_policy> beta;
	const boost::random::uniform_real_distribution<double> unif01;
	const double o_restrict;

	public:
	BetaPrior(const double aShape1, const double aShape2,
		const double anO_restrict):
		beta(aShape1, aShape2), unif01(0.0, 1.0),
		o_restrict(anO_restrict) {}

	double pdf(double x) const override {
		if (!isValid(x)) {
			return 0;
		}
		return boost::math::pdf(beta, x);
	}
	double sample(std::mt19937_64& rng) const override  {
		double uniform = unif01(rng);
		double val = quantile(beta, uniform);
		if (!isValid(val)) {
			return std::numeric_limits<double>::signaling_NaN();
		}
		return val;
	}
	/** Returns true if x is in the range between lowerBound and upperBound, false otherwise */
	bool isValid(double x) const override {
		return x > o_restrict && x < (1 - o_restrict);
	}
	void print(std::ostream &out) const override;
};

class LaplacePrior: public GenericPrior {
	const boost::math::laplace_distribution<double> laplace;
	const boost::random::uniform_real_distribution<double> unif01;

	public:
	LaplacePrior(const double aScale):
		laplace(0, aScale), unif01(0.0, 1.0) {}
	double pdf(double x) const override {
		return boost::math::pdf(laplace, x);
	}
	double sample(std::mt19937_64& rng) const override {
		double uniform = unif01(rng);
		return quantile(laplace, uniform);
	}
	void print(std::ostream &out) const override;
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
	friend std::ostream& operator<< (std::ostream&, const ParamPriors<PriorType>&);
};

template <typename PriorType>
std::ostream& operator<< (std::ostream &out, const ParamPriors<PriorType> &someParamPriors) {
	out << "<ParamPriors@" << &someParamPriors;
	out << ": " << someParamPriors.param << " (" << *(someParamPriors.param) << "), " << someParamPriors.priorValue;
	out << someParamPriors.prior;
	out << ">";
}


/** This class complements the likelihood in the Monte Carlo estimation and computes the prior
 * It takes a reference to c and o, and calculates their priors
 *
 */
class Prior {
	//cParams &c;
	//oParams &o;
	LaplacePrior laplace;
	BetaPrior beta;
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
		const LaplacePrior &aLaplacePrior, const BetaPrior &aBetaPrior):
			laplace(aLaplacePrior), beta(aBetaPrior),
			cPriors(makeCPriors(aC, laplace)), oPriors(makeOPriors(anO, beta)) {
		updateAll();
	}
	Prior(cParams &aC, oParams &anO,
		const double aScale, const double aShape1, const double aShape2,
		const double anO_restrict):
			Prior(aC, anO, LaplacePrior(aScale), BetaPrior(aShape1, aShape2, anO_restrict)) {}

	LaplacePrior& getLaplacePrior() {return laplace;}
	BetaPrior& getBetaPrior() {return beta;}

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

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const Prior&);
};