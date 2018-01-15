#pragma once

#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <iostream>
#include "Parameters.h"
//#include <cmath>
#include <random>
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

//class LaplacePrior: public GenericPrior {
//	const boost::math::laplace_distribution<double> laplace;
//	const boost::random::uniform_real_distribution<double> unif01;
//
//	public:
//	LaplacePrior(const double aScale):
//		laplace(0, aScale), unif01(0.0, 1.0) {}
//	double pdf(double x) const override {
//		return boost::math::pdf(laplace, x);
//	}
//	double sample(std::mt19937_64& rng) const override {
//		double uniform = unif01(rng);
//		return quantile(laplace, uniform);
//	}
//	void print(std::ostream &out) const override;
//};

class NormalPrior: public GenericPrior {
	const boost::math::normal_distribution<double> normal;
	const boost::random::uniform_real_distribution<double> unif01;

	public:
	NormalPrior(const double aSd):
		normal(0, aSd), unif01(0.0, 1.0) {}
	double pdf(double x) const override {
		return boost::math::pdf(normal, x);
	}
	double sample(std::mt19937_64& rng) const override {
		double uniform = unif01(rng);
		return quantile(normal, uniform);
	}
	void print(std::ostream &out) const override;
};

template <typename PriorType>
class ParamPriors {
	double *param; // link directly into cParams/oParams
	PriorType& prior;
	double priorValue;
	double temptativePriorValue; // prior of a proposed changed
	
	double calcPrior() {
		return std::log(prior.pdf(*param));
	}

	public:
	/** Constructor */
	ParamPriors(double *aParam, PriorType& aPrior):
		param(aParam), prior(aPrior),
		priorValue(calcPrior()), temptativePriorValue(0) {};

	/** MC functions */
	double changed() { /** Called just after a parameter was changed */
		double previousPriorValue = priorValue;
		priorValue = calcPrior();
		return priorValue - previousPriorValue;
	}
	/** The parameter just changed temptatively, update temptativePriorValue */
	double temptativeChanged() {
		temptativePriorValue = calcPrior();
		return temptativePriorValue - priorValue;
	}
	double accept() {
		priorValue = temptativePriorValue;
		return priorValue;
	}

	double update() { // just update and return
		priorValue = calcPrior();
		return priorValue;
	}

	double getValue() const {
		return priorValue;
	}

	/** Output */
	template <typename T>
	friend std::ostream& operator<< (std::ostream&, const ParamPriors<T>&);
};

template <typename T>
std::ostream& operator<< (std::ostream &out, const ParamPriors<T> &someParamPriors) {
	out << "<ParamPriors@" << &someParamPriors;
	out << ": " << someParamPriors.param << " (" << *(someParamPriors.param) << "), " << someParamPriors.priorValue;
	out << someParamPriors.prior;
	out << ">";
	return out;
}


/** This class complements the likelihood in the Monte Carlo estimation and computes the prior
 * It takes a reference to c and o, and calculates their priors
 *
 */
class Prior {
	//cParams &c;
	//oParams &o;
	NormalPrior normal;
	//LaplacePrior laplace;
	BetaPrior beta;
	std::vector<ParamPriors<NormalPrior>> cPriors;
	std::vector<std::vector<ParamPriors<BetaPrior>>> oPriors;
	double priorTotal, temptativePriorTotal;

	static std::vector<ParamPriors<NormalPrior>> makeCPriors(cParams &aC, NormalPrior &aNormalPrior);
	static std::vector<std::vector<ParamPriors<BetaPrior>>> makeOPriors(oParams &anO, BetaPrior &aBetaPrior);

	public:
	double updateAll();
	double updateC(const size_t);
	double updateO(const size_t, const size_t);

	Prior(cParams &aC, oParams &anO,
		const NormalPrior &aNormalPrior, const BetaPrior &aBetaPrior):
			normal(aNormalPrior), beta(aBetaPrior),
			cPriors(makeCPriors(aC, normal)), oPriors(makeOPriors(anO, beta)) {
		updateAll();
	}
	Prior(cParams &aC, oParams &anO,
		const double aCPriorSd, const double anOPriorShape1, const double anOPriorShape2,
		const double anO_restrict):
			Prior(aC, anO, NormalPrior(aCPriorSd), BetaPrior(anOPriorShape1, anOPriorShape2, anO_restrict)) {}

	NormalPrior& getNormalPrior() {return normal;}
	BetaPrior& getBetaPrior() {return beta;}

	double getPriorTotal() {return priorTotal;}

	double changedO(const size_t sample, const size_t site);
	double temptativeChangedO(const size_t sample, const size_t site);
	double acceptO(const size_t sample, const size_t site);
	double changedC(const size_t i);
	double temptativeChangedC(const size_t i);
	double acceptC(const size_t i);

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const Prior&);
};