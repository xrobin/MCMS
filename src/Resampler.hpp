#pragma once

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <iostream>
#include "Parameters.hpp"

class MoveSpec {
	static const boost::random::uniform_real_distribution<double> unif01;

	public:
	const double oldParam, newParam;
	const double logPreviousBias, logNewBias;

	MoveSpec(const double anOldParam, const double aNewParam, const double logPreviousBias, const double logNewBias):
		oldParam(anOldParam), newParam(aNewParam), logPreviousBias(logPreviousBias), logNewBias(logNewBias) {}
	bool accept(const double, std::mt19937_64&);

	/** Output */
	friend std::ostream& operator<< (std::ostream &out, const MoveSpec &aMoveSpec);
};

class Resampler {
	const Constants& constants;
	std::mt19937_64& rng;
	const boost::random::uniform_real_distribution<double> unif01;
	boost::random::normal_distribution<double> c_normal;
	boost::random::normal_distribution<double> o_normal;
	//const boost::random::normal_distribution<double>::param_type o_normal_params;

	public:
	Resampler(std::mt19937_64& aRNG, const Constants& someConstants):
		constants(someConstants), rng(aRNG),
		unif01(0.0, 1.0), c_normal(0, constants.c_sd), o_normal(0, constants.o_sd) {}

	/** Randomly choose if we're doing a prior move */
	bool priorMove() {
	return unif01(rng) > (1 - constants.priorMoveProportion);
}

    MoveSpec resampleC(const ParamSpecs& paramSpec, const double param);
    MoveSpec resampleO(const ParamSpecs& paramSpec, const double param);

    MoveSpec resampleCFromPrior(const ParamSpecs& paramSpec, const double param);
    MoveSpec resampleOFromPrior(const ParamSpecs& paramSpec, const double param);
    MoveSpec resampleCStandard(const ParamSpecs& paramSpec, const double param);
    MoveSpec resampleOStandard(const ParamSpecs& paramSpec, const double param);
};