#pragma once

#include  <boost/random/normal_distribution.hpp>
#include "Parameters.hpp"

class Resampler {
	const Constants& constants;
	std::mt19937_64& rng;
	const boost::random::uniform_int_distribution<size_t> unif01;
	const boost::random::normal_distribution<double> c_normal;
	const boost::random::normal_distribution<double> o_normal;
	//const boost::random::normal_distribution<double>::param_type o_normal_params;

	Resampler(std::mt19937_64& aRNG, const Constants& someConstants):
		constants(someConstants), rng(aRNG),
		unif01(0, 1), c_normal(0, constants.c_sd), o_normal(0, constants.o_sd) {}

	public:
	/** Randomly choose if we're doing a prior move */
	bool priorMove() {
		return unif01(rng) > (1 - constants.priorMoveProportion);
	}

    double resampleC(const ParamSpecs& paramSpec, const double param);
    double resampleO(const ParamSpecs& paramSpec, const double param);

    double resampleCFromPrior(const ParamSpecs& paramSpec, const double param);
    double resampleOFromPrior(const ParamSpecs& paramSpec, const double param);
    double resampleCStandard(const ParamSpecs& paramSpec, const double param);
    double resampleOStandard(const ParamSpecs& paramSpec, const double param);
};