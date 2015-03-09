
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include "Parameters.hpp"
#include "Resample.hpp"

double calcSigma(const BetaPrior& prior, const double o, const double o_sd, const double k_scale) {
	double k = prior.pdf(o) * (1 / (2 * (1 - o)) - 1 / (2 * o)) * k_scale;
	return 1 / (std::abs(k) + 1 / o_sd);
}

static const boost::random::uniform_int_distribution<size_t> MoveSpec::unif01(0, 1);

inline MoveSpec Resampler::resampleC(const ParamSpecs& paramSpec, const double param) {
	if(priorMove()) {
		resampleCFromPrior(paramSpec, param);
	}
	else {
		resampleCStandard(paramSpec, param);
	}
}


inline MoveSpec Resampler::resampleO(const ParamSpecs& paramSpec, const double param) {
	if(priorMove()) {
		resampleOFromPrior(paramSpec, param);
	}
	else {
		resampleOStandard(paramSpec, param);
	}
}

inline MoveSpec resampleCFromPrior(const ParamSpecs& paramSpec, const double oldC) {
	double logPreviousBias = std::log(paramSpec.prior.pdf(oldC));
	double newC = paramSpec.prior.pdf(oldC);
	double logNewBias = std::log(paramSpec.prior.pdf(newC));
	return MoveSpec(oldC, newC, logPreviousBias, logNewBias);
}
inline MoveSpec resampleOFromPrior(const ParamSpecs& paramSpec, const double oldO) {
	double logPreviousBias = std::log(paramSpec.prior.pdf(oldO));
	double newO = paramSpec.prior.pdf(oldO);
	double logNewBias = std::log(paramSpec.prior.pdf(newO));
	return MoveSpec(oldO, newO, logPreviousBias, logNewBias);
}

inline MoveSpec resampleCStandard(const ParamSpecs& paramSpec, const double oldC) {
	double newC = oldC + c_normal(rng);
	return MoveSpec(oldC, newC, 0, 0); // there is no change in the bias when we resample a c, still have a normal(0, c_sd)
}

inline MoveSpec resampleOStandard(const ParamSpecs& paramSpec, const double oldO) {
	// Compute the sd for normal
	double oldSigma = calcSigma(paramSpec.prior, oldO, constants.o_sd, constants.k_scale)
	// Resample a new O
	boost::random::normal_distribution<double>::param_type local_params(0, oldSigma);
	double newO =  o + o_normal(rng, local_params);
	// Calculate the biases
	// Previous
	boost::math::normal_distribution<double> oldNormal(oldO, oldSigma);
	double logPreviousBias = std::log(pdf(oldNormal, newO));
	// New
	boost::math::normal_distribution<double> newNormal(newO, newSigma);
	double logNewBias = std::log(pdf(newNormal, oldO));

	return MoveSpec(oldO, newOlogPreviousBias, logNewBias);
}

bool MoveSpec::accept(const double lpChange, std::mt19937_64& rng) {
	double nominator = likelihoodChange + logNewBias;
	double denominator = logPreviousBias;

	return (nominator > denominator || unif01(rng) < std::exp(nominator - denominator)));
}