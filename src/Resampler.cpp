
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include "Parameters.hpp"
#include "Priors.hpp"
#include "Resampler.hpp"

double calcSigma(const GenericPrior& prior, const double o, const double o_sd, const double k_scale) {
	double k = prior.pdf(o) * (1 / (2 * (1 - o)) - 1 / (2 * o)) * k_scale/o_sd;
	return 1 / (std::abs(k) + 1 / o_sd);
}

const boost::random::uniform_real_distribution<double> MoveSpec::unif01(0.0, 1.0); // static initialization

MoveSpec Resampler::resampleC(const ParamSpecs& paramSpec, const double param) {
	if(priorMove()) {
		return resampleCFromPrior(paramSpec, param);
	}
	else {
		return resampleCStandard(paramSpec, param);
	}
}


MoveSpec Resampler::resampleO(const ParamSpecs& paramSpec, const double param) {
	if(priorMove()) {
		return resampleOFromPrior(paramSpec, param);
	}
	else {
		return resampleOStandard(paramSpec, param);
	}
}

MoveSpec Resampler::resampleCFromPrior(const ParamSpecs& paramSpec, const double oldC) {
	double newC = paramSpec.prior.sample(rng);
	double logPreviousBias = std::log(paramSpec.prior.pdf(newC));
	double logNewBias = std::log(paramSpec.prior.pdf(oldC));
	return MoveSpec(oldC, newC, logPreviousBias, logNewBias);
}

MoveSpec Resampler::resampleOFromPrior(const ParamSpecs& paramSpec, const double oldO) {
	double newO = paramSpec.prior.sample(rng);
	if (std::isnan(newO)) { // The resampled
		return (MoveSpec(false));
	}
	double logPreviousBias = std::log(paramSpec.prior.pdf(newO));
	double logNewBias = std::log(paramSpec.prior.pdf(oldO));
	return MoveSpec(oldO, newO, logPreviousBias, logNewBias);
}

MoveSpec Resampler::resampleCStandard(const ParamSpecs& paramSpec, const double oldC) {
	double newC = oldC + (c_normal(rng) / paramSpec.sd_corr);
	return MoveSpec(oldC, newC, 0, 0); // there is no change in the bias when we resample a c, still have a normal(0, c_sd)
}

MoveSpec Resampler::resampleOStandard(const ParamSpecs& paramSpec, const double oldO) {
	// Compute the sd for normal
	double oldSigma = calcSigma(paramSpec.prior, oldO, constants.o_sd / paramSpec.sd_corr, constants.o_k_scale);
	//Rcpp::Rcout << "oldSigma = " << oldSigma << std::endl;
	//Rcpp::Rcout << "paramSpec.sd_corr = " << paramSpec.sd_corr << std::endl;
	// Resample a new O
	boost::random::normal_distribution<double>::param_type local_params(0, oldSigma);
	double newO =  oldO + o_normal(rng, local_params);
	if (! paramSpec.prior.isValid(newO)) {
		return (MoveSpec(false));
	}
	double newSigma = calcSigma(paramSpec.prior, newO, constants.o_sd / paramSpec.sd_corr, constants.o_k_scale);
	//Rcpp::Rcout << "newSigma = " << newSigma << std::endl;
	// Calculate the biases
	// Previous
	boost::math::normal_distribution<double> oldNormal(oldO, oldSigma);
	double logPreviousBias = std::log(pdf(oldNormal, newO));
	// New
	boost::math::normal_distribution<double> newNormal(newO, newSigma);
	double logNewBias = std::log(pdf(newNormal, oldO));

	return MoveSpec(oldO, newO, logPreviousBias, logNewBias);
}

bool MoveSpec::accept(const double lpChange, std::mt19937_64& rng) {
	double nominator = lpChange + logNewBias;
	double denominator = logPreviousBias;

	return (nominator > denominator || unif01(rng) < std::exp(nominator - denominator));
}