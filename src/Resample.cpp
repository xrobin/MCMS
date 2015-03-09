#include "Parameters.hpp"
#include "Resample.hpp"


inline double Resampler::resampleC(const ParamSpecs& paramSpec, const double param) {
	if(priorMove()) {
		resampleCFromPrior(paramSpec, param);
	}
	else {
		resampleCStandard(paramSpec, param);
	}
}


inline double Resampler::resampleO(const ParamSpecs& paramSpec, const double param) {
	if(priorMove()) {
		resampleOFromPrior(paramSpec, param);
	}
	else {
		resampleOStandard(paramSpec, param);
	}
}

double resampleCFromPrior(const ParamSpecs& paramSpec, const double param) {
	return paramSpec.prior.pdf(param);
}
double resampleOFromPrior(const ParamSpecs& paramSpec, const double param) {
	return paramSpec.prior.pdf(param);
}

double resampleCStandard(const ParamSpecs& paramSpec, const double param) {
	return param + c_normal(rng);
}

double resampleOStandard(const ParamSpecs& paramSpec, const double param) {
	double k = paramSpec.prior.pdf(param) * (1 / (2 * (1 - param)) - 1 / (2 * param)) * constants.k_scale;
	boost::random::normal_distribution<double>::param_type local_params(0, k);
	return o + o_normal(rng, local_params);
}