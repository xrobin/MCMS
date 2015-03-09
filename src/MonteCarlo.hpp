#pragma once

#include <boost/random/uniform_int_distribution.hpp>
#include "Helpers.hpp"
#include "Likelihood.hpp"
#include "Priors.hpp"
#include "Resampler.hpp"

class MonteCarlo {
	const Constants constants;
	cParams c;
	oParams o;
	Prior p;
	Likelihood l;
	const ParamSpecsVector paramSpecs; // a vector with all the parameters
	Resampler resampler;
	std::mt19937_64& rng; // The random numbers for the whole MC run

	static ParamSpecsVector makeParamSpecsVector(const cParams& c, const oParams& o, Prior& p);

	public:
	MonteCarlo(const std::vector<Peptide>& peptides,
		const cParams::c_type& aCMap, const oParams::o_type& anOMap,
		const Constants& someConstants,
		std::mt19937_64& anRNG):
		constants(someConstants),
		c(aCMap, constants.sampleDependenceMatrix),
		o(anOMap),
		p(c, o, constants.scale, constants.shape1, constants.shape2),
		l(peptides, c, o, constants),
		paramSpecs(makeParamSpecsVector(c, o, p)),
		resampler(anRNG, constants),
		rng(anRNG){}

	void iterate(unsigned long);
	void iterate();
};