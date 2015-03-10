#pragma once

#include <boost/random/uniform_int_distribution.hpp>
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
	std::mt19937_64 rng; // The random numbers for the whole MC run
	Resampler resampler;

	static ParamSpecsVector makeParamSpecsVector(cParams& c, oParams& o, Prior& p);

	public:
	MonteCarlo(const std::vector<Peptide>& peptides,
		const cParams::c_type& aCMap, const oParams::o_type& anOMap,
		const Constants& someConstants,
		const std::mt19937_64& anRNG):
		constants(someConstants),
		c(aCMap, constants.sampleDependenceMatrix),
		o(anOMap),
		p(c, o, constants.scale, constants.shape1, constants.shape2),
		l(peptides, c, o, constants),
		paramSpecs(makeParamSpecsVector(c, o, p)),
		rng(anRNG),
		resampler(rng, constants){}
	/** Delete dangerous copy and copy assignment constructors */
	MonteCarlo(const MonteCarlo&) = delete;
	MonteCarlo& operator=(const MonteCarlo&) = delete;
	/** Use a move constructor instead */
	MonteCarlo(MonteCarlo&& old) :
		constants(std::move(old.constants)),
		c(std::move(old.c)),
		o(std::move(old.o)),
		p(std::move(old.p)),
		l(std::move(old.l)),
		paramSpecs(std::move(old.paramSpecs)),
		rng(std::move(old.rng)),
		resampler(std::move(old.resampler)){}

	void iterate(unsigned long);
	void iterate();

	oParams& getOByReference() {
		return o;
	}
	cParams& getCByReference() {
		return c;
	}

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const MonteCarlo&);
};