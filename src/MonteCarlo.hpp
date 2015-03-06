#pragma once

#include <boost/random/uniform_int_distribution.hpp>
#include "Helpers.hpp"
#include "Likelihood.hpp"
#include "Priors.hpp"

class ParamSpecs {
	public:
	enum ParamCategory { c, o };
	const ParamCategory category;
	const size_t index1, index2; // 1: c or sample (o); 2, o only: site

	ParamSpecs(const ParamCategory& aCategory, const size_t anIndex1): // c
		category(aCategory), index1(anIndex1), index2(0) {}
	ParamSpecs(const ParamCategory& aCategory, const size_t anIndex1, const size_t anIndex2): // o
		category(aCategory), index1(anIndex1), index2(anIndex2) {}
};

/** A linear representation of all the parameters */
typedef RandomizingConstVector<ParamSpecs> ParamSpecsVector;

class MonteCarlo {
	cParams c;
	oParams o;
	Prior p;
	Likelihood l;
	const ParamSpecsVector paramSpecs; // a vector with all the parameters
	std::mt19937_64& rng; // The random numbers for the whole MC run

	static ParamSpecsVector makeParamSpecsVector(const cParams& c, const oParams& o);

	public:
	MonteCarlo(const std::vector<Peptide>& peptides,
		const cParams::c_type& aCMap, const oParams::o_type& anOMap,
		const LikelihoodConstants& constants,
		std::mt19937_64& anRNG):
		c(aCMap, constants.sampleDependenceMatrix),
		o(anOMap),
		p(c, o, constants.scale, constants.shape1, constants.shape2),
		l(peptides, c, o, constants),
		paramSpecs(makeParamSpecsVector(c, o)),
		rng(anRNG){}

	void iterate(unsigned long);
	void iterate();
};