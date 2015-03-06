#pragma once

#include "Likelihood.hpp"
#include "Priors.hpp"

class ParamSpecs {
	public:
	enum ParamCategory { c, o };
	const ParamCategory category;
	const size_t index1, const index2; // 1: c or sample (o); 2, o only: site

	ParamSpecs(const ParamCategory& aCategory, const size_t anIndex1): // c
		category(aCategory), index1(anIndex1) {}
	ParamSpecs(const ParamCategory& aCategory, const size_t anIndex1, const size_t anIndex2): // o
		category(aCategory), index1(anIndex1), index2(anIndex2) {}
};

/** A linear representation of all the parameters */
class ParamSpecsVector {
	std::mt19937_64& rng;
	std::uniform_int_distribution<size_t> dist;
	const std::vector<ParamSpecs> paramSpecs;

	public:
	ParamSpecsVector(const std::vector<ParamSpecs>& someParamSpecs, std::mt19937_64& anRNG):
		rng(anRNG), dist(0, someParamSpecs.size()), paramSpecs(someParamSpecs) {}

	ParamSpecs& getRandomParamSpecByReference() {
		return paramSpecs[dist(rng)];
	}
};

class MonteCarlo {
	cParams c;
	oParams o;
	Prior p;
	Likelihood l;
	const ParamSpecsVector paramSpecs; // a vector with all the parameters
	std::mt19937_64& rng; // The random numbers for the whole MC run

	static ParamSpecsVector makeParamSpecsVector(const cParams& c, const oParams& o);

	public:
	MonteCarlo(const Likelihood& aLikelihood,
		const cParams& aC, const oParams& anO,
		std::mt19937_64& anRNG):
		c(aC),
		o(anO),
		p(aC, anO),
		l(aLikelihood),
		paramSpecs(makeParamSpecs(aC, anO), anRNG),
		rng(anRNG){}

	void iterate(unsigned long);
	void iterate();
};