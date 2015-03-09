#pragma once

#include <boost/random/uniform_int_distribution.hpp>
#include "Helpers.hpp"
#include "Likelihood.hpp"
#include "Priors.hpp"
#include <stdexcept>

class ParamSpecs {
	public:
	enum ParamCategory { c, o };
	const ParamCategory category;
	const size_t index1, index2; // 1: c or sample (o); 2, o only: site
	GenericPrior& prior;

	ParamSpecs(const ParamCategory& aCategory, const size_t anIndex1, GenericPrior& aPrior): // c
		category(aCategory), index1(anIndex1), index2(0), prior(aPrior) {
			if (category != c) {
				throw std::invalid_argument("One argument constructor expects to receive a 'c' category param");
			}
		}
	ParamSpecs(const ParamCategory& aCategory, const size_t anIndex1, const size_t anIndex2, GenericPrior& aPrior): // o
		category(aCategory), index1(anIndex1), index2(anIndex2), prior(aPrior) {
			if (category != o) {
				throw std::invalid_argument("Two argument constructor expects to receive an 'o' category param");
			}
		}
};

/** A linear representation of all the parameters */
//typedef RandomizingConstVector<ParamSpecs> ParamSpecsVector; // can't be const because of the distributions :(
typedef RandomizingVector<ParamSpecs> ParamSpecsVector;

class MonteCarlo {
	const Constants constants;
	cParams c;
	oParams o;
	Prior p;
	Likelihood l;
	const ParamSpecsVector paramSpecs; // a vector with all the parameters
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
		rng(anRNG){}

	void iterate(unsigned long);
	void iterate();
};