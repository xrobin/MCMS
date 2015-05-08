#pragma once

#include <boost/random/uniform_int_distribution.hpp>
#include "Likelihood.hpp"
#include "Priors.hpp"
#include <Rcpp.h>
#include "Resampler.hpp"
#include <string>
#include <vector>

class MonteCarlo {
	const Constants constants;
	cParams c;
	oParams o;
	Prior p;
	Likelihood l;
	const ParamSpecsVector paramSpecs; // a vector with all the parameters
	std::mt19937_64 rng; // The random numbers for the whole MC run
	Resampler resampler;

	static ParamSpecsVector makeParamSpecsVector(cParams& c, oParams& o, Prior& p, Likelihood& l);

	public:
	MonteCarlo(const std::vector<Peptide>& peptides,
		const cParams::c_type& aCMap, const oParams::o_type& anOMap,
		const Constants& someConstants,
		const std::mt19937_64& anRNG):
		constants(someConstants),
		c(aCMap, constants.sampleDependenceMatrix),
		o(anOMap),
		p(c, o, constants.c_prior_sd, constants.o_prior_shape1, constants.o_prior_shape2, constants.o_restrict),
		l(peptides, c, o, constants),
		paramSpecs(makeParamSpecsVector(c, o, p, l)),
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

	Rcpp::NumericMatrix iterate(const unsigned long, const unsigned long, const unsigned long);
	void iterate(const double cooling_rate = 1);
	Rcpp::NumericVector recordState();

	std::vector<std::string> getIterateNames();
	std::vector<std::string> getParamNames();

	const oParams& getOByReference() const {
		return o;
	}
	const cParams& getCByReference() const {
		return c;
	}

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const MonteCarlo&);
};