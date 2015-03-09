#include "MonteCarlo.hpp"
#include "Resample.hpp"

void MonteCarlo::iterate() {
	// Choose parameter to change
	const ParamSpecs& randomParam = paramSpecs.getRandomElementByReference(rng);
	if (randomParam.category == ParamSpecs::ParamCategory::c) {
		// 1. resample a new c
		double newC = resampleC(randomParam, constants);
		// 2. Update the parameter
		double oldC = c.setC(randomParam.index1, newC);
		// 3. Query new likelihood
		// 4. Query new prior
		// 5 decide acceptance
		if (accept) {
			l.accept(randomParam.index1);
			prior.accept(randomParam.index1);
		}
		else {
			c.setC(randomParam.index1, oldC);
		}

	}
	else { // we have an o
		double newO = resampleO(randomParam, constants);
		if (newO is out of range) {
			...
		}

	}
}


void MonteCarlo::iterate(unsigned long i) {
	while (i--) {
		iterate();
	}
}

ParamSpecsVector MonteCarlo::makeParamSpecsVector(const cParams& c, const oParams& o, Prior& p) {
	std::vector<ParamSpecs> paramSpecs;
	// Fill the paramSpecs vector with o and c
	for (size_t i = 0; i < c.size(); ++i) {
		paramSpecs.push_back(ParamSpecs(ParamSpecs::ParamCategory::c, i, p));
	}
	for (size_t i = 0; i < o.size(); ++i) {
		for (size_t j = 0; j < o.size(i); ++j) {
			paramSpecs.push_back(ParamSpecs(ParamSpecs::ParamCategory::o, i, j, p));
		}
	}
	return paramSpecs;

}