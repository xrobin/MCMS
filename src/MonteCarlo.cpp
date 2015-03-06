#include "MonteCarlo.hpp"

void MonteCarlo::iterate() {
	// Choose parameter to change
	ParamSpecs& randomParam = paramsSpecs.getRandomElementByReference();
	if (randomParam.category == randomParam::ParamCategory::c) {
		// 1. resample a new c
		double newC = resample(c);
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

	}
}


void MonteCarlo::iterate(unsigned int i) {
	while (i--) {
		iterate();
	}
}

static ParamSpecsVector makeParamSpecsVector(const cParams& c, const oParams& o) {
	std::vector<ParamSpecs> paramSpecs;
	// Fill the paramSpecs vector with o and c
	for (size_t i = 0; i < c.size(); ++i) {
		paramSpecs.push_back(ParamSpecs(ParamSpecs::ParamCategory::c, i));
	}
	for (size_t i = 0; i < o.size(); ++i) {
		for (size_t j = 0; j < o.size(i); ++j) {
			paramSpecs.push_back(ParamSpecs(ParamSpecs::ParamCategory::o, i, j));
		}
	}
	return paramSpecs;

}