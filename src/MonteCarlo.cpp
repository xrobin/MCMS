#include "MonteCarlo.hpp"
#include "Resample.hpp"

void MonteCarlo::iterate() {
	// Choose parameter to change
	const ParamSpecs& randomParam = paramSpecs.getRandomElementByReference(rng);
	if (randomParam.category == ParamSpecs::ParamCategory::c) {
		// 1. Get previous C
		double oldC = c.getC(randomParam.index1);
		// 2. resample a new c
		MoveSpec move = resampler.resampleC(randomParam, oldC);
		double newC =
		// 3. Update the parameter
		c.setC(randomParam.index1, newC);
		// 4. Query new likelihood
		double likelihoodChange = l.temptativeChangedC(randomParam.index1);
		// 5. Query new prior
		double priorChange = p.temptativeChangedC(randomParam.index1);
		// 6 decide acceptance
		if (move.accept(likelihoodChange + priorChange, rng)) {
			l.accept(randomParam.index1);
			p.accept(randomParam.index1);
		}
		else {
			c.setC(randomParam.index1, oldC);
			// There is no need for explicit reject
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