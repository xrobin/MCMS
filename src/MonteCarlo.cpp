#include <iostream>
#include "MonteCarlo.hpp"
#include "Resampler.hpp"

using std::cout;
using std::endl;

void MonteCarlo::iterate() {
	// Choose parameter to change
	const ParamSpecs& randomParam = paramSpecs.getRandomElementByReference(rng);
	cout << randomParam << endl;
	if (randomParam.category == ParamSpecs::c) {
		// 1. Get previous C
		double oldC = c.getC(randomParam.index1);
		// 2. resample a new c
		MoveSpec move = resampler.resampleC(randomParam, oldC);
		cout << "Move: " << move << endl;
		double newC = move.newParam;
		// 3. Update the parameter
		c.setC(randomParam.index1, newC);
		// 4. Query new likelihood
		double likelihoodChange = l.temptativeChangedC(randomParam.index1);
		cout << "Likelihood change: " << likelihoodChange << endl;
		// 5. Query new prior
		double priorChange = p.temptativeChangedC(randomParam.index1);
		cout << "Prior change: " << priorChange << endl;
		// 6 decide acceptance
		if (move.accept(likelihoodChange + priorChange, rng)) {
			cout << "Accepted!" << endl;
			l.acceptC(randomParam.index1);
			p.acceptC(randomParam.index1);
		}
		else {
			cout << "Rejected!" << endl;
			c.setC(randomParam.index1, oldC);
			// There is no need for explicit reject
		}

	}
	else { // we have an o
		double oldO = o.getO(randomParam.index1, randomParam.index2);
		MoveSpec move = resampler.resampleO(randomParam, oldO);
		cout << "Move: " << move << endl;
		double newO = move.newParam;
		o.setO(randomParam.index1, randomParam.index2, newO);
		double likelihoodChange = l.temptativeChangedO(randomParam.index1, randomParam.index2);
		cout << "Likelihood change: " << likelihoodChange << endl;
		double priorChange = p.temptativeChangedO(randomParam.index1, randomParam.index2);
		cout << "Prior change: " << priorChange << endl;
		if (move.accept(likelihoodChange + priorChange, rng)) {
			l.acceptO(randomParam.index1, randomParam.index2);
			p.acceptO(randomParam.index1, randomParam.index2);
		}
		else {
			o.setO(randomParam.index1, randomParam.index2, oldO);
		}
	}
	//throw std::runtime_error("Save the parameters & likelihood!");
}


void MonteCarlo::iterate(unsigned long i) {
	Rcpp::Rcout << i << "\n";
	while (i--) {
		iterate();
	}
}

ParamSpecsVector MonteCarlo::makeParamSpecsVector(cParams& c, oParams& o, Prior& p) {
	std::vector<ParamSpecs> paramSpecs;
	LaplacePrior& laplacePrior = p.getLaplacePrior();
	BetaPrior& betaPrior = p.getBetaPrior();

	// Fill the paramSpecs vector with o and c
	for (size_t i = 0; i < c.size(); ++i) {
		paramSpecs.push_back(ParamSpecs(ParamSpecs::ParamCategory::c, c.getPointerToC(i), i, laplacePrior));
	}
	for (size_t i = 0; i < o.size(); ++i) {
		for (size_t j = 0; j < o.size(i); ++j) {
			paramSpecs.push_back(ParamSpecs(ParamSpecs::ParamCategory::o, o.getPointerToO(i, j), i, j, betaPrior));
		}
	}
	return paramSpecs;

}