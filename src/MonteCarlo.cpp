#include <iostream>
#include <iterator> // std::back_inserter
#include "MonteCarlo.hpp"
#include "Resampler.hpp"
#include <string>
#include <vector>

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using std::cout;
using std::endl;
using std::string;
using std::vector;

NumericVector MonteCarlo::iterate() {
	// Choose parameter to change
	const ParamSpecs& randomParam = paramSpecs.getRandomElementByReference(rng);
//	cout << endl << randomParam << endl;
	if (randomParam.category == ParamSpecs::c) {
		// 1. Get previous C
		double oldC = c.getC(randomParam.index1);
		// 2. resample a new c
		MoveSpec move = resampler.resampleC(randomParam, oldC);
//		cout << "Move: " << move << endl;
		if (move.valid) {
			// 3. Update the parameter
			c.setC(randomParam.index1, move.newParam);
			// 4. Query new likelihood
			double likelihoodChange = l.temptativeChangedC(randomParam.index1);
//			cout << "Likelihood change: " << likelihoodChange << endl;
			// 5. Query new prior
			double priorChange = p.temptativeChangedC(randomParam.index1);
//			cout << "Prior change: " << priorChange << endl;
			// 6 decide acceptance
			if (move.accept(likelihoodChange + priorChange, rng)) {
//				cout << "Accepted!" << endl;
				l.acceptC(randomParam.index1);
				p.acceptC(randomParam.index1);
			}
			else {
//				cout << "Rejected!" << endl;
				c.setC(randomParam.index1, oldC);
				// There is no need for explicit reject
			}
		}
	}
	else { // we have an o
		double oldO = o.getO(randomParam.index1, randomParam.index2);
		MoveSpec move = resampler.resampleO(randomParam, oldO);
//		cout << "Move: " << move << endl;
		if (move.valid) {
			o.setO(randomParam.index1, randomParam.index2, move.newParam);
			double likelihoodChange = l.temptativeChangedO(randomParam.index1, randomParam.index2);
//			cout << "Likelihood change: " << likelihoodChange << endl;
			double priorChange = p.temptativeChangedO(randomParam.index1, randomParam.index2);
//			cout << "Prior change: " << priorChange << endl;
			if (move.accept(likelihoodChange + priorChange, rng)) {
//				cout << "Accepted!" << endl;
				l.acceptO(randomParam.index1, randomParam.index2);
				p.acceptO(randomParam.index1, randomParam.index2);
			}
			else {
//				cout << "Rejected!" << endl;
				o.setO(randomParam.index1, randomParam.index2, oldO);
			}
		}
	}

	// Save the parameters
	NumericVector McResult(paramSpecs.size() + 2);
	int i = 0;
	for (; i < paramSpecs.size(); ++i) {
		McResult[i] = *(paramSpecs[i].param);
	}
	McResult[i++] = l.getLikelihoodValue();
	McResult[i++] = p.getPriorTotal();
	//McResult.reserve(ParamSpecs.size() + 2)
	//for (const ParamSpecs& val: paramSpecs) {
	//	McResult.push_back(*(val.param));
	//}
	//transform(paramSpecs.begin(), paramSpecs.end(), std::back_inserter(McResult),
    //    [](ParamSpecs const& val) {return *(val.param);});
	// Also likelihood and prior
	//McResult.push_back(l.getLikelihoodValue());
	//McResult.push_back(p.getPriorTotal());
	return McResult;
}


NumericMatrix MonteCarlo::iterate(unsigned long i) {
	NumericMatrix McResult(i, paramSpecs.size() + 2);

	for (unsigned long j = 0; j < i; ++j) {
		McResult(j, Rcpp::_) = iterate();
	}
	return McResult;
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



vector<string> MonteCarlo::getIterateNames() {
	vector<string> names = getParamNames();
	names.push_back("Likelihood");
	names.push_back("Prior");
	return names;
}

vector<string> MonteCarlo::getParamNames() {
	vector<string> names;
	names.reserve(paramSpecs.size() + 2);
	for (const ParamSpecs& par: paramSpecs) {
		if (par.category == ParamSpecs::c) {
			names.push_back(std::string("c.") + c.getName(par.index1));
		}
		else {
			string sample = o.getSampleName(par.index1);
			string site = o.getSiteName(par.index1, par.index2);
			names.push_back(std::string("o.") + sample + "." + site);
		}
	}
	return names;
}