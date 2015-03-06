#include "Priors.hpp"

std::vector<ParamPriors<LaplacePrior>> Prior::makeCPriors(cParams &aC, LaplacePrior &aLaplacePrior) {
	std::vector<ParamPriors<LaplacePrior>> newCPriors;
	for (size_t i = 0; i < aC.size(); ++i) {
		newCPriors.push_back(ParamPriors<LaplacePrior>(aC.getPointerToC(i), aLaplacePrior));
	}
	return newCPriors;
}

std::vector<std::vector<ParamPriors<BetaPrior>>> Prior::makeOPriors(oParams &anO, BetaPrior &aBetaPrior) {
	std::vector<std::vector<ParamPriors<BetaPrior>>> newOPriors;
	for (size_t sample = 0; sample < anO.size(); ++sample) {
		std::vector<ParamPriors<BetaPrior>> newSampleOPriors;
		for (size_t site = 0; site < anO.size(sample); ++ site) {
			newSampleOPriors.push_back(ParamPriors<BetaPrior>(anO.getPointerToO(sample, site), aBetaPrior));
		}
		newOPriors.push_back(newSampleOPriors);
	}
	return newOPriors;
}

double Prior::updateAll() {
	priorTotal = temptativePriorTotal = 0;
	// Loop over the Cs
	for (size_t i = 0; i < cPriors.size(); ++i) {
		priorTotal += cPriors.at(i).update();
	}
	// Loop over the Os
	for (size_t sample = 0; sample < oPriors.size(); ++sample) {
		std::vector<ParamPriors<BetaPrior>>& oSamplePrior = oPriors.at(sample);
		for (size_t site = 0; site < oSamplePrior.size(); ++ site) {
			priorTotal += oSamplePrior.at(site).update();
		}
	}
	return priorTotal;
}