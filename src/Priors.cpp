#include "Priors.h"

std::vector<ParamPriors<NormalPrior>> Prior::makeCPriors(cParams &aC, NormalPrior &aNormalPrior) {
	std::vector<ParamPriors<NormalPrior>> newCPriors;
	for (size_t i = 0; i < aC.size(); ++i) {
		newCPriors.push_back(ParamPriors<NormalPrior>(aC.getPointerToC(i), aNormalPrior));
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

double Prior::changedO(const size_t sample, const size_t site) {
	double previousPriorTotal = priorTotal;
	priorTotal += oPriors.at(sample).at(site).changed();
	return priorTotal - previousPriorTotal;
}
double Prior::temptativeChangedO(const size_t sample, const size_t site) {
	temptativePriorTotal = priorTotal + oPriors.at(sample).at(site).temptativeChanged();
	return temptativePriorTotal - priorTotal;
}
double Prior::acceptO(const size_t sample, const size_t site) {
	oPriors.at(sample).at(site).accept();
	priorTotal = temptativePriorTotal;
	return priorTotal;
}
double Prior::changedC(const size_t i) {
	double previousPriorTotal = priorTotal;
	priorTotal += cPriors.at(i).changed();
	return priorTotal - previousPriorTotal;
}
double Prior::temptativeChangedC(const size_t i) {
	temptativePriorTotal = priorTotal + cPriors.at(i).temptativeChanged();
	return temptativePriorTotal - priorTotal;
}
double Prior::acceptC(const size_t i) {
	cPriors.at(i).accept();
	priorTotal = temptativePriorTotal;
	return priorTotal;
}