#include <cmath> // abs, pow
#include "VarianceModel.hpp"

double VarianceModel::calcPowerPriorRate(const double ratio) const {
	return rate0 + std::pow((std::abs(ratio) / rate1), rate2);
}

double VarianceModel::calcConstantPriorRate(const double ratio) const {
	return rate0;
}

double VarianceModel::calcPriorRate(const double ratio) const {
	if (rateType == RateType::constant) {
		return calcConstantPriorRate(ratio);
	}
	else {
		return calcPowerPriorRate(ratio);
	}
}