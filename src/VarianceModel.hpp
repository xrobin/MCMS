#pragma once

#include <cmath> // abs, pow

class VarianceModel {
	double rate0, rate1, rate2, shape0; // the parameters

	public:
	VarianceModel(const double aRate0, const double aRate1, const double aRate2, const double aShape0):
		rate0(aRate0), rate1(aRate1), rate2(aRate2), shape0(aShape0) {};

	double calcPriorRate(const double ratio) const {
		return rate0 + std::pow((std::abs(ratio) / rate1), rate2);
		//return rate0 + exp(rate2*log(fabs(ratio) / rate1));
		//return 2.1*3.1;
	}

	double calcPriorShape() const {
		return shape0;
	}

	double calcGamma(const double ratio, const double n, const double q) const {
		double b_s = 0.5 * q + calcPriorRate(ratio);
		return n / (2 * b_s);
	}

	double calcNuTilde(const double n) const {
		double nu = n - 1;
		double nu_s = nu + 2 * calcPriorShape();
		return (nu_s + 1) / 2;
	}
};