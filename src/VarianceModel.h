#pragma once

#include <ostream>

class VarianceModel {
	public:
		enum RateType { constant, power };
	private:
		const RateType rateType;
		double rate0, rate1, rate2, shape0; // the parameters

	public:
	VarianceModel(const RateType aRateType, const double aRate0, const double aRate1, const double aRate2, const double aShape0):
		rateType(aRateType), rate0(aRate0), rate1(aRate1), rate2(aRate2),
		shape0(aShape0) {};


	private:
	double calcPowerPriorRate(const double ratio) const;
	double calcConstantPriorRate(const double ratio) const;

	public:
	double calcPriorRate(const double ratio) const;

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

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const VarianceModel&);
};