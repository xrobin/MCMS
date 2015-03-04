#include <boost/math/distributions/beta.hpp>
#include <cmath>
//#include <boost/math/distributions/normal.hpp>
//#include <random>


class BetaPrior {
	boost::math::beta_distribution<double> bd;

	public:
	BetaPrior(const double aShape1, const double aShape2):
		bd(aShape1, aShape2) {}
	double operator(const double x) {
		return bd(x);
	}
}

class NormalPrior {
	public:
	NormalPrior(): {}
	double operator(const double x) {
		return - std::abs(x) * 2;
	}
}

