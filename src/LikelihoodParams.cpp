#include "LikelihoodParams.hpp"

void LikelihoodParams::updateRedundantC() {
	const Rcpp::NumericMatrix& sampleDependence = constants.sampleDependence;
	for (int row = 0; row < sampleDependence.nrow(); ++row) {
		redundantC[row] = 0;
		for (int col = 0; col < sampleDependence.ncol(); ++col) {
			int elem = sampleDependence(row, col);
			if (elem != 0) {
				redundantC[row] += elem * c[col];
			}
		}
	}
}