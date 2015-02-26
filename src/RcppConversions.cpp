#include <Rcpp.h>


using namespace Rcpp;

//namespace Rcpp {
//	// LikelihoodParams
//	template <> LikelihoodParams as(Rcpp::S4& someParams) {
//		NumericMatrix sampleDependency = someParams.slot("sample.dependency");
//	}
//
//	template <> Rcpp::List wrap(const LikelihoodParams<> &someParams) {
//		stop("Not implemented");
//	}
//
//}