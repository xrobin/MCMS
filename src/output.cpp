#include <iostream> // cout
#include "Likelihood.hpp"
#include "MonteCarlo.hpp"
#include "Parameters.hpp"
#include "Peptide.hpp"
#include "prettyprint.hpp"
#include "Priors.hpp"
#include "Resampler.hpp"
#include <vector>

using std::endl;
using std::ostream;
using std::string;
using std::vector;


ostream& operator<< (ostream &out, const Likelihood &Likelihood) {
	out << Likelihood.peptideLikelihoods << endl;
	out << "onpudate_c = " << Likelihood.onupdate_c << endl;
	out << "onpudate_o = " << Likelihood.onupdate_o << endl;
	return out;
}

ostream& operator<< (ostream &out, const PeptideLikelihood &aPeptideLikelihood) {
	out << "<PeptideLikelihood@" << &aPeptideLikelihood << ": {" << aPeptideLikelihood.likelihoodValue << ", ";
	out << aPeptideLikelihood.peptide;
	out << "}";
	return out;
}

ostream& operator<< (ostream &out, const Peptide &aPeptide) {
	out << "<Peptide@" << &aPeptide << ": " << aPeptide.pairName << " " << aPeptide.nEff << " x " << aPeptide.ratio << " ± " << aPeptide.q << "; ";
	out << "c = " << *(aPeptide.c) << " (" << aPeptide.c << ")>" << endl;
	for (const auto& SiteSpecs: aPeptide.siteSpecs) {
		out << SiteSpecs;
	}
	return out;
}

ostream& operator<< (ostream &out, const SiteSpecs &aSiteSpec) {
	out << "<SiteSpecs@" << &aSiteSpec << ": " << aSiteSpec.siteName << " (" << (aSiteSpec.siteActivity ? "on" : "off") << "); ";
	out << "o = {" << *(aSiteSpec.params.sample) << " (" << aSiteSpec.params.sample << "), ";
	out << *(aSiteSpec.params.reference) << " (" << aSiteSpec.params.reference << ")}>" << endl;
	return out;
}

ostream& operator<< (ostream &out, const GenericPrior &aPrior) {
	aPrior.print(out);
	return out;
}

void BetaPrior::print(ostream &out) const {
	out << "<BetaPrior@" << this << "(" << beta.alpha() << ", "  << beta.beta() << ")>";
}

//void LaplacePrior::print(ostream &out) const {
//	out << "<LaplacePrior@" << this << "(" << laplace.location() << ", "  << laplace.scale() << ")>";
//}

void NormalPrior::print(ostream &out) const {
	out << "<NormalPrior@" << this << "(" << normal.mean() << ", "  << normal.standard_deviation() << ")>";
}

ostream& operator<< (ostream& out, const VarianceModel& aVarianceModel) {
	out << "<VarianceModel@" << &aVarianceModel << " (";
	out << "rate0 = " << aVarianceModel.rate0;
	out << ", rate1 = " << aVarianceModel.rate1;
	out << ", rate2 = " << aVarianceModel.rate2;
	out << ", shape0 = " << aVarianceModel.shape0;
	out << ")>";
	return out;
}


ostream& operator<< (ostream& out, const Constants& someConstants) {
	out << "<Constants@" << &someConstants << " (";
	out << "varianceModel = " << someConstants.varianceModel << endl;
	out << "sampleDependenceMatrix = " << someConstants.sampleDependenceMatrix << endl;
	out << ", o_prior_shape1 = " << someConstants.o_prior_shape1;
	out << ", o_prior_shape2 = " << someConstants.o_prior_shape2;
	out << ", c_prior_sd = " << someConstants.c_prior_sd;
	out << ", o_restrict = " << someConstants.o_restrict;
	out << ", priorMoveProportion = " << someConstants.priorMoveProportion;
	out << ", c_sd = " << someConstants.c_sd;
	out << ", o_sd = " << someConstants.o_sd;
	out << ", o_k_scale = " << someConstants.o_k_scale;
	out << ")>";
	return out;
}

ostream& operator<< (ostream &out, const ParamSpecs &aParamSpecs) {
	out << "<ParamSpecs@" << &aParamSpecs << ": ";
	// Print the category as text (default for an enum is integer)
	switch(aParamSpecs.category) {
		case ParamSpecs::c: out << "c"; break;
		case ParamSpecs::o: out << "o"; break;
	}
	out << ", " << aParamSpecs.index1 << ", " << aParamSpecs.index2;
	out << ", sd_corr = " << aParamSpecs.sd_corr;
	out << "; param@" << aParamSpecs.param << " (" << *(aParamSpecs.param) << "), ";
	out << aParamSpecs.prior << ">";
	return out;
}

ostream& operator<< (ostream &out, const MoveSpec &aMoveSpec) {
	out << "<MoveSpec@" << &aMoveSpec << ": (";
	out << aMoveSpec.oldParam << " -> " << aMoveSpec.newParam;
	out << "; bias = " << aMoveSpec.logPreviousBias << " -> " << aMoveSpec.logNewBias;
	out << ")>";
	return out;
}

void printAC(ostream &out, const vector<double> &aC, const std::unordered_map<std::string, size_t> &aCNames) {
	for (size_t i = 0; i < aC.size(); ++i) {
		// get the cName
		auto it = aCNames.begin();
		while (it->second != i) {
			++it;
		};
		out << i << " (" << it->first << "): " << aC[i] << " (" << &(aC[i]) << ")";
		if (i < aC.size() - 1) {
			out << ", ";
		}
	}
}

ostream& operator<< (ostream &out, const cParams &aCParams) {
	out << "c = [";
	printAC(out, aCParams.c, aCParams.cNames);
	out << "]" << endl;

	out << "redundantC = [";
	printAC(out, aCParams.redundantC, aCParams.redundantCNames);
	out << "]" << endl;

	out << "redundantCToC = " << aCParams.redundantCToC << endl;

//	const Rcpp::NumericMatrix sampleDependence;
//	std::vector<double> c, redundantC;
//	std::unordered_map<std::string, size_t> cNames, redundantCNames; // so from a name we know which c to update
//	std::vector<std::vector<size_t>> redundantCToC; // Maps a redundantC to one or more c on

	return out;
}

ostream& operator<< (ostream &out, const oParams &anOParams) {
	out << "o = [";
	auto &o = anOParams.o;
	for (size_t i = 0; i < o.size(); ++i) {
		auto &oSample = o[i];
		// Get sample name
		auto it = anOParams.sampleNames.begin();
		while (it->second != i) {
			++it;
		};
		string sampleName = it->first;
		out << i << " (" << sampleName << "): {";
		for (size_t j = 0; j < oSample.size(); ++j) {
			// Get site name
			auto it = anOParams.siteNames[i].begin();
			while (it->second != j) {
				++it;
			};
			string siteName = it->first;
			out << j << " (" << siteName << "): ";
			out << oSample[j] << " (" << &(oSample[j]) << ")";
			if (j < oSample.size() - 1) {
				out << ", ";
			}
		}
		out << "}";
		if (i < o.size() - 1) {
			out << ", ";
		}
	}
	out << "]" << endl;
	return out;
}


ostream& operator<< (ostream& out, const Prior& aPrior) {
	out << "<Prior@" << &aPrior << ":" << endl;
	//out << aPrior.laplace << endl;
	out << aPrior.normal << endl;
	out << aPrior.beta << endl;
	out << "cPriors = " << aPrior.cPriors << endl;
	out << "oPriors = " << aPrior.oPriors;
	out << ">" << endl;
	return out;
}


ostream& operator<< (ostream &out, const MonteCarlo &aMonteCarlo) {
	out << "<MonteCarlo@" << &aMonteCarlo << ":" << endl;
	out << aMonteCarlo.constants << endl;
	out << aMonteCarlo.c << endl;
	out << aMonteCarlo.o << endl;
	//out << aMonteCarlo.p << endl;
	out << aMonteCarlo.l << endl;
	//out << aMonteCarlo.paramSpecs << endl;
	out << "<Resampler@" << &(aMonteCarlo.resampler)  << ">" << endl;
	out << ">" << endl;
	return out;
}