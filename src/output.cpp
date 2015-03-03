#include <iostream> // cout
#include "Likelihood.hpp"
#include "LikelihoodParams.hpp"
#include "Peptide.hpp"
#include "prettyprint.hpp"
#include <vector>

using std::endl;
using std::ostream;
using std::string;
using std::vector;


ostream& operator<< (ostream &out, const Likelihood &Likelihood) {
	out << Likelihood.peptideLikelihoods << endl;
	out << Likelihood.c;
	out << Likelihood.o;
	out << "onpudate_c = " << Likelihood.onupdate_c << endl;
	out << "onpudate_o = " << Likelihood.onupdate_o << endl;
	return out;
}

ostream& operator<< (ostream &out, const PeptideLikelihood &aPeptideLikelihood) {
	out << "<PeptideLikelihood@" << &aPeptideLikelihood << ": {" << aPeptideLikelihood.likelihood << ", ";
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