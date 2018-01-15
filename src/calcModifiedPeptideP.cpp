#include <Rcpp.h>
using namespace Rcpp;

// This function pushes back an NA (if currentProb is the empty string "") or the currentString converted to double, into probs
// Modifies prob in place
void push_backNumberOrZero(std::vector<double>& probs, std::string& currentProb) {
	if (currentProb == "") {
		probs.push_back(0);
	}
	else {
		probs.push_back(std::stod(currentProb));
	}
}


// Given a sequence probability string of the form "AGDS(0.995)DS(0.998)WDADAFS(0.007)VEDPVRK",
// This function returns a vector of probabilites per position in the peptide sequence
std::vector<double> getModificationProbabilities(std::string& ProbabilitySequence) {
	const char *numberChars = "0.123456789";
	std::string currentProb = "";
	std::vector<double> probs;

	bool firstChar = true; // Do not push at first character

	if (ProbabilitySequence.size() > 0) {
		for(char& ps : ProbabilitySequence) {
			if (isupper(ps)) {
				if (! firstChar) {
					push_backNumberOrZero(probs, currentProb);
					currentProb = "";
				}
				else {
					firstChar = false;
				}
			}
			else if (strchr(numberChars, ps))  {
				currentProb.append(1, ps);
			}
		}
		// push_back the final number
		push_backNumberOrZero(probs, currentProb);
	}

	return probs;
}

double calcModifiedPeptideP(std::string ModifiedSequence, std::string PhosphoProbabilitySequence, std::string OxProbabilitySequence) {
	double p = 1;

	// Get Phospho Probabilites
	std::vector<double> phosphoProbabilities = getModificationProbabilities(PhosphoProbabilitySequence);
	std::vector<double> oxProbabilities = getModificationProbabilities(OxProbabilitySequence);

	// Walk the Modified sequence
	std::string modifications = "";
	int positionInPeptideSequence = 0;
	for(char& ms : ModifiedSequence) {
		if (isupper(ms)) {
			//Rcout << positionInPeptideSequence << ": " << ms;
			//if (modifications != "") {
			//	Rcout << " (" << modifications << ")";
			//}

			if (modifications == "(ph)") {
				if (positionInPeptideSequence > phosphoProbabilities.size()) {
					//Rcout << "skipping missing ph: " << positionInPeptideSequence << " > " << phosphoProbabilities.size() << std::endl;
					return 0;
				}
				//Rcout << " p[" << positionInPeptideSequence - 1 << "] = " << phosphoProbabilities[positionInPeptideSequence - 1];
				p = p * phosphoProbabilities[positionInPeptideSequence - 1];
			}
			else if (modifications == "(ox)") {
				if (positionInPeptideSequence > oxProbabilities.size()) {
					//Rcout << "skipping missing ox: " << positionInPeptideSequence << " > " << oxProbabilities.size() << std::endl;
					return 0;
				}
				//Rcout << " o[" << positionInPeptideSequence - 1 << "] = " << oxProbabilities[positionInPeptideSequence - 1];
				p = p * oxProbabilities[positionInPeptideSequence - 1];
			}
			//Rcout << std::endl;

			++positionInPeptideSequence;
			modifications = "";
		}
		else {
			modifications.append(1, ms);
		}
	}

	//Rcout << "Last: ";
	if (modifications == "(ph)_") {
		if (positionInPeptideSequence > phosphoProbabilities.size()) {
			//Rcout << "skipping final p: " << positionInPeptideSequence << " > " << phosphoProbabilities.size() << std::endl;
			return 0;
		}
		//Rcout << " p[" << positionInPeptideSequence - 1 << "] = " << phosphoProbabilities[positionInPeptideSequence - 1];
		p = p * phosphoProbabilities[positionInPeptideSequence - 1];
	}
	else if (modifications == "(ox)_") {
		if (positionInPeptideSequence > oxProbabilities.size()) {
			//Rcout << "skipping final ox: " << positionInPeptideSequence << " > " << oxProbabilities.size() << std::endl;
			return 0;
		}
		//Rcout << " o[" << positionInPeptideSequence - 1 << "] = " << oxProbabilities[positionInPeptideSequence - 1];
		p = p * oxProbabilities[positionInPeptideSequence - 1];
	}
	//Rcout << std::endl;
	return p;
}


// [[Rcpp::export]]
std::vector<double> calcModifiedPeptidePCpp(std::vector<std::string> ModifiedSequences, std::vector<std::string> PhosphoProbabilitySequences, std::vector<std::string> OxProbabilitySequences) {
	size_t s = ModifiedSequences.size();
	if (s != PhosphoProbabilitySequences.size() || s != OxProbabilitySequences.size()) {
		stop("Incompatible sizes");
	}

	std::vector<double> p(s);

	for (size_t i = 0; i < s; i++) {
		checkUserInterrupt();
		p[i] = calcModifiedPeptideP(ModifiedSequences[i], PhosphoProbabilitySequences[i], OxProbabilitySequences[i]);
	}
	return p;
}
