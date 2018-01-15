// [[Rcpp::depends(BH)]]
#include <boost/algorithm/string/join.hpp>
#include <Rcpp.h>
using namespace Rcpp;

// This function takes a Modified sequence and optional origin (typically 1-based) return the modifications (T_12, a, M_72) as a vector<string>.
std::vector<std::string> constructModifiedPeptide(std::string ModifiedSequence, int origin = 1) {

	std::vector<std::string> modifications;
	std::string buffer = "";
	std::string previousAA = "";
	int positionInPeptideSequence = origin;

	// Walk the Modified sequence
	for(char& ms : ModifiedSequence) {
		if (isupper(ms)) {

			if (buffer == "(ph)" || buffer == "(ox)") {
				modifications.push_back(std::string(previousAA) + "_" + std::to_string(positionInPeptideSequence - 1));
			}
			else if (buffer == "_(ac)") {
				modifications.push_back("a");
			}

			++positionInPeptideSequence;
			buffer = "";
			previousAA = ms;
		}
		else {
			buffer.append(1, ms);
		}
	}

	if (buffer == "(ph)_" || buffer == "(ox)_") {
		modifications.push_back(std::string(previousAA) + "_" + std::to_string(positionInPeptideSequence - 1));
	}


	return modifications;
}

// This function acts like constructModifiedPeptide but returns the modifications as a string, delimited by 'delimiter'.
std::string constructModifiedPeptideAsString(std::string ModifiedSequence, int origin = 1, std::string delimiter = ";") {
	return boost::algorithm::join(constructModifiedPeptide(ModifiedSequence, origin), delimiter);
}


//' Construct the modified peptide sequence as a string
//' @description Takes the Modified Sequence from MaxQuant and turns it into a string of the form a;T_12;M_72.
//' @param ModifiedSequences the "Modified Sequence" column in evidence.txt
//' @param positions the peptides start positions within the protein. Recycled according to the standard rules if necessary, with a warning.
//' @param delimiter an optional delimiter to separate the modifications within one peptide
//' @return the modifications string
//' @importFrom Rcpp evalCpp
//' @examples
//' constructModifiedPeptide("_(ac)AGDS(ph)DSWDADAFSVEDPVRK_", 1)
//' constructModifiedPeptide("_AAFNSGKVDIVAINDPFIDLNYM(ox)VYM(ox)FQYDSTHGK_", 20)
//' constructModifiedPeptide("_AAEM(ox)CY(ph)RK_", 10, ":")
//'
//' # Also vectorized:
//' constructModifiedPeptide(c("_(ac)AGDS(ph)DSWDADAFSVEDPVRK_",
//'                        "_AAFNSGKVDIVAINDPFIDLNYM(ox)VYM(ox)FQYDSTHGK_",
//'                        "_AAEM(ox)CY(ph)RK_"), c(1, 20, 10), "+")
//' @export
// [[Rcpp::export]]
CharacterVector constructModifiedPeptide(CharacterVector ModifiedSequences, IntegerVector positions = IntegerVector(), std::string delimiter = ";") {
	R_xlen_t s = ModifiedSequences.size();
	if (positions.size() == 0) {
		positions.push_back(1);
	}
	if (positions.size() < s) {
		if (positions.size() > 1) {
			warning("Recycling shorter 'positions' vector");
		}
		positions = rep_len(positions, s);
	}
	else if (positions.size() > s) {
		stop("Incompatible sizes");
	}

	CharacterVector p(s);

	for (R_xlen_t i = 0; i < s; i++) {
		checkUserInterrupt();
		p[i] = constructModifiedPeptideAsString(std::string(ModifiedSequences[i]), positions[i], delimiter);
	}
	return p;
}
