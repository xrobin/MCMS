#' Calculate modified peptide probability
#' @description Integrates the Modification Probabilities from MaxQuant into a localization probability of the peptide
#' @param ModifiedSequences the "Modified Sequence" column in evidence.txt
#' @param PhosphoProbabilitySequences the "Phospho (STY) Probabilities" column in evidence.txt
#' @param OxProbabilitySequences the "Oxidation (M) Probabilities" column in evidence.txt
#' @return the numeric probabilities
#' @useDynLib MCMS
#' @importFrom Rcpp evalCpp
#' @examples
#' calcModifiedPeptideP("_(ac)AGDS(ph)DSWDADAFSVEDPVRK_", "AGDS(1)DSWDADAFSVEDPVRK", "")
#' calcModifiedPeptideP("_AAFNSGKVDIVAINDPFIDLNYM(ox)VYM(ox)FQYDSTHGK_", "",
#'                      "AAFNSGKVDIVAINDPFIDLNYM(1)VYM(1)FQYDSTHGK")
#' calcModifiedPeptideP("_AAEM(ox)CY(ph)RK_", "AAEMCY(1)RK", "AAEM(1)CYRK")
#'
#' # Also vectorized:
#' calcModifiedPeptideP(
	#'                      c("_AAEM(ox)CY(ph)RK_", "_(ac)AGDS(ph)DSWDADAFSVEDPVRS(ph)_",
#'                        "_(ac)AGDS(ph)DSWDADAFSVEDPVRM(ox)_"),
#'                      c("AAEMCY(0.99)RK", "AGDS(0.995)DSWDADAFS(0.007)VEDPVRS(0.998)",
#'                        "AGDS(0.995)DSWDADAFS(0.007)VEDPVRM"),
#'                      c("AAEM(0.87)CYRK", "", "AGDSDSWDADAFSVEDPVRM(0.998)"))
#' @export
calcModifiedPeptideP <- function(ModifiedSequences, PhosphoProbabilitySequences = NULL, OxProbabilitySequences = NULL) {
	# Make sure we have both Phospho and Ox ProbabilitySequences
	if (missing(PhosphoProbabilitySequences) || is.null(PhosphoProbabilitySequences)) {
		PhosphoProbabilitySequences <- rep("", length(ModifiedSequences))
	}
	if (missing(OxProbabilitySequences) || is.null(OxProbabilitySequences)) {
		OxProbabilitySequences <- rep("", length(ModifiedSequences))
	}
	return(calcModifiedPeptidePCpp(ModifiedSequences, PhosphoProbabilitySequences, OxProbabilitySequences))
}
