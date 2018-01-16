#' Calculate the number of iterations for the MCMS run
#' @description The calculation is based on the number of parameters (c and o)
#' @param Peptides the \code{Peptides}
#' @param round whether to round the number up to the next power of 10 (default). If \code{FALSE} the value is only rounded up to the next integer.
#' @return an integer
#' @examples 
#' data(ENSTest, var.model)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- Peptides(ENSTestProtein)
#' calcIterations(ENSTestModel)
#' @export
calcIterations <- function(Peptides, round = TRUE) {
  n.params <- Peptides@num.c + Peptides@num.o
  
  computedTotalNumberIterations <-
    # Computed in MQ_1.3_patient_data_MC.R
    100 / exp(9.227 + -1.898 * log(n.params)) * 1E7
  
  if (round) {
    return(10^ceiling(log10(computedTotalNumberIterations)))
  }
  else {
    return(ceiling(computedTotalNumberIterations))
  }
}
