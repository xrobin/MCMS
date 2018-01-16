#' Number of parameters
#' @param object the \code{Peptides} model
#' @return an integer
#' @include PeptidesClass.R
setGeneric("nParams", function(object) {
  standardGeneric("nParams")
})


#' @describeIn nParams This corresponds to the sum of the concentration and occupancy parameters
setMethod("nParams", signature("Peptides"), function(object) {
  object@num.c + object@num.o
})

