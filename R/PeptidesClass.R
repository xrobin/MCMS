check_peptides <- function(object) {

	errors <- character()
	if (length(object@c) + 1 != length(object@o)) {
		# TODO: make sure this is the general case and not a special case that we're having now!
		msg <- paste0("Expected to have n concentrations and n+1 occupancy numbers, not n + ", length(object@o) - length(object@c))
		errors <- c(errors, msg)
	}
	if (length(object@c) != object@num.c || length(object@c) != length(object@names.c) ) {
		msg <- paste0("Inconsistent number of parameters c + ", length(object@c), ", not ", object@num.c, " or ", length(object@names.c))
		errors <- c(errors, msg)
	}
	if (length(object@o) != object@num.o || length(object@o) != length(object@names.o) ) {
		msg <- paste0("Inconsistent number of parameters o + ", length(object@o), ", not ", object@num.o, " or ", length(object@names.o))
		errors <- c(errors, msg)
	}

	if (length(errors) == 0) TRUE else errors
}

#' @name Peptides-class
#' @rdname Peptides-class
#' @exportClass Peptides
#' @title A class that describe a set of peptides. Use the \code{\link{PeptidesModel}} function for easy object creation
#' @slot c the concentration ratios (per sample pair) as a named numeric (name is sampleX_sampleY, ...)
#' @slot o the occupancy ratios (per sample) as a named numeric (name is sampleX, ...)
#' @slot num.c,num.o number of o and c parameters
#' @slot names.c,names.o names of the pair to which each c, and of the sample to which each o applies.
#' @include ProteinClass.R
setClass("Peptides",
		 representation(
		 	c = "numeric",
		 	o = "numeric",
		 	num.c = "numeric",
		 	num.o = "numeric",
		 	names.c = "character",
		 	names.o = "character",
		 	protein = "Protein"),
		 validity = check_peptides)

#' Creates a model with concentration and occupancies for a Protein.
#' @param proteins a \code{\link{Protein-class}} object
#' @import dplyr
#' @examples
#' data(ENSTest)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- PeptidesModel(ENSTestProtein)
#' @export
PeptidesModel <- function(protein) {
	# Compute the number of occupancies and concentrations
	names.c <- colnames(protein@sample.dependency)
	num.c <- length(names.c)

	# Take the median ratios as starting point for the c
	median.ratios <- protein@data %>%
		filter(modifications == "") %>%
		group_by(pair) %>%
		summarize(ratio = median(ratio)) %>%
		slice(match(names.c, pair))
	c.initial <- median.ratios$ratio

	# Set the o's to 0.5 for now
	o.initial <- matrix(0.5,
						nrow = length(protein@reference.sample.intersect),
						ncol = length(protein@modifications),
						dimnames = list(
							sample = newdata@reference.sample.intersect,
							site = newdata@modifications
							)
						)
	names.o <- dimnames(o.initial)
	num.o <- prod(dim(o.initial))
	# o's with no covering a sample should be NA

	peptides <- new("Peptides",
		c = c.initial,
		o = o.initial,
		num.c = num.c,
		num.o = num.o,
		names.c = names.c,
		names.o = names.o,
		protein = protein)

	return(peptides)
}