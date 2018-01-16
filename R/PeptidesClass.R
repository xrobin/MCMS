check_peptides <- function(object) {
	errors <- character()
	if (length(object@c) != object@num.c || length(object@c) != length(object@names.c) ) {
		msg <- paste0("Inconsistent number of parameters c + ", length(object@c), ", not ", object@num.c, " or ", length(object@names.c))
		errors <- c(errors, msg)
	}
	# Ensure o and its name have the same dimensions
	o.dim <- sapply(object@o, length, simplify = FALSE)
	o.names.dim <- sapply(object@o, length, simplify = FALSE)
	if (!identical(o.dim, o.names.dim)) {
		msg <- paste0("Inconsistent number of parameters o.")
	}

	if (length(errors) == 0) TRUE else errors
}

#' @name Peptides-class
#' @rdname Peptides-class
#' @exportClass Peptides
#' @aliases Peptides
#' @title Set of peptides
#' @description A class that describe a set of peptides. Use the \code{\link{PeptidesModel}} function for easy object creation
#' @slot c the concentration ratios (per sample pair) as a named numeric (name is sampleX_sampleY, ...)
#' @slot o the occupancy ratios (per sample) as a named numeric (name is sampleX, ...)
#' @slot num.c,num.o number of o and c parameters
#' @slot names.c,names.o names of the pair to which each c, and of the sample to which each o applies.
#' @slot protein the \code{\link{Protein-class}} object
#' @include ProteinClass.R
setClass("Peptides",
		 representation(
		 	c = "numeric",
		 	o = "list",
		 	num.c = "numeric",
		 	num.o = "numeric",
		 	names.c = "character",
		 	names.o = "list",
		 	protein = "Protein"),
		 validity = check_peptides)

#' Creates a model with concentration and occupancies for a Protein.
#' @param protein a \code{\link{Protein-class}} object
#' @import dplyr
#' @import methods
#' @importFrom stats median
#' @examples
#' data(ENSTest)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- PeptidesModel(ENSTestProtein)
#' @export
PeptidesModel <- function(protein) {
	# Compute the number of occupancies and concentrations
	names.c <- colnames(protein@sample.dependency)
	num.c <- length(names.c)
	c.initial <- rep(0, num.c) # Start from 0
	names(c.initial) <- names.c

	# Take the median ratios as starting point for the c
	median.ratios <- protein@data %>%
		filter(modifications == "") %>%
		filter(pair %in% names.c) %>% # Don't compute other pairs
		group_by(pair) %>%
		summarize(ratio = median(ratio))# %>%
		#slice(match(names.c, pair))

	c.guess <- median.ratios$ratio
	names(c.guess) <- median.ratios$pair

	# Replace the initial values with the estimated one
	c.initial[names(c.guess)] <- c.guess

	# Set the o's to 0.5 for now
	o.initial <- sapply(protein@reference.sample.intersect, function(sample) {
		idx <- protein@data$sample == sample | protein@data$reference == sample
		sample.coverage <- protein@sites.coverage[idx,, drop = FALSE]
		sample.modifications <- colSums(sample.coverage) > 0
		sample.modifications[sample.modifications] <- 0.5
		return(sample.modifications[sample.modifications == 0.5])
	}, simplify = FALSE)
	names.o <- sapply(o.initial, names, simplify = FALSE)

	# How many o do we have?
	num.o <- length(unlist(names.o, use.names = FALSE))

	#o.initial <- matrix(0.5,
	#					nrow = length(protein@reference.sample.intersect),
	#					ncol = length(protein@modifications),
	#					dimnames = list(
	#						sample = protein@reference.sample.intersect,
	#						site = protein@modifications
	#						)
	#					)
	#names.o <- dimnames(o.initial)
	#num.o <- prod(dim(o.initial))
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
