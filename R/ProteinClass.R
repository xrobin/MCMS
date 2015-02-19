check_protein <- function(object) {
	data.columns <- c("protein", "sequence", "modifications", "start", "end", "length", "sample", "reference", "ratio", "q", "n")

	errors <- character()
	if (any(missing.cols <- ! data.columns %in% colnames(object@data))) {
		msg <- paste0("Missing columns in data: ", paste(data.columns[missing.cols], collapse = ", "))
		errors <- c(errors, msg)
	}

	if (length(unique(object@data$protein)) != 1) {
		msg <- paste0("Only one protein is expected at a time, not ", paste(unique(object@data$protein), collapse = ", "))
		errors <- c(errors, msg)
	}

	if (!identical(sort(unique(object@data$sample)), object@samples)) {
		msg <- "Inconsistency in the samples"
		errors <- c(errors, msg)
	}

	if (!identical(sort(unique(object@data$reference)), object@references)) {
		msg <- "Inconsistency in the references"
		errors <- c(errors, msg)
	}

	if (length(errors) == 0) TRUE else errors
}


#' @slot data a data.frame with the data for one protein.
#' @slot modifications vector listing all the modification sites
#' @slot samples,references vector with the name of all the samples, reference or sample (unique and sorted)
#' @slot reference.sample.overlap those sample that appear both as reference and as sample so one can "close the loop"
#' @slot sample.dependency a matrix with all the sample_reference pairs in rows, and the columns over all non-redundant pairs. Useful to compute concentration ratios
#' @slot site.coverage indicator matrix of which peptide covers which site
#' @slot site.on.off indicator matrix of which sites are on or off. Meaningful only on the indices where \code{site.coverage} is 1
setClass("Protein",
		 representation(
		 	data = "data.frame",
		 	modifications = "character",
		 	samples = "character", references = "character",
		 	reference.sample.overlap = "character",
			sample.dependency = "matrix",
			sites.coverage = "matrix",
			sites.activation = "matrix"
			),
		 validity = check_protein)


#' End-user function to create a Protein object
#' @param data the \code{data.frame}
#' @description
#' Only acetylation (a) is supported as a non-positional modification. It is added to the N-term of the peptide
#' @examples
#' data(ENSTest)
#' Protein(ENSTest)
#' @export
Protein <- function(data) {
	browser()
	modifications <- unique(unlist(str_split(unique(data$modifications), ";")))
	modifications <- modifications[modifications != ""]

	samples <-  sort(unique(data$sample))
	references <- sort(unique(data$reference))
	reference.sample.overlap = sort(intersect(samples, references))

	protein <- new("Protein",
		data = data,
		modifications = modifications,
		samples = samples, references = references,
		reference.sample.overlap = reference.sample.overlap
		)

	protein@sample.dependency = make.sample.dependency.matrix(protein@samples, protein@references),
	protein@sites.coverage = make.sites.coverage.matrix()
	protein@sites.activation <- make.sites.activation.matrix()
}