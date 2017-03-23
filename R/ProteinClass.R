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

# Make a list of unique modifications by cleaning up the "modifications" column of the data
make_unique_modifications <- function(data, modifications) {
	# Make the list of all the modifications we have
	unique.mods <- unique(unlist(str_split(unique(modifications), ";"))) # And split the multiply-modified peptides
	unique.mods <- unique.mods[unique.mods != ""] # Remove unmodified
	has.acetylations <- any(unique.mods == "a")

	# Where do we have acetylations?
	if (has.acetylations) {
		# Remove acetylation first
		unique.mods <- unique.mods[unique.mods != "a"]
		# Detect where it occurs
		acetylation.pos <- data %>%
			filter(str_detect(modifications, "^([A-Z]_[0-9]+;)*a(;[A-Z]_[0-9]+)*$")) %>%
			select(start) %>%
			distinct(start)
		acetylation.pos <- unique(acetylation.pos$start)
		# Add it to the unique mods
		unique.mods <- c(unique.mods, paste("a", acetylation.pos, sep = "_"))
	}

	return(sort_modifications(unique.mods))
}

# sorts the modifications by site order (a simple sort() would sort by modification type)
sort_modifications <- function(modifications) {
	tmp <- modifications
	if (any(modifications == "a")) {
		warning("Assuming acetylation is at position 0")
		tmp <- str_replace(tmp, "^a$", "a_0") # Assume A is at the very beginning...
	}
	num <- as.numeric(str_split_fixed(tmp, "_", 2)[,2])
	if (any(nas <- is.na(num))) {
		stop(paste0("Modifications with no positional information can't be sorted: ", paste0(tmp[nas], collapse = ", ")))
	}
	return(modifications[order(num)])
}


#' @name Protein-class
#' @rdname Protein-class
#' @exportClass Protein
#' @title A class that describe a protein. Use the \code{\link{Protein}} function for easy object creation
#' @slot data a data.frame with the data for one protein.
#' @slot modifications vector listing all the modification sites
#' @slot samples,references vector with the name of all the samples, reference or sample (unique and sorted)
#' @slot reference.sample.overlap those sample that appear both as reference and as sample so one can "close the loop"
#' @slot reference.sample.intersect all the sample and reference names. Will for instance correspond to all the occupancy ratios that need to be computed.
#' @slot sample.dependency a matrix with all the sample_reference pairs in rows, and the columns over all non-redundant pairs. Useful to compute concentration ratios
#' @slot site.coverage indicator matrix of which peptide covers which site
#' @slot site.on.off indicator matrix of which sites are on or off. Meaningful only on the indices where \code{site.coverage} is 1
setClass("Protein",
		 representation(
		 	data = "data.frame",
		 	modifications = "character",
		 	samples = "character", references = "character",
		 	reference.sample.overlap = "character",
		 	reference.sample.intersect = "character",
			sample.dependency = "matrix",
			sites.coverage = "matrix",
			sites.activation = "matrix"
			),
		 validity = check_protein)


#' End-user function to create a Protein object
#' @param data the \code{data.frame}
#' @param sample.dependency.matrix a matrix encoding how samples (rows) depend on the parameters (columns). It can be omitted for simple problems with no complex dependencies.
#' @param decimate.sample.dependency.matrix whether to remove samples that are not observed and the parameters that are not
#' @param na.action either \code{\link{na.fail}} or \code{\link{na.omit}}. Other actions
#' @description
#' Only acetylation (a) is supported as a non-positional modification. It is added to the N-term of the peptide
#' @import methods
#' @examples
#' data(ENSTest)
#' Protein(ENSTest)
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_replace
#' @export
Protein <- function(data, sample.dependency.matrix, decimate.sample.dependency.matrix = TRUE, na.action = na.fail) {

	### Remove NAs in data
	if (any(which.ones <- is.na(data$q) & data$n == 1)) {
		warning(sprintf("Filling %d missing q values with 0 where n == 1", sum(which.ones)))
		data$q[which.ones] <- 0
	}
	# Process the na.action
	data <- na.action(data)

	modifications <- make_unique_modifications(data, data$modifications)

	samples <-  sort(unique(data$sample))
	references <- sort(unique(data$reference))
	reference.sample.overlap = sort(intersect(samples, references))
	reference.sample.intersect <- sort(unique(c(samples, references)))

	data$pair <- paste(data$sample, data$reference, sep = "_")

	protein <- new("Protein",
		data = data,
		modifications = modifications,
		samples = samples, references = references,
		reference.sample.overlap = reference.sample.overlap,
		reference.sample.intersect = reference.sample.intersect
		)

	if (missing(sample.dependency.matrix)) {
		protein@sample.dependency <- make.sample.dependency.matrix(data$sample, data$reference)
	}
	else {
		if (decimate.sample.dependency.matrix) {
			if (any(which <- ! data$pair %in% rownames(sample.dependency.matrix))) {
				stop(paste0("Some pairs were not encoded in the dependency matrix: ", paste(data$pair[which], collapse = ", ")))
			}
			# Remove rows that do not appear in a data$pair:
			sample.dependency.matrix <- sample.dependency.matrix[rownames(sample.dependency.matrix) %in% data$pair,, drop=FALSE]
			# Remove columns that have only 0s (= parameter is not used)
			sample.dependency.matrix <- sample.dependency.matrix[, colSums(abs(sample.dependency.matrix)) > 0, drop = FALSE]
		}
		protein@sample.dependency <- sample.dependency.matrix
	}

	protein@sites.coverage <- make.sites.coverage.matrix(data, modifications)
	protein@sites.activation <- make.sites.activation.matrix(data, modifications)

	return(protein)
}
