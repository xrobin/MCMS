#' The \code{sample.dependency.matrix} defines the dependency between the concentration and samples.
#' It will solve simple circular dependencies and encode the last link as a function of other samples
#' @param sample,reference the aligned vectors of sample names that represent the sample itself and its reference
#' @examples
#' data(ENSTest)
#' make.sample.dependency.matrix(ENSTest$sample, ENSTest$reference)
#' @export
make.sample.dependency.matrix <- function(sample, reference) {
	if (length(sample) != length(reference)) {
		stop("'sample' and 'reference' must have the same length")
	}
	# List all samples that are reference at some point
	ref.as.sample <- intersect(sample, reference)
	# List unique sample-reference sets
	unique.sample.pairs <- unique(data.frame(sample, reference, stringsAsFactors = FALSE))
	# If no intersection, easy
	if (length(ref.as.sample) == 0) {
		m <- diag(nrow(unique.sample.pairs))
		colnames(m) <- rownames(m) <- paste0(unique.sample.pairs$sample, "_", unique.sample.pairs$reference)
	}
	else { # Construct the matrix
		# Can we explain all samples with this overlap?
		samples.with.sampleref <- subset(unique.sample.pairs, reference %in% ref.as.sample)
		samples.with.samples <- subset(unique.sample.pairs, !reference %in% ref.as.sample)
		samples.with.samples.only <- subset(unique.sample.pairs, !reference %in% ref.as.sample & ! sample %in% ref.as.sample)
		# If we can't, bail out for now
		if (!identical(sort(samples.with.samples.only$sample), sort(samples.with.sampleref$sample))) {
			stop("Inconsistence or over-complicated triangular dependencies cannot be solved for now")
		}
		# Support multiple samples appearing in reference only...
		reference.sample <- unique(samples.with.samples.only$reference)
		if (length(reference.sample) > 1) {
			stop("Multiple samples that appear only as reference is not supported yet")
		}
		# Construct the matrix
		m <- matrix(0, nrow = nrow(unique.sample.pairs), ncol = nrow(samples.with.samples))
		colnames(m) <- paste0(samples.with.samples$sample, "_", samples.with.samples$reference)
		rownames(m) <- c(colnames(m), paste0(samples.with.sampleref$sample, "_", samples.with.sampleref$reference))
		diag(m) <- 1
		# List the samples that must be reversed
		# Right now we can do this only
		for (line in seq(nrow(samples.with.sampleref))) {
			x <- samples.with.sampleref[line,]
			linename <- paste0(x["sample"], "_", x["reference"])
			ref1 <- paste0(x["sample"], "_", reference.sample)
			ref2 <- paste0(x["reference"], "_", reference.sample)
			# Find the sample-reference pair
			m[linename, c(ref1)] <- 1
			m[linename, c(ref2)] <- -1
		}
		# Ensure we get either 1 or 0 in the rows (+C_ac - C_ab = 0)
		if (! all(rowSums(m) %in% c(0, 1))) {
			stop("Some rows in the matrix seem incorrect (sum not +1 or -2)")
		}
	}
	return(m)
}

#' The \code{sites.coverage} matrix indicates which observed ratios (rows) cover which modification site (column)
#' @param data the peptide data, as for the \code{\link{Protein}} function
#' @param modifications the list of modifications to be considered
#' @examples
#' data(ENSTest)
#' make.sites.coverage.matrix(ENSTest, unique(ENSTest$modifications))
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom stringr str_replace
#' @export
make.sites.coverage.matrix <- function(data, modifications) {
	# Get the numeric positions
	modification.positions <- as.integer(str_replace(modifications, "[aA-Z]_", ""))

	# Matrix of peptide coverage
	siteAtOrAfterStart <- outer(data$start, modification.positions, FUN = "<=")
	siteAtOrBeforeEnd <- outer(data$end, modification.positions, FUN = ">=")
	sites.coverage <- siteAtOrAfterStart & siteAtOrBeforeEnd

	colnames(sites.coverage) <- modifications
	rownames(sites.coverage) <- str_replace(paste(data$sequence, data$modifications, sep=";"), ";$", "")
	return(sites.coverage + 0) # We want a 0/1 indicator matrix, not TRUE/FALSE
}

#' The \code{sites.activation} indicates which observed ratios (rows) cover have the given modification (column) active
#' @param data the peptide data, as for the \code{\link{Protein}} function
#' @param modifications the list of unique modifications to be considered
#' @examples
#' data(ENSTest)
#' make.sites.activation.matrix(ENSTest, unique(ENSTest$modifications))
#' @importFrom stringr str_detect
#' @export
make.sites.activation.matrix <- function(data, modifications) {
	if (length(modifications) > 0) {
		modification.regexes = paste("^((a|[A-Z]_[0-9]+);)*", modifications, "(;(a|[A-Z]_[0-9]+))*$", sep="")
		mod.on <- outer(data$modifications, modification.regexes, str_detect)
	}
	else {
		mod.on <- matrix(nrow = length(data$peptide), ncol = 0)
	}


	rownames(mod.on) <- str_replace(paste(data$sequence, data$modifications, sep=";"), ";$", "")
	colnames(mod.on) <- modifications
	return(mod.on + 0) # We want a 0/1 indicator matrix, not TRUE/FALSE
}
