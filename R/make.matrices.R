#' The \code{sample.dependency.matrix} defines the dependency between the concentration and samples.
#' It will solve simple circular dependencies and encode the last link as a function of other samples
#' @param data the MS data, specifically columns \code{}
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
	unique.sample.pairs <- unique(data.frame(sample, reference))
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
			m[linename, c(ref1, ref2)] <- -1
		}
		# Ensure we get either 1 or -2 in the rows
		if (! all(rowSums(m) %in% c(-2, 1))) {
			stop("Some rows in the matrix seem incorrect (sum not +1 or -2)")
		}
	}

	return(m)
}