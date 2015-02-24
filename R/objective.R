# Compute the best guess from the model on the ratios
# @param c the concentration
# @param o.ref,o.sample the occupancy ratios in the reference and positive sample
# @param sites.coverage a \code{\link[=make.sites.coverage.matrix]{sites.coverage.matrix}}
# @param sites.activation a \code{\link[=make.sites.activation.matrix]{sites.activation.matrix}}
predicted.ratios <- function(c, o.ref, o.sample, sites.coverage, sites.activation) {
	res <- rep(NA, nrow(sites.coverage))
	for (row in seq(nrow(sites.coverage))) {
		res[row] <- c + sum(sites.coverage[row,] * (log((1 - o.sample) / (1 - o.ref)) + sites.activation[row,] * (log(o.sample / (1 - o.sample)) - log(o.ref / (1 - o.ref)))))
	}
	return(res)
}