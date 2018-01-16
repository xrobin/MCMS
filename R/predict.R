#' Predicts ratios from a Peptides model
#' @param object the \code{\link{Peptides}} model to predict from
#' @param newdata an optional \code{\link{Protein}} model to use as data for the predictions
#' @param ... ignored.
#' @examples
#' data(ENSTest)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- Peptides(ENSTestProtein)
#' predict(ENSTestModel)
#' @importFrom stats predict
#' @include PeptidesClass.R
#' @export
predict.Peptides <- function(object, newdata = object@protein, ...) {
	if (!is(newdata, "Protein"))
		stop("Protein object expected as newdata")

	# Get a c for all the pairs
	all.c <- tcrossprod(object@c, newdata@sample.dependency)

	# Copy the data
	data <- newdata@data
	# Get the relevant c's and o's
	c <- all.c[1, match(data$pair, colnames(all.c))]
	#o.ref <- object@o[match(data$reference, object@names.o)]
	#o.sample <- object@o[match(data$sample, object@names.o)]

	all.pairs <- unique(data$pair)
	mu <- rep(NA, nrow(data))
	for (current.pair in all.pairs) {
		print(current.pair)
		row.numbers <- which(data$pair == current.pair)
		data.rows <- data[row.numbers,]

		o.ref <- object@o[data.rows[1, 'reference'],]
		o.sample <- object@o[data.rows[1, 'sample'],]

		mu[row.numbers] <- predicted.ratios(c[current.pair], o.ref, o.sample, newdata@sites.coverage[row.numbers,, drop=FALSE], newdata@sites.activation[row.numbers,, drop=FALSE])
	}

	return(mu)
}

#setMethod("predict", signature(object = "Peptides"), predict.Peptides)