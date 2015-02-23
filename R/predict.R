#' Predicts ratios from a Peptides model
#' @examples
#' data(ENSTest)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- PeptidesModel(ENSTestProtein)
#' predict(ENSTestModel)
#' @importFrom stats predict
#' @include PeptidesClass.R
#' @export
predict.Peptides <- function(object, newdata = object@protein) {
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


	for (i in 1:seq_along(newdata@data$ratio)) {
		row <- newdata@data[i, ]
		# Get o.ref and o.sample
		o.ref <- object@o[row$reference,]
		o.sample <- object@o[row$sample,]

	}
	browser()
}

#setMethod("predict", signature(object = "Peptides"), predict.Peptides)