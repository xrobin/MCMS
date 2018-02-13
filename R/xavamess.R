#' Remove missing values with a warning
#' This function behaves exactly like \code{\link{na.omit}}, except that it prints a warning if any missing value was removed
#' @param object an R object
#' @param ... additional arguments
#' @importFrom stats na.omit
#' @examples
#' DF <- data.frame(x = c(1, 2, 3), y = c(0, 10, NA))
#' na.warn(DF)
#' @export
na.warn <- function(object, ...) {
	object <- na.omit(object, ...)
	if (!is.null(missing <- attr(object, "na.action"))) {
		warning(sprintf("removed %d missing values in object", length(missing)))
	}
	return(object)
}


#' Maps with a left_join, ensuring constant number of rows
#'
#' This avoids generating left duplicates if entries are duplicated in the right \code{by}.
#' If duplicates are generated, the function will stop with an error.
#' @param x the \code{\link{data.frame}} with fixed size
#' @param map a \code{\link{data.frame}} with the new columns.
#' @param by see \code{\link{left_join}}
#' @importFrom dplyr left_join
#' @examples
#' a <- data.frame(A = LETTERS[1:5], B = 1:5)
#' b <- data.frame(A = LETTERS[1:5], C = 15:11)
#' dplyr::left_join(a, b)
#' safe.mapping(a, b)
#'
#' # Now with duplicates
#' \dontrun{
#' c <- data.frame(A = LETTERS[1:5], C = 20:11)
#' dplyr::left_join(a, c) # creates duplicates of x
#' safe.mapping(a, c) # stops with an error
#' }
#'
#' @export
safe.mapping <- function(x, map, by = NULL) {
	nrow.orig <- nrow(x) # Keep nrow before mapping
	x <- left_join(x, map, by = by)
	if (nrow.orig != nrow(x)) {
		stop(sprintf("More rows after mapping: %d -> %d", nrow.orig, nrow(x)))
	}
	return(x)
}
