.onLoad <- function(libname, pkgname) {
	op <- options()
	op.MCMS <- list(
		MCMS.beta.shape1 = 0.5,
		MCMS.beta.shape2 = 0.5
	)
	toset <- !(names(op.MCMS) %in% names(op))
	if(any(toset)) options(op.MCMS[toset])

	invisible()
}