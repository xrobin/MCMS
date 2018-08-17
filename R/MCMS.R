#' The main function here.
#' Receives the data and the variance model, prepares all the elements and then performs the simulation
#' @param Peptides the Peptides model data
#' @param var.model the variance model
#' @param n number of iterations
#' @param n.out number of iterations to record. Should be a divisor of (1 - burn.in) * n for best results.
#' @param burn.in the fraction of the time \code{n} to use as burn-in (a proportion, not integer time)
#' @param max.repeats the maximum number of repeats of 1E9 iterations to run
#' @param c.prior.sd,o.prior.shape1,o.prior.shape2 parameters for the prior
#' @param prior_move_proportion the proportion of moves that are sampled from the prior
#' @param c_sd,o_sd standard deviations for the moves on o and c
#' @param o_k_scale scaling factor for the moves on o when they are close to 0 or 1
#' @param o_restrict reject any move that is closer than o_restrict to 0 or 1
#' @param verbose if TRUE, will print internal details of the Monte Carlo object (elements and pointers)
#' @param cooling whether to perform cooling during the MonteCarlo. If \code{TRUE}, the simulation will start with
#' @param seed a seed from the sampling, taken from \code{\link{runif}} by default.
#' @importFrom stats runif
#' @examples
#' data(ENSTest)
#' data(var.model)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- Peptides(ENSTestProtein)
#' MCMS(ENSTestModel, var.model,
#' 	prior_move_proportion = .02, c_sd = 0.05, o_sd = 0.05, o_k_scale = 1/100)
#' @export
#' @useDynLib MCMS
MCMS <- function(Peptides, var.model, n = calcIterations(Peptides), n.out = 7E3, burn.in = 0.3, max.repeats = 20,
				 c.prior.sd = sqrt(2), o.prior.shape1 = .5, o.prior.shape2 = .5,
				 prior_move_proportion = .02, c_sd = 0.05,
				 o_sd = 0.05, o_k_scale = 1/100, o_restrict = 1e-3,
				 seed = round(runif(10, -.Machine$integer.max, .Machine$integer.max)),
				 verbose = FALSE, cooling = TRUE) {

	n.params <- nParams(Peptides)

	if (n > 1E9) {
		# Do 1E9 10 times with less output...
		repeats <- min(ceiling(n / 1E9), 20) # max 10 repeats
		n.out <- n.out / repeats # 1E4 - 30%
		n <- 1E9
	}
	else {
		repeats <- 1
	}

	if (burn.in > 1) {
		stop("burn.in should be a proportion, not the time itself.")
	}

	mc.result <- sapply(seq_len(repeats), function(i) {
		if (verbose) {
			message(sprintf("Doing repeat %d/%d", i, repeats))
		}

		run_MCMC_Cpp(Peptides, var.model, n, n.out, burn.in * n, c_prior_sd = c.prior.sd, o_prior_shape1 = o.prior.shape1,
					 o_prior_shape2 = o.prior.shape2, o_restrict = o_restrict,
					 prior_move_proportion = prior_move_proportion, c_sd = c_sd, o_sd = o_sd, o_k_scale = o_k_scale,
					 seed = seed, verbose = verbose, cooling = cooling)
	}, simplify = FALSE)
	mc.result <- do.call(rbind, mc.result)

	if (! requireNamespace("coda", quietly = TRUE)) {
		stop("coda package unavailable, effictive sample size not available")
	}
	else {
		attr(mc.result, "effectiveSize") <-
			coda::effectiveSize(coda::mcmc(mc.result[,seq_len(n.params)], start = (n * burn.in) + 1, end = n * repeats, thin = (n * repeats - (n * burn.in) - 1) / (n.out * repeats)))
	}

	return(mc.result)
}
