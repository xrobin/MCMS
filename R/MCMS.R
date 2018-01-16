#' The main function here.
#' Receives the data and the variance model, prepares all the elements and then performs the simulation
#' @param Peptides the Peptides model data
#' @param var.model the variance model
#' @param n number of iterations
#' @param n.out number of iterations to record
#' @param burn.in the burn-in time (integer, NOT a proportion)
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
#' \dontrun{
#' # In a shell:
#' R -e "library(MCMS); m <- MCMS(ENSTestModel, var.model, 10000000, 10000, 300000); \
#' summary(as.data.frame(m))";
#' }
#' @export
#' @useDynLib MCMS
MCMS <- function(PeptidesModel, var.model, n, n.out, burn.in,
				 c.prior.sd = sqrt(2), o.prior.shape1 = .5, o.prior.shape2 = .5,
				 prior_move_proportion = .02, c_sd = 0.05,
				 o_sd = 0.05, o_k_scale = 1/100, o_restrict = 1e-3,
				 seed = round(runif(10, -.Machine$integer.max, .Machine$integer.max)),
				 verbose = FALSE, cooling = TRUE) {
	run_MCMC_Cpp(PeptidesModel, var.model, n, n.out, burn.in, c_prior_sd = c.prior.sd, o_prior_shape1 = o.prior.shape1,
				 o_prior_shape2 = o.prior.shape2, o_restrict = o_restrict,
				 prior_move_proportion = prior_move_proportion, c_sd = c_sd, o_sd = o_sd, o_k_scale = o_k_scale,
				 seed = seed, verbose = verbose, cooling = cooling)
}