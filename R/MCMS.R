#' The main function here.
#' Receives the data and the variance model, prepares all the elements and then performs the simulation
#' @param PeptidesModel the PeptidesModel data
#' @param var.model the variance model
#' @param n number of iterations
#' @param n.out number of iterations to record
#' @param burn.in the burn-in time (integer, NOT a proportion)
#' @param c.prior.sd,o.prior.shape1,o.prior.shape2 parameters for the prior
#' @param prior_move_proportion the proportion of moves that are sampled from the prior
#' @param c_sd,o_sd standard deviations for the moves on o and c
#' @param o_k_scale scaling factor for the moves on o when they are close to 0 or 1
#' @param o_restrict reject any move that is closer than o_restrict to 0 or 1
#' @examples
#' data(ENSTest)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- PeptidesModel(ENSTestProtein)
#' MCMS(ENSTestModel, var.model, 10000, 1000
#' 	scale = 1, shape1 = .5, shape2 = .5,
#' 	prior_move_proportion = .02, c_sd = 0.05, o_sd = 0.05, o_k_scale = 1/100)
#' \dontrun{
#' # In a shell:
#' # R -e "library(MCMS); m <- MCMS(ENSTestModel, var.model, 10000000, 10000, 300000); summary(as.data.frame(m))";
#' }
#' @export
#' @useDynLib MCMS
MCMS <- function(PeptidesModel, var.model, n, n.out, burn.in,
				 c.prior.sd = sqrt(2), o.prior.shape1 = .5, o.prior.shape2 = .5,
				 prior_move_proportion = .02, c_sd = 0.05,
				 o_sd = 0.05, o_k_scale = 1/100, o_restrict = 1e-3) {
	run_MCMC_Cpp(PeptidesModel, var.model, n, n.out, burn.in, c_prior_sd = c.prior.sd, o_prior_shape1 = o.prior.shape1,
				 o_prior_shape2 = o.prior.shape2, o_restrict = o_restrict,
				 prior_move_proportion = prior_move_proportion, c_sd = c_sd, o_sd = o_sd, o_k_scale = o_k_scale)
}