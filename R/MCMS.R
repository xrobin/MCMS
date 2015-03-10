#' The main function here.
#' Receives the data and the variance model, prepares all the elements and then performs the simulation
#' @param ENSTestModel
#' @param var.model the variance model
#' @param n number of iterations
#' @param number of iterations to record
#' @examples
#' data(ENSTest)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- PeptidesModel(ENSTestProtein)
#' MCMS(ENSTestModel, var.model, 10000, 1000
#' 	scale = 1, shape1 = .5, shape2 = .5,
#' 	prior_move_proportion = .02, c_sd = 0.05, o_sd = 0.05, o_k_scale = 1/100)
#' @export
#' @useDynLib MCMS
MCMS <- function(ENSTestModel, var.model, n, n.out,
				 scale = 1, shape1 = .5, shape2 = .5,
				 prior_move_proportion = .02, c_sd = 0.05,
				 o_sd = 0.05, o_k_scale = 1/100, o_restrict = 1e-5) {
	run_MCMC_Cpp(ENSTestModel, var.model, i, scale = scale, shape1 = shape1, shape2 = shape2,
			prior_move_proportion = prior_move_proportion, c_sd = c_sd, o_sd = o_sd, o_k_scale = o_k_scale, o_restrict = o_restrict)
}