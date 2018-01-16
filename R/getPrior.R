#' Compute the prior of a set of parameters
#' @param Peptides the Peptides model data
#' @param var.model the variance model
#' @param c,o values of the parameters for which to estimate the prior
#' @param c.prior.sd,o.prior.shape1,o.prior.shape2 parameters for the prior
#' @param verbose if TRUE, will print internal details of the Monte Carlo object (elements and pointers)
#' @examples
#' data(ENSTest, var.model)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- Peptides(ENSTestProtein)
#' getPrior(ENSTestModel)
#' # With custom c
#' someC <- ENSTestModel@@c
#' someC[] <- 0
#' getPrior(ENSTestModel, var.model, c = someC, o = ENSTestModel@@o)
#' @export
#' @useDynLib MCMS
getPrior <- function(Peptides, var.model, c = Peptides@c, o = Peptides@o,
						  c.prior.sd = sqrt(2), o.prior.shape1 = .5, o.prior.shape2 = .5,
						  verbose = FALSE) {
	getPrior_MCMC_Cpp(Peptides, var.model, c, o,
						   c_prior_sd = c.prior.sd, o_prior_shape1 = o.prior.shape1, o_prior_shape2 = o.prior.shape2,
						   verbose = verbose)
}