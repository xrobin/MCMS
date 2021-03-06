#' Compute the likelihood of a Peptides model
#' @param Peptides the Peptides model data
#' @param var.model the variance model
#' @param c,o values of the parameters for which to estimate the likelihood
#' @param c.prior.sd,o.prior.shape1,o.prior.shape2 parameters for the prior
#' @param verbose if TRUE, will print internal details of the Monte Carlo object (elements and pointers)
#' @examples
#' data(ENSTest, var.model)
#' ENSTestProtein <- Protein(ENSTest)
#' ENSTestModel <- Peptides(ENSTestProtein)
#' getLikelihood(ENSTestModel, var.model)
#' # With custom c
#' someC <- ENSTestModel@c
#' someC[2] <- 15
#' getLikelihood(ENSTestModel, var.model, c = someC)
#' @export
#' @useDynLib MCMS
getLikelihood <- function(Peptides, var.model, c = Peptides@c, o = Peptides@o,
						  c.prior.sd = sqrt(2), o.prior.shape1 = .5, o.prior.shape2 = .5,
						  verbose = FALSE) {
	getLikelihood_MCMC_Cpp(Peptides, var.model, c, o,
						   c_prior_sd = c.prior.sd, o_prior_shape1 = o.prior.shape1, o_prior_shape2 = o.prior.shape2,
						   verbose = verbose)
}