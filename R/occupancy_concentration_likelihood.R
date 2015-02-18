#' Compute the best guess from the model on the ratios
#' @param c the 
#' @param o,o.prime the occupancy ratios, as matrices of the same size as Iis and tis.
#' @keywords internal
mu.coo <- function(c, o, o.prime, Iis, tis, npeptides, nsites) {
	c + rowSums( Iis * (log((1 - o.prime) / (1 - o)) + tis * (log(o.prime / (1 - o.prime)) - log(o / (1 - o)))))
}

#' K_i
#' @keywords internal
K.coo <- function(gammas, xbar, mu) {
	1 + gammas * (xbar - mu) ^2
}

#' l;i
#' @param xbar observed ratios
#' @param mu predicted ratios
#' @keywords internal
li.coo <- function(xbar, mu, gammas, nu.tildes, K = K.coo(gammas, xbar, mu)) {
	(2 * gammas * nu.tildes) / K * (xbar - mu)
}

#' l;ii
#' @keywords internal
lii.coo <- function(xbar, mu, gammas, nu.tildes, K = K.coo(gammas, xbar, mu)) {
	(2 * gammas * nu.tildes) / K *  (2 * gammas / K * (xbar - mu)^2 - 1)
}

#' mu_i;c
#' @keywords internal
muic.coo <- function(c) {
	return(1)
}

#' mu_i;s'
#' @keywords internal
muis.prime.coo <- function(o.prime, Iis, tis) {
	Iis * (-1 / (1 - o.prime) + tis * ( 1 / o.prime + 1 / (1 - o.prime)))
}

#' mu_i;s'
#' @keywords internal
muis.coo <- function(o, Iis, tis) {
	Iis * (1 / (1 - o) - tis * ( 1 / o + 1 / (1 - o)))
}

#' mu_i;s's'
#' @keywords internal
muis.prime.2.coo <- function(o.prime, Iis, tis) {
	Iis * (-1 / (1 - o.prime)^2 + tis * ( -1 / o.prime^2 + 1 / (1 - o.prime)^2))
}
#' mu_i;ss
#' @keywords internal
muis.2.coo <- function(o, Iis, tis) {
	Iis * (1 / (1 - o)^2 - tis * ( -1 / o^2 + 1 / (1 - o)^2))
}


#' Objective function (log likelihood) and gradient functions for the computation of the concentration and occupancy ratios
#' @param coo the parameters: c, o and o'.
#' @param x the measured data, in particular the ratios in column \code{x}, and the \code{n} and \code{q} measures
#' @param ratios the measured ratios
#' @param gammas,nu.tildes the bayesian estimates of the ratios precisions. TODO: change the function so that it uses the proposed ratio precisions instead!
#' @param Iis the I_is a logical matrix that maps which peptide (rows) spans which modification sites (columns)
#' @param tis the t_is a logical matrix that define whether the modification sites (columns) are active or not for a given peptide (row)
#' @param on.observed,off.observed counts of the # of times the site has been seen on and off - for the beta prior
#' @param npeptides,nsites the number of peptides and sites under consideration
#' @param use.prior the weight of the prior, typically 0 (ignore prior) or 1 (use prior). 
#' @export
fit.concentration.occupancy.objective <- function(coo, data, model, var.model,
												  Iis, tis, 
												  #on.observed, off.observed,
												  npeptides, nsites
												  ) {
	c <- coo[1]
	o <- matrix(coo[seq_len(nsites) + 1], nrow(tis), ncol(tis), byrow = TRUE)
	o.prime <- matrix(coo[seq_len(nsites) + nsites + 1], nrow(tis), ncol(tis), byrow = TRUE)
	
	# Getting new ratios from the occupancy ratios
	current.mu <- mu.coo(c, o, o.prime, Iis, tis, npeptides, nsites)
	# Put it into the "model" data.frame and re-compute the bayesian variance estimate if possible
	if (missing(var.model)) {
		if (missing(model)) {
			warning("Re-using variance from data instead of computing it from the model")
			model <- data
			model$x <- current.mu
		}
		else { # Assume 'model' contains the updated variance already - but replace x just to be sure
			model$x <- current.mu
		}
	} else if (missing(model)) {
		model <- data
		model$x <- current.mu
		model <- addBayesPrec(var.model, model, mean="x", q = "q", n = "n")
	} else {
		model$x <- current.mu
		model <- addBayesPrec(var.model, model, mean="x", q = "q", n = "n")
	}

	# Rating the fit
	sum(- model$nu.tilde * log(K.coo(model$gamma, data$x, current.mu)) + log(model$gamma) / 2) + 
			#Beta prior
			#sum(log(dbeta(o[1,], .5+on.observed*0, .5+off.observed*0))) +
			#sum(log(dbeta(o.prime[1,], .5+on.observed*0, .5+off.observed*0))) +
			#Beta prior without observations
			sum(log(dbeta(o[1,], beta.shape1, beta.shape2))) +
			sum(log(dbeta(o.prime[1,], beta.shape1, beta.shape2))) +
			# Prior on the concentration
			(- abs(c) * 2)
}
