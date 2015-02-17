#' Generates a new coo by sampling from the prior
#' @param coo the previous coo
#' @param c.sd,o.sd standard deviations for the normal distribution
#' @param beta.shape1, beta.shape2 the shape parameters for the beta distribution (prior)
#' @return a list with 3 elements:
#' \itemize{
#' 	\item coo the new coo
#' 	\item log.previous.bias the bias of the new coo within the previous prior
#' 	\item log.new.bias the bias of the previous coo within the new prior
#' }
draw.new.coo.from.prior <- function(coo, c.sd = 0.05, o.sd = 0.05, beta.shape1 = getOption("MCMS.beta.shape1"), beta.shape2 = getOption("MCMS.beta.shape2")) {
	which.changes <- floor(runif(1, min=1, max = length(coo) + 1))
	if (which.changes == 1) {
		# Draw concentration from double laplace distribution
		log.new.bias <- log(1/2*dexp(abs(coo[1]), rate = 2))
		coo[1] <- rexp(1, rate = 2) * sign(round(runif(1)) - .5)
		log.previous.bias <- log(1/2*dexp(abs(coo[1]), rate = 2))
	}
	else {
		log.new.bias <- log(dbeta(coo[which.changes], beta.shape1, beta.shape2))
		coo[which.changes] <- rbeta(1, beta.shape1, beta.shape2)
		log.previous.bias <- log(dbeta(coo[which.changes], beta.shape1, beta.shape2))
	}
	return(list(coo = coo, log.previous.bias = log.previous.bias, log.new.bias = log.new.bias))
}

#' Generates a new coo by a random normal number
#' @param coo the previous coo
#' @param c.sd,o.sd standard deviations for the normal distribution
#' @param k.scale a scaling factor that will reduce the changes on o and o' when close to 0 or 1
#' @param beta.shape1, beta.shape2 the shape parameters for the beta distribution (prior)
#' @param o.range accept o and o' only within that range (non inclusive)
#' \itemize{
#' 	\item coo the new coo
#' 	\item log.previous.bias the bias of the new coo within the previous prior
#' 	\item log.new.bias the bias of the previous coo within the new prior
#' }
draw.new.coo <- function(coo, c.sd = 0.05, o.sd = 0.05, k.scale = 1/100, beta.shape1 = getOption("MCMS.beta.shape1"), beta.shape2 = getOption("MCMS.beta.shape2"), o.range = c(1e-5, 1-1e-5)) {
	log.previous.bias <- log.new.bias <- 0
	which.changes <- floor(runif(1, min=1, max = length(coo) + 1))
	if (which.changes == 1) {
		coo[1] <- coo[1] + rnorm(1, mean=0, sd = c.sd)
	}
# 	else if (which.changes <= length(coo)) {
	else {
		# Make the move smaller for o close to 0 or 1
		previous.o <- coo[which.changes]
		k <- dbeta(previous.o, beta.shape1, beta.shape2) * (1 / (2 * (1 - previous.o)) - 1 / (2 * previous.o)) * k.scale
		sigma <- 1 / (abs(k) + 1 / o.sd)
		new.o <- previous.o + rnorm(1, mean = 0, sd = sigma)
		if (new.o > o.range[1] && new.o < o.range[2]) {
			coo[which.changes] <- new.o
			log.previous.bias <- log(dnorm(new.o, mean = previous.o, sd = sigma))
			new.k <- dbeta(new.o, beta.shape1, beta.shape2) * (1 / (2 * (1 - new.o)) - 1 / (2 * new.o)) * k.scale
			new.sigma <- 1 / (abs(new.k) + 1 / o.sd)
			log.new.bias <- log(dnorm(previous.o, mean = new.o, sd = new.sigma))
		}
	}
# 	else { # Move both o and o'
# 		move <- rnorm(1, mean=0, sd = .05)
# 		idx.o <- which.changes - length(coo)
# 		idx.o.prime <- which.changes - length(coo) + nsites
# 		temp.o <- coo[idx.o] + move
# 		temp.o.prime <- coo[idx.o.prime] + move
# 		if (temp.o > 0 && temp.o.prime > 0 && temp.o < 1 && temp.o.prime < 1) {
# 			coo[idx.o] <- temp.o
# 			coo[idx.o.prime] <- temp.o.prime
# 		}
# 	}
	return(list(coo = coo, log.previous.bias = log.previous.bias, log.new.bias = log.new.bias))
}

#' Estimate occupancy ratios and uncertainties with Monte Carlo sampling
#' @param coo.initial initial guess on c, o and o.prime
#' @param x the peptides data.frame with gammas and nu.tildes
#' @param var.model the variance model (from \code{globalBayesVar})
#' @param Iis matrix informing whether a peptide spans the modification site
#' @param tis matrix informing whether the peptide has the modification on
#' @param steps number of monte carlo steps (per parameter)
#' @param burn.in.time proporition of \code{steps} that is used for burn-in
#' @param prior.move.proportion proportion of non-standard moves based on the prior
#' @param npeptides,nsites number of peptide and phosphosites in the estimation
#' @param beta.shape1,beta.shape2 shapes for the beta prior on occupancy ratios
#' @param on.observed,off.observed prior on the number of observations
#' @return a matrix with ncol as the dimension of the problem + 1. Each row is one step of the simulation.
#' The first column contains the c, the last the likelihood, and the middle ones the occupancy ratios corresponding to the coo.
mc.occupancy <- function(coo.initial, x, var.model, Iis, tis, steps = 100000, burn.in.proportion = 0.3, prior.move.proportion = .02, npeptides = npeptides, nsites = nsites,
						 beta.shape1 = getOption("MCMS.beta.shape1"), beta.shape2 = getOption("MCMS.beta.shape2")
						 #on.observed, off.observed
						 ) {
	print(unique(x$Longest.Isoform))
	print(length(coo.initial))

	burn.in.time <- ceiling(burn.in.proportion * steps)
	allocate <- try(mc.result <- matrix(NA, nrow=steps - burn.in.time, ncol=length(coo.initial) + 1))
	if (is(allocate, "try-error")) {
		print("Cannot allocate the results.")
		print(allocate)
		return(NULL)
	}
	# Having model as an already shaped data.frame might speed up the simulation... or not
	model <- x

	previous.coo <- coo.initial
	previous.likelihood <-
		fit.concentration.occupancy.objective(coo.initial, x, model = model, var.model,
											  Iis = Iis, tis = tis, #on.observed = on.observed, off.observed = off.observed,
											  npeptides = npeptides, nsites = nsites)

	for (t in 1:steps) {
		log.previous.bias <- 0
		log.new.bias <- 0
		if (runif(1) > (1 - prior.move.proportion)) { # Resample from the prior every 50 steps
			new.coo.and.biases <- draw.new.coo.from.prior(previous.coo)
			new.coo <- new.coo.and.biases[["coo"]]
			log.previous.bias <- new.coo.and.biases[["log.previous.bias"]]
			log.new.bias <- new.coo.and.biases[["log.new.bias"]]
			new.likelihood <-
				fit.concentration.occupancy.objective(new.coo, data = x, model = model, var.model = var.model,
													  Iis = Iis, tis = tis, #on.observed = on.observed, off.observed = off.observed,
													  npeptides = npeptides, nsites = nsites)
		}
		else {
			new.coo.and.biases <- draw.new.coo(previous.coo)
			new.coo <- new.coo.and.biases[["coo"]]
			log.previous.bias <- new.coo.and.biases[["log.previous.bias"]]
			log.new.bias <- new.coo.and.biases[["log.new.bias"]]
			new.likelihood <-
				fit.concentration.occupancy.objective(new.coo, data = x, model = model, var.model = var.model,
													  Iis = Iis, tis = tis, #on.observed = on.observed, off.observed = off.observed,
													  npeptides = npeptides, nsites = nsites)
		}
		if (is.na(new.likelihood)) {
			print("NA in likelihood")
			browser()
		}
		nominator <- new.likelihood + log.new.bias
		denominator <- previous.likelihood + log.previous.bias
		if (nominator > denominator || runif(1) < exp(nominator - denominator)) {
		#if (new.likelihood + log.previous.bias > previous.likelihood + log.new.bias || runif(1) < exp(new.likelihood - previous.likelihood+ log.previous.bias - log.new.bias)) {
			previous.coo <- new.coo
			previous.likelihood <- new.likelihood
		}
		if (t > burn.in.time) {
			mc.result[t - burn.in.time,] <- c(previous.coo, previous.likelihood)
		}
	}
	return(mc.result)
}

library(matrixStats)

#' Reads the monte carlo output.
#' @param dir the folder containing the output .RData files
#' @param burn.in.proportion the proportion of the MC simulation to remove. Should be 0 if it was removed before saving
#' @importFrom matrixStats colSds
#' @importFrom stringr str_match
#' @importFrom stringr str_replace
read.mc.results <- function(dir = mc.dir, burn.in.proportion = 0.1) {
	output.files <- list.files(dir, "\\.RData$", full.names = TRUE)
	ENSPs <- str_match(output.files, ".+(ENSP\\d+)\\.RData")[, 2]
	names(output.files) <- ENSPs

	#c.list <- c.sd.list <- rep(NA, length(output.files))
	#p.list <- pbar.list <- p.sd.list <- pbar.sd.list <- vector("list", length(output.files))
	#names(c.list) <- names(c.sd.list) <- names(p.list) <- names(pbar.list) <- names(p.sd.list) <- names(pbar.sd.list) <- ENSPs
	all.data <- sapply(ENSPs, function(ENSP) {
		print(ENSP)
		file <- output.files[ENSP]
		loaded <- try(load(file))
		if (! is(loaded, "try-error")) {
			nsteps <- nrow(mc.result)
			nsites <- (ncol(mc.result) - 2) / 2
			#browser()
			#ENSP <- str_match(file, ".+/(ENSP\\d+)\\.RData")[1, 2]
			#burnt.in <- tail(mc.result, n = (1 - burn.in.proportion) * nsteps)
			burnt.in <- mc.result[ -seq(1, burn.in.proportion * nsteps),]
			c <- mean(burnt.in[, 1])
			c.sd <- sd(burnt.in[, 1])
			#c.list[ENSP] <- c
			#c.sd.list[ENSP] <- sd(burnt.in[, 1])

			if (nsites > 0) {
				o.cols <- (1:nsites) + 1
				o.prime.cols <- o.cols + nsites
				o <- burnt.in[, o.cols, drop=FALSE]
				o.prime <- burnt.in[, o.prime.cols, drop=FALSE]
				delta.o <- o - o.prime

				on <- c + log(o.prime) - log(o)
				off <- c + log(1 - o.prime) - log(1 - o)

				site.names <- str_replace(names(means)[o.cols], "_480", "")
				site.names.on <- paste(site.names, "+", sep = "")
				site.names.off <- paste(site.names, "-", sep = "")

				p <- colMeans(on)
				pbar <- colMeans(off)
				p.sd <- colSds(on)
				pbar.sd <- colSds(off)
				final.o <- colMeans(o)
				final.o.prime <- colMeans(o.prime)
				final.delta.o <- colMeans(delta.o)
				final.o.sd <- colSds(o)
				final.o.prime.sd <- colSds(o.prime)
				final.delta.o.sd <- colSds(delta.o)
				names(p) <- site.names.on
				names(pbar) <- site.names.off
				names(p.sd) <- paste0(site.names.on, ".sd")
				names(pbar.sd) <- paste0(site.names.off, ".sd")
				names(final.o) <- paste0(site.names, ".o_480")
				names(final.o.prime) <- paste0(site.names, ".o_620")
				names(final.o.sd) <- paste0(site.names, ".o_480_sd")
				names(final.o.prime.sd) <- paste0(site.names, ".o_620_sd")

				names(final.delta.o) <- paste0(site.names, ".do")
				names(final.delta.o.sd) <- paste0(site.names, ".do_sd")
				#p.list[[ENSP]] <- colMeans(on)
				#pbar.list[[ENSP]] <- colMeans(off)
				#p.sd.list[[ENSP]] <- colSds(on)
				#pbar.sd.list[[ENSP]] <- colSds(off)
				#names(p.list[[ENSP]]) <- site.names.on
				#names(p.list[[ENSP]]) <- site.names.off
				#names(p.sd.list[[ENSP]]) <- site.names.on
				#names(p.sd.list[[ENSP]]) <- site.names.off
				return(c(c = c, c.sd = c.sd, p, pbar, p.sd, pbar.sd, final.o, final.o.prime, final.o.sd, final.o.prime.sd, final.delta.o, final.delta.o.sd))
			}
			else {
				return(c(c = c, c.sd = c.sd))
			}
		}
		else {
			warning("Error reading " %+% ENSP %+% ", skipping")
			return(c(NA, NA))
		}
	})

	c <- sapply(all.data, function(x) {x[1]}, simplify=TRUE)
	c.sd <- sapply(all.data, function(x) {x[2]})
	phosphos <- unlist(sapply(all.data, function(x) {x[-c(1, 2)]}))

	p <- phosphos[str_detect(names(phosphos), "\\+$")]
	pbar <- phosphos[str_detect(names(phosphos), "-$")]
	p.sd <- phosphos[str_detect(names(phosphos), "\\+\\.sd$")]
	pbar.sd <- phosphos[str_detect(names(phosphos), "-\\.sd$")]
	o <- phosphos[str_detect(names(phosphos), "\\.o_480$")]
	o.prime <- phosphos[str_detect(names(phosphos), "\\.o_620$")]
	o.sd <- phosphos[str_detect(names(phosphos), "\\.o_480_sd$")]
	o.prime.sd <- phosphos[str_detect(names(phosphos), "\\.o_620_sd$")]
	final.delta.o <- phosphos[str_detect(names(phosphos), ".do$")]
	final.delta.o.sd <- phosphos[str_detect(names(phosphos), ".do_sd$")]

	# Remove useless stuff from names
	names(c) <- str_replace(names(c), "\\.c$", "")
	names(c.sd) <- str_replace(names(c.sd), "\\.c.sd$", "")
	names(p.sd) <- str_replace(names(p.sd), "\\.sd$", "")
	names(pbar.sd) <- str_replace(names(pbar.sd), "\\.sd$", "")
	names(final.delta.o) <- str_replace(names(final.delta.o), "\\.do$", "")
	names(final.delta.o.sd) <- str_replace(names(final.delta.o.sd), "\\.do_sd$", "")

	names(pbar) <- str_replace(names(pbar), "_", "")
	names(pbar) <- str_replace(names(pbar), "\\.", "_")
	names(p.sd) <- str_replace(names(p.sd), "_", "")
	names(p.sd) <- str_replace(names(p.sd), "\\.", "_")

	return(list(
			concentration = list(mean = c, sd = c.sd),
			phospho = list(mean = c(p, pbar),
						   prec = 1 / c(p.sd, pbar.sd)),
			occupancies = list(mean = c(o, o.prime),
							   sd = c(o.sd, o.prime.sd)),
			delta.occupancies = list(mean = final.delta.o,
									 sd = final.delta.o.sd)
		))
}
# opt <- optim(coo.initial, fn = fit.concentration.occupancy.objective,
# 			 gr = fit.concentration.occupancy.gr,
# 			 control = list(fnscale = -1, maxit = 50000),
# 			 ratios = mu, gammas = gammas, nu.tildes = nu.tildes,
# 			 Iis = Iis, tis = tis, on.observed = on.observed, off.observed = off.observed,
# 			 npeptides = npeptides, nsites = nsites#,
# )
# View the results
# print("Result:")
# print(opt$par[1])
# print((1 + tanh(opt$par[2:4])) / 2)
# print((1 + tanh(opt$par[5:7])) / 2)
# print(opt$value)
#
# # Likelihood of the optimized parameters
# fit.concentration.occupancy.objective(opt$par, ratios = mu, gammas = gammas, nu.tildes = nu.tildes,
# 									  Iis = Iis, tis = tis, on.observed = on.observed, off.observed = off.observed,
# 									  npeptides = npeptides, nsites = nsites)
#
# # Likelihood for the true parameters:
# fit.concentration.occupancy.objective(c(c, o.to.x(o), o.to.x(o.prime)), ratios = mu, gammas = gammas, nu.tildes = nu.tildes,
# 									  Iis = Iis, tis = tis, on.observed = on.observed, off.observed = off.observed,
# 									  npeptides = npeptides, nsites = nsites)
#
#
# mu <- mu.coo(c,
# 			 matrix((1 + tanh(opt$par[2:4])) / 2, 7, 3, byrow = TRUE),
# 			 matrix((1 + tanh(opt$par[5:7])) / 2, 7, 3, byrow = TRUE), Iis, tis, npeptides, nsites)
#

# par(mfcol=c(3, 3))
# par(mar=c(2, 2, 0, 0))
# plot(mc.result[,2], type="l", ylim = 0:1)
# abline(h = o[1], col="red", lty=2)
# plot(mc.result[,3], type="l", ylim = 0:1)
# abline(h = o[2], col="red", lty=2)
# plot(mc.result[,4], type="l", ylim = 0:1)
# abline(h = o[3], col="red", lty=2)
# plot(mc.result[,5], type="l", ylim = 0:1)
# abline(h = o.prime[1], col="red", lty=2)
# plot(mc.result[,6], type="l", ylim = 0:1)
# abline(h = o.prime[2], col="red", lty=2)
# plot(mc.result[,7], type="l", ylim = 0:1)
# abline(h = o.prime[3], col="red", lty=2)
# #plot(mc.result[,8], type="l")
# plot(mc.result[,1], type="l")
# abline(h = c, col="red", lty=2)
# plot(mc.result[,8], type="l")
# plot(-mc.result[,8], type="l", log="y", ylim = rev(range(-mc.result[,8])))
#
# par(mfcol=c(3, 3))
# hist(tail(mc.result[,1], n = mc.steps / 2))
# hist(tail(mc.result[,2], n = mc.steps / 2))
# hist(tail(mc.result[,3], n = mc.steps / 2))
# hist(tail(mc.result[,4], n = mc.steps / 2))
# hist(tail(mc.result[,5], n = mc.steps / 2))
# hist(tail(mc.result[,6], n = mc.steps / 2))
# hist(tail(mc.result[,7], n = mc.steps / 2))
# hist(tail(mc.result[,8], n = mc.steps / 2))
# print(previous.coo)
