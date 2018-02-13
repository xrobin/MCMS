library(MCMS)

context("The overall workflow works")

test_that("read.labelfree reads some data", {
	labelfree <<- read.labelfree(ups.dir, reference.experiment = "2500amol", plot = TRUE)

	expect_equal(dim(labelfree), c(47643L, 11L))
	expect_equal(names(labelfree), c("protein", "sequence", "modifications", "sample", "reference", "start", "end", "length", "ratio", "q", "n"))
	# Check one line
	one <- labelfree %>% filter(protein == "P02768ups|ALBU_HUMAN_UPS", sample == "12500amol", sequence == "LCTVATLR")

	expect_equivalent(as.data.frame(one), structure(list(protein = "P02768ups|ALBU_HUMAN_UPS", sequence = "LCTVATLR",
						  modifications = "", sample = "12500amol", reference = "2500amol",
						  start = 73L, end = 80L, length = 8L, ratio = 1.58730090501295,
						  q = 0.00994451379532787, n = 3L), class = c("data.frame")))
})


test_that("variance.model runs", {
	var.mdl <<- variance.model(labelfree)
	expect_equal(class(var.mdl), "list")

	expect_equivalent(var.mdl$mdl$rate[,"Estimate"], 0.02472515, tolerance = .001)
	expect_equivalent(var.mdl$mdl$shape[,"Estimate"], 0.8705470249, tolerance = .001)
})

test_that("MCMS can sample one protein", {
	library(dplyr)
	albumin <- Peptides(Protein(labelfree %>% filter(protein == "P02768ups|ALBU_HUMAN_UPS")))
	albumin.results <- MCMS(albumin, var.mdl, n = calcIterations(albumin) / 10)

	expect_equal(dim(albumin.results), c(7000L, 10L))
	# Test column names
	expect_true(all(c("c.12500amol_2500amol", "c.125amol_2500amol", "c.25000amol_2500amol", "c.250amol_2500amol", "c.50000amol_2500amol", "c.5000amol_2500amol", "c.500amol_2500amol", "c.50amol_2500amol", "Likelihood",     "Prior") %in% colnames(albumin.results)))
	# Some basic tests that should always hold true...
	mean.50000 <- mean(albumin.results[,"c.50000amol_2500amol"])
	mean.25000 <- mean(albumin.results[,"c.25000amol_2500amol"])
	mean.12500 <- mean(albumin.results[,"c.12500amol_2500amol"])
	mean.5000 <- mean(albumin.results[,"c.5000amol_2500amol"])
	expect_true(mean.50000 > mean.25000)
	expect_true(mean.25000 > mean.12500)
	expect_true(mean.12500 > 1)
	expect_true(mean.5000 < 1)

})


