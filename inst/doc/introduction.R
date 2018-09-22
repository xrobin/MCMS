## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE--------------------------------------------------------
#  library(devtools)
#  devtools::install_github("xrobin/MCMS")

## ------------------------------------------------------------------------
library(MCMS)

## ------------------------------------------------------------------------
dl.dst <- file.path(tempdir(), "UPS.tar.bz2")
download.file("http://lindinglab.org/downloads/UPS.tar.bz2", dl.dst)
ups.dir <- file.path(tempdir(), "UPS")
untar(dl.dst, exdir=ups.dir, compressed = "bzip2")

## ------------------------------------------------------------------------
labelfree <- read.labelfree(ups.dir, reference.experiment = "2500amol", plot = TRUE)

## ------------------------------------------------------------------------
var.model <- variance.model(labelfree)

## ------------------------------------------------------------------------
library(dplyr)
albumin <- Peptides(Protein(labelfree %>% filter(protein == "P02768ups|ALBU_HUMAN_UPS")))
albumin.results <- MCMS(albumin, var.model, n = calcIterations(albumin) / 10)
# Results is a matrix 
head(albumin.results)

## ------------------------------------------------------------------------
palette <- rainbow(8)
matplot(albumin.results[,c(8, 2, 4, 7, 6, 1, 3, 5)], type="l", col = palette)

## ------------------------------------------------------------------------
attr(albumin.results, "effectiveSize")
if (any(attr(albumin.results, "effectiveSize") < 100)) {
	warning("You should consider sampling again with an increased n.")
}

