% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeptidesClass.R
\name{Peptides}
\alias{Peptides}
\title{Creates a Peptides model with concentration and occupancies for a Protein.}
\usage{
Peptides(protein)
}
\arguments{
\item{protein}{a \code{\link{Protein-class}} object}
}
\description{
Creates a Peptides model with concentration and occupancies for a Protein.
}
\examples{
data(ENSTest)
ENSTestProtein <- Protein(ENSTest)
ENSTestModel <- Peptides(ENSTestProtein)
}
