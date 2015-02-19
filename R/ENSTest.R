#' A dummy sample data to test the stuff
#' @docType data
#' @name ENSTest
#' @examples
#' data(ENSTest)
#' \dontrun{
#' # Refresh from the source file
#' ENSTest <- read.csv(system.file("data/ENSTest.csv", package="MCMS", mustWork = TRUE))
#' ENSTest$reference <- as.character(ENSTest$reference)
#' save(ENSTest, file="data/ENSTest.RData")
#' }
#' 
NULL