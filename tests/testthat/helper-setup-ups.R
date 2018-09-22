# Prepare the tests by downloading some data
dl.dst <- file.path(tempdir(), "UPS.tar.bz2")
ups.dir <- file.path(tempdir(), "UPS")
if (! all(file.exists(file.path(ups.dir, c("evidence.txt", "msms.txt", "peptides.txt"))))) {
	if (! file.exists(dl.dst)) {
		download.file("https://www.lindinglab.org/downloads/UPS.tar.bz2", dl.dst)
	}
	untar(dl.dst, exdir=ups.dir, compressed = "bzip2")
}
