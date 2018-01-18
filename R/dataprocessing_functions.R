#' Reads a MQ file with the proper format: header, tab separator, no stringsAsFactor, no quotes
#' @return the MQ file as a \cpode{\link{data.frame}}
#' @param MQ.file the file name
#' @param ... additional arguments for \code{\link{read.table}} such as \code{nrow}, etc.
#' @export
read.MQ <- function(MQ.file, ...) {
	read.table(MQ.file, header = TRUE, sep="\t", stringsAsFactors = FALSE, quote = "", nrows = nrows, ...)
}


annotate_Experiment <- function(data, reference.experiment) {
	raw.file.breakdown <- str_match(data$Raw.file, "UPS1_((\\d+)amol)_R(\\d)$")
	data$replicate <- raw.file.breakdown[,3]
	data$Experiment <- raw.file.breakdown[,2]
	data$concentration <- raw.file.breakdown[,3]

	if (any(is.na(data$replicate))) {
		stop("Couldn't assign some replicates")
	}

	# Some common annotations
	data$run <- "Proteome"
	data$reference <- data$Experiment == reference.experiment

	return(data)
}

# Read an evidence.txt file from folder
# @param MQ.dir the directory containing the MQ files to process
# @param nrows see \code{read.table}
read.evidence <- function(MQ.dir, nrows = -1) {
	evidence.file <- file.path(MQ.dir, "evidence.txt")
	evidence <- read.MQ(evidence.file, nrows = nrows) %>%
		mutate(
			Reverse = Reverse == "+",
			Potential.contaminant = Potential.contaminant == "+",
			Score = ifelse(Type == "MULTI-MATCH", Match.score, Score)
		) %>%
		select(Raw.file, id, Leading.Razor.Protein,
			   Sequence, Modified.sequence, Score, PEP, Mass, Type,
			   Modifications, Oxidation..M..Probabilities,
			   Acetyl..Protein.N.term., Oxidation..M.,
			   Intensity, Retention.length, Calibrated.retention.time,
			   Number.of.data.points, Reverse, Potential.contaminant,
			   Peptide.ID) %>%
	rename(Evidence.ID = id)

	# Check reverse and contaminents match with REV__ and CON__
	print(table(evidence$Reverse, grepl("^REV_", evidence$Leading.Razor.Protein)))
	print(table(evidence$Potential.contaminant, grepl("^CON_", evidence$Leading.Razor.Protein)))

	return(evidence)

}

# Get the peptide positions on protein from peptides.txt
# @param MQ.dir the directory containing the MQ files to process
# @param nrows see \code{read.table}
read.peptides <- function(MQ.dir, nrows = -1) {
	peptides.file <- file.path(MQ.dir, "peptides.txt")
	peptides <- read.MQ(peptides.file, nrows = nrows) %>%
		select(Sequence, id, Leading.razor.protein,
			   Length, Start.position, End.position) %>%
		rename(Peptide.ID = id)
	return(peptides)
}

# Read an msms.txt file from folder
# @param MQ.dir the directory containing the MQ files to process
# @param nrows see \code{read.table}
read.msms <- function(MQ.dir, nrows = -1) {
	msms.file <- file.path(MQ.dir, "msms.txt")
	msms <- read.MQ(msms.file, nrows = nrows) %>%
		select(Raw.file, Scan.number, id, Evidence.ID, Sequence, Modified.sequence, Score, PEP, Mass, Type,
			   Modifications, Oxidation..M..Probabilities,
			   Acetyl..Protein.N.term., Oxidation..M.) %>%
		rename(MSMS.id = id)
}


# Reads the msms.txt and evidence.txt files and returns it after some basic processing:
# - removal of spectras where there is a suspicion of co-fragmentation (same Raw.file, Scan.number, Scan.index; OR MULTI-SECPEP Type)
# - Suppression of Pyro-Glu/n modifications that are technical artefacts
# @param MQ.dir the directory containing the MQ files to process
# @param limit.nrow.evidence,limit.nrow.msms read no more than this many rows from the files. Useful for quick superficial debugging
# @param filter whether to filter reverse matches/contaminants, dirty peaks and low score peaks
# @param mod.threshold,score.threshold threshold on the modification site probability and score, if \code{filter = TRUE}.
# @import dplyr
read.data <- function(MQ.dir, nrows.evidence = -1, nrows.msms = -1, filter = FALSE, mod.threshold = 0.9, score.threshold = 40) {
	stop("Function is deprecated. Use read.evidence and read.msms instead")
	evidence.file <- file.path(MQ.dir, "evidence.txt")
	msms.file <- file.path(MQ.dir, "msms.txt")

	evidence <- read.MQ(evidence.file, nrows = nrows.evidence)
	#evidence <- evidence %>%
	#	select(Raw.file, id, Leading.Razor.Protein, Experiment,
	#		   Intensity, Retention.length, Number.of.data.points,
	#		   Reverse) %>%
	#	mutate(
	#		Reverse = Reverse == "+",
	#		Potential.contaminant = Potential.contaminant == "+"
	#	)

	msms <- read.MQ(msms.file, nrows = nrows.msms)
	save(evidence, msms, file = "mmap.mq.RData")

	# Verify multi-secpep peptides
	# Isolate MULTI-SECPEP
	#multi.secpep <- msms %>% filter(Type == "MULTI-SECPEP") %>%
	##	select(Scan.number, Raw.file) %>%
	#	left_join(msms) %>%
	#	group_by(Scan.number, Raw.file) %>%
	#	summarize(mzrange = diff(range(m.z)), n = n())
	#multi.secpep %>% filter(n > 1)
	#hist(multi.secpep$mzrange)


	# Get the corresponding primary peptide
	#scan.summary <- msms %>%
	#	group_by(Scan.number, Raw.file) %>%
	#	filter(n() > 1) %>%
	#	summarize(mzrange = diff(range(m.z)))
	#hist(scan.summary$mzrange)

	#test <- data %>% group_by(Scan.number, Raw.file) %>% summarize(n = n())
	# Here is an OK case: mass and m/z actually differ
	#msms %>% filter(Scan.number == 1963, Raw.file == "20161202_PE_FV_MM_phospho_1min_rep2")
	# Here is a bad one: same MSMS used in 2 different evidence IDs
	#msms %>% filter(Scan.number == 1577, Raw.file == "20161202_PE_FV_MM_phospho_30min_rep1")

	results <- msms %>%
		select(Raw.file, Scan.number, id, Evidence.ID, Sequence, Modified.sequence, Score, PEP, Mass, Type,
			   Modifications, Oxidation..M..Probabilities, Phospho..STY..Probabilities,
			   Acetyl..Protein.N.term., Oxidation..M., Phospho..STY.) %>%
		rename(MSMS.id = id) %>%
		#drop.pyro.Glu %>%
		# Merge msms.txt with evidence
		#stop("Don't left_join evidence. Just process it separately, list Evidence.IDs to keep or remove")
		left_join(evidence, by = c("Evidence.ID" = "id", "Raw.file")) %>% # Join by Evidence ID AND Raw File - somehow an evidence can have multiple Raw file
		mutate(file = evidence.file) %>% # Add an ID unique to this file set, otherwise Evidence.IDs are void...
		add.group.indices
		# Filtering
		if (filter) {
			results <- results %>%
				filter.reverse.contaminants() %>%
				filter.dirty.peaks(msms = broken, mod.threshold = mod.threshold) %>%
				filter(Score >= score.threshold)
		}
		return(results)
}


# Drops x$Modification == "Glu->pyro-Glu" || x$Modification == "Gln->pyro-Glu".
# Rationale: this is a purely technical artefact during the sample prep/processing. Not biologically relevant
# @param x the peptides (with column \code{Modifications}) in which to drop Pyro-Glu/n
drop.pyro.Glu <- function(x) {
	x %>% mutate(
		Modifications = str_replace(Modifications, ",?Gl[un]->pyro-Glu", ""),
		Modifications = ifelse(Modifications == "", "Unmodified", Modifications)
	)
}

# Removes contaminants and matches in the reverse database...
filter.reverse.contaminants <- function(data) {
	data %>% filter(
		! str_detect(Leading.Razor.Protein, "^CON__"),
		! str_detect(Leading.Razor.Protein, "^REV__")
	) %>%
		return
}


# Get the probability that the modified sequence was observed, from a MQ sequence probability string.
# @param probabilities.sequence the sequence with probabilties, such as "S(0.027)DDY(0.969)MPMS(0.002)PAS(0.001)VS(0.001)APK" directly from the Phospho..STY..Probabilities column
# @param phospho.Modified.sequence The Modified.sequence from MQ, with the leading/trailing underscores and other modifications removed
# @param n.on number of modifications that are on - used only to control the processing, could be removed in the future
# @return The probability of the modified sequence
# @importFrom stringr str_locate_all
# get.modification.peptide.probability <- function(probabilities.sequence, Modified.sequence, n.on) {
# 	# Get the position of the sites selected by MQ
# 	chosen.sites <- str_locate_all(Modified.sequence, "\\([a-z]+\\)")[[1]]
# 	location.differences <- cumsum(chosen.sites[, 2] - chosen.sites[, 1] + 1)
# 	chosen.sites.positions <- c(chosen.sites[1], chosen.sites[-1, 1] - location.differences[-length(location.differences)]) - 1
#
# 	# Get the positions of the quantified sites
# 	available.sites <- str_locate_all(probabilities.sequence, "\\([0-9.]+\\)")[[1]]
# 	location.differences <- cumsum(available.sites[, 2] - available.sites[, 1] + 1)
# 	available.positions <- c(available.sites[1], available.sites[-1, 1] - location.differences[-length(location.differences)]) - 1
#
# 	# Get per-site probabilities
# 	probabilities <- as.numeric(str_sub(probabilities.sequence, available.sites[,1] + 1, available.sites[,2] - 1))
# 	chosen.sites.probabilities <- probabilities[available.positions %in% chosen.sites.positions]
#
# 	if (length(chosen.sites.probabilities) != n.on) {
# 		stop(sprintf("Wrong number of sites detected: %s, %s, expected %s, got %s", probabilities.sequence, Modified.sequence, n.on, length(chosen.sites.probabilities)))
# 	}
#
# 	p <- prod(chosen.sites.probabilities)
# 	if (length(p) != 1) {
# 		browser()
# 	}
# 	return(p)
# }

# Removes evidence peaks that are matched with multiple different identifications in the MSMS
# Specifically, it looks for a unique Modified.Sequence across the peak, and high confidence of phospho/oxidation localization
# @param data an msms.txt data.frame
# @param mod.threshold the probability threshold for phospho/oxidation localization. Peptides with probability lower than this will be removed
filter.dirty.peaks <- function(evidence, msms, mod.threshold) {
	# Find evidence peaks with more than 1 ID
	dirty.peaks <- msms %>%
		group_by(Raw.file,Evidence.ID) %>%
		summarize(lu = length(unique(Modified.sequence))) %>%
		ungroup() %>%
		filter(lu > 1) %>%
		select(Evidence.ID)

		#filter(Phospho..STY. == 2) %>% head %>%
		#select(Modified.sequence, Phospho..STY..Probabilities, Phospho..STY.) %>%

	# Some evidence peptides have Phospho..STY. > 0 yet Phospho..STY..Probabilities == "". Remove them
	if ("Phospho..STY." %in% names(evidence)) {
		weird.missing.probs <- (evidence$Phospho..STY. > 0 & evidence$Phospho..STY..Probabilities == "" & evidence$Type != "MULTI-MATCH")
		if (any(weird.missing.probs)) {
			warning(sprintf("Filtering %s evidence IDs with Phospho but no probability: %s",
							sum(weird.missing.probs),
							paste(evidence$Evidence.ID[weird.missing.probs], collapse = ", "))
			)
			evidence <- evidence %>% filter(! weird.missing.probs)
		}
	}
	evidence %>%
		mutate(
			# With phospho:
			# modified.peptide.p = ifelse(Type == "MULTI-MATCH", 1, calcModifiedPeptideP(Modified.sequence, Phospho..STY..Probabilities, Oxidation..M..Probabilities))
			# Without phospho:
			modified.peptide.p = ifelse(Type == "MULTI-MATCH", 1, calcModifiedPeptideP(Modified.sequence, NULL, Oxidation..M..Probabilities))
		) %>%
		ungroup %>%
		filter(
			modified.peptide.p >= mod.threshold,
			! Evidence.ID %in% dirty.peaks$Evidence.ID) %>%
		return
}



# Adds indices for the groups (Evidence, MSMS and Scan number)
add.group.indices <- function(data) {
	data %>% mutate(
		evidence = group_indices_(data, .dots=c("Evidence.ID", "file", "Raw.file")),
		msms = group_indices_(data, .dots=c("MSMS.id", "file", "Raw.file")),
		scan = group_indices_(data, .dots=c("Scan.number", "file", "Raw.file"))
	)
}

# Adds pair indices for evidence and msms IDs, named evidence.pair and msms.pair.
# After create.pairs.from.triplets, we'll have several source.masses for a given evidence or msms ID and must keep track of those
# Otherwise source.masses is not present, and evidence.pair and msms.pair will simply correspond (but not equal) to the evidence and msms IDs.
add.pair.indices <- function(pairs) {
	if ("source.masses" %in% names(pairs)) {
		pairs %>% mutate(
			evidence.pair = group_indices_(pairs, .dots=c("evidence", "source.masses")),
			msms.pair = group_indices_(pairs, .dots=c("msms", "source.masses"))
		) %>% return
	}
	else {
		pairs %>% mutate(
			evidence.pair = group_indices_(pairs, .dots=c("evidence")),
			msms.pair = group_indices_(pairs, .dots=c("msms"))
		) %>% return
	}
}


# Aggregate peptide over a Raw file
# @param x data.frame out of read.evidence and normalized
aggregate.raw.file.intensities <- function(x) {
	x %>% group_by(Modified.sequence, Peptide.ID, Raw.file, replicate, run, Experiment) %>%
		summarize(
			peptide.Intensity = sum(Intensity),
			log.peptide.Intensity = log(peptide.Intensity),
			n.evidence = n(),
			Retention.length = sum(Retention.length),
			Number.of.data.points = sum(Number.of.data.points)
		)
}

# Aggregate peptide over the replicates
# Calculates sums and means
# @param x data.frame out of aggregate.raw.file.intensities
aggregate.replica.intensities <- function(x) {
	x %>% group_by(Modified.sequence, Peptide.ID, run, Experiment) %>%
		summarize(
			# General statistics
			n.evidence = sum(n.evidence),
			n.replicate = n(),
			Retention.length = sum(Retention.length),
			Number.of.data.points = sum(Number.of.data.points),
			# Unnormalized
			I.tot = sum(peptide.Intensity),
			I.mean = mean(peptide.Intensity),
			I.sd = sd(peptide.Intensity),

			# lm normalized
			norm.I.mean = mean(norm.I),
			norm.I.sd = sd(norm.I)
		)
}
