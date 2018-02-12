#' Read a label-free MaxQuant project
#' This function reads, normalises and calculates ratios from a MaxQuant project ("txt" result folder)
#' @param dir the \dQuote{txt} folder containing the MaxQuant results
#' @param reference.experiment the name of the reference experiment
#' @param score.threshold minimal identification score. Peptides with values below this threshold will be filtered out.
#' @param mod.threshold minimal site localisation probability accepted. Peptides with values below this threshold will be filtered out.
#' @param raw.file.filter an optional filter for the "Raw file" column passed to \code{\link{str_detect}}
#' @param plot show diagnostic plots
#' @import dplyr ggplot2
#' @importFrom stringr str_replace str_detect str_match
#' @importFrom xavamess safe.mapping
#' @export
read.labelfree <- function(dir, reference.experiment,
						   score.threshold = 40, mod.threshold = .9,
						   raw.file.filter = NULL,
						   plot = FALSE) {
	evidence <- read.evidence(dir) %>%
		annotate_Experiment(reference.experiment = reference.experiment)

	# Make sure we have some reference
	if (sum(evidence$Experiment == reference.experiment) == 0) {
		print(table(evidence$Experiment))
		stop("Reference experiment not found. Please select an experiment from the table above.")
	}

	if (plot) {
		ggplot(evidence, aes(Score)) + geom_histogram(aes(fill = Reverse), binwidth=0.01) + facet_grid(Experiment~.) + geom_vline(xintercept = score.threshold, linetype="dashed", alpha = .5) + scale_x_log10()  + ggtitle("Scores per experiment")
		ggplot(evidence, aes(Score)) + geom_histogram(aes(fill = Reverse)) + facet_grid(Type~run) + geom_vline(xintercept = score.threshold, linetype="dashed", alpha = .5) + scale_x_log10() + ggtitle("Scores per run and type")
	}

	msms <- read.msms(dir)

	sum.intensities.per.peptide.per.raw.file <- evidence %>%
		filter.reverse.contaminants() %>%
		filter.dirty.peaks(msms, mod.threshold = mod.threshold) %>%
		filter(Score >= score.threshold) %>%
		aggregate.raw.file.intensities


	# Optional filter
	if (!is.null(raw.file.filter)) {
		sum.intensities.per.peptide.per.raw.file <- sum.intensities.per.peptide.per.raw.file %>% filter(str_detect(Raw.file, raw.file.filter))
	}

	if (plot) {
		ggplot(sum.intensities.per.peptide.per.raw.file, aes(peptide.Intensity)) + geom_histogram(aes(fill = Raw.file), binwidth=0.1) + scale_x_log10() + facet_grid(Experiment~run) + ggtitle("Un-normalized peptide Intensities")
	}

	normalized.sum.intensities.per.peptide.per.raw.file <- sum.intensities.per.peptide.per.raw.file %>%
		normalize.labelfree()

	# Make sure intensities are > 0
	dev.null <- lapply("norm.I", function(x) {
		if (any(normalized.sum.intensities.per.peptide.per.raw.file[[x]] < 0, na.rm = TRUE)) {
			stop(sprintf("Negative %s", x))
		}
	})


	if (plot) {
		ggplot(normalized.sum.intensities.per.peptide.per.raw.file, aes(norm.I)) + geom_histogram(aes(fill = Raw.file), binwidth=0.1) + scale_x_log10() + facet_grid(Experiment~run) + ggtitle("Normalized peptide Intensities")
	}

	summarized.normalized.mq.data <- normalized.sum.intensities.per.peptide.per.raw.file %>%
		aggregate.replica.intensities() %>% ungroup %>% print

	summarized.normalized.control <- summarized.normalized.mq.data %>%
		filter(Experiment == reference.experiment) %>%
		rename(
			n.evidence.control = n.evidence,
			n.replicate.control = n.replicate,
			Retention.length.control = Retention.length,
			Number.of.data.points.control = Number.of.data.points,
			I.tot.control = I.tot,
			I.mean.control = I.mean,
			I.sd.control = I.sd,
			norm.I.mean.control = norm.I.mean,
			norm.I.sd.control = norm.I.sd
		) %>%
		select(-Experiment) %>% print
	ratios <- left_join(summarized.normalized.mq.data, summarized.normalized.control, by = c("Modified.sequence", "run")) %>%
		mutate(
			# log ratio
			ratio = log(I.mean.control / I.mean),
			norm.ratio = log(norm.I.mean / norm.I.mean.control),
			# sd
			ratio.sd = sqrt((I.sd / I.mean)^2 + (I.sd.control / I.mean.control)^2),
			norm.ratio.sd = sqrt((norm.I.sd / norm.I.mean)^2 + (norm.I.sd.control / norm.I.mean.control)^2),
			# n
			n.replicate.total = n.replicate.control + n.replicate,
			n.eff = pmin(n.replicate.control, n.replicate),
			# q
			ratio.q = (n.replicate.total-1) * norm.ratio.sd^2,
			norm.ratio.q = (n.replicate.total-1) * norm.ratio.sd^2,

			# Fix peptide.ID
			Peptide.ID = ifelse(is.na(Peptide.ID.x), Peptide.ID.y, Peptide.ID.x)
		) %>%
		select(-Peptide.ID.x, -Peptide.ID.y) %>% print


	# Map peptides to the original protein
	peptides <- read.peptides(dir)

	n.eff <- ratios %>% safe.mapping(peptides, by = "Peptide.ID") %>%
		filter(!is.na(n.eff), # remove if n is missing
			   !is.na(norm.ratio),
			   Experiment != control.experiment) %>%
		mutate(
			reference = control.experiment,
			# Calculate modifications
			modifications = constructModifiedPeptide(Modified.sequence, Start.position)#,
		) %>%
		select(Leading.razor.protein, Sequence, modifications, Experiment, reference,
			   Start.position, End.position, Length,
			   norm.ratio, norm.ratio.q, n.eff) %>%
		rename(
			protein = Leading.razor.protein,
			sequence = Sequence,
			start = Start.position,
			end = End.position,
			length = Length,
			sample = Experiment,
			ratio = norm.ratio,
			q = norm.ratio.q,
			n = n.eff
		)

	# Replace q with 0 when n == 1 (instead of NA)
	n.eff$q[n.eff$n == 1] <- 0

	# Ensure sequence matches length...
	stopifnot(identical(n.eff$length, nchar(n.eff$sequence)))
	stopifnot(identical(n.eff$length + n.eff$start - 1L, n.eff$end))

	if (plot) {
		n.eff.ggplot <- n.eff[sample(seq_len(nrow(n.eff))),]
		g <- ggplot(n.eff.ggplot, aes(ratio, q)) + geom_point(aes(color = sample, size = n)) + ggtitle("Effective observations")
		print(g)
	}

	return(n.eff)
}


#' Normalize the sum intensities data and a linear model
#' @param data the MS data to normalize
#' @param lm.model an optional \code{\link{lm}} model of the form lm(log.Intensity ~ Raw.file). If provided, a column named norm.I will be calculated. Alternative models can be used if they provide compatible methods for \code{\link{coef}} and \code{\link{predict}}
#' @export
normalize.labelfree <- function(data, lm.model) {
	sum.intensities <- data %>%
		# Integrate all Modified.sequences
		group_by(Raw.file, replicate, run, Experiment) %>%
		summarize(I.tot = sum(peptide.Intensity, na.rm = TRUE),
				  I.tot.log = sum(log(peptide.Intensity), na.rm = TRUE)
				  )

	# Apply the normalization
	data <- data %>%
		left_join(sum.intensities %>% ungroup %>% select(-replicate, -run, -Experiment), by = "Raw.file")

	# Check if we got an lm model
	if (missing(lm.model) || is.null(lm.model)) {
		lm.model <- lm(log.peptide.Intensity ~ Raw.file, data = data)
	}

	normalization.coefs <- coef(lm.model)
	names(normalization.coefs) <- str_replace(names(normalization.coefs), "^Raw.file", "")
	# Get the lm normalization column name
	norm.colname <- all.vars(terms(lm.model))[attr(terms(lm.model), "response")]
	# Make sure we have all the Raw.files
	if (any(which <- ! data$Raw.file %in% names(normalization.coefs))) {
		missing <- unique(data$Raw.file[which])
		# Maybe we can fix that?
		# Let's make sure we have only 1 missing, and that if we predict we get the intercept term...
		if (length(missing == 1) && predict(lm.model, newdata = data.frame(Raw.file = missing)) == normalization.coefs["(Intercept)"]) {
			# Then don't even warn anymore, it is safe now
			# warning(sprintf("Missing %s Raw file in coefficients: %s. Using the intercept term", length(missing), paste(missing, collapse=", ")))
			normalization.coefs[missing] <- 0
		}
		else {
			stop(sprintf("Missing %s Raw file in coefficients: %s", length(missing), paste(missing, collapse=", ")))
		}

	}
	data$norm.I = exp(data[[norm.colname]] - normalization.coefs[data$Raw.file])
	return(data)
}
