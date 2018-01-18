#' Read a label-free MaxQuant project
#' This function reads, normalises and calculates ratios from a MaxQuant project ("txt" result folder)
#' @param dir the \dQuote{txt} folder containing the MaxQuant results
#' @param reference.experiment the name of the reference experiment
#' @param score.threshold minimal identification score. Peptides with values below this threshold will be filtered out.
#' @param mod.threshold minimal site localisation probability accepted. Peptides with values below this threshold will be filtered out.
#' @param raw.file.filter an optional filter for the "Raw file" column passed to \code{\link{str_detect}}
#' @import dplyr
#' @importFrom stringr str_replace str_detect
#' @importFrom xavamess safe.mapping
#' @export
read.labelfree <- function(dir, reference.experiment,
						   score.threshold = 40, mod.threshold = .9) {
	evidence <- read.evidence(dir) %>%
		annotate_Experiment(reference.experiment = reference.experiment)

	msms <- read.msms(dir)

	sum.intensities.per.peptide.per.raw.file <- evidence %>%
		filter.reverse.contaminants() %>%
		filter.dirty.peaks(msms, mod.threshold = mod.threshold) %>%
		filter(Score >= score.threshold) %>%
		aggregate.raw.file.intensities

	# Optional filter
	if (!is.null(raw.file.filter)) {
		sum.intensities.per.peptide.per.raw.file <- sum.intensities.per.peptide.per.raw.file %>% filter(str_detect(Raw.file, "UPS1_\\d{4,5}amol"))
	}

	normalized.sum.intensities.per.peptide.per.raw.file <- sum.intensities.per.peptide.per.raw.file %>%
		normalize()

	# Make sure intensities are > 0
	#dev.null <- lapply(c("normmedian.I", "normgmean.I", "normlm.I"), function(x) {
	dev.null <- lapply("normlm.I", function(x) {
		if (any(normalized.sum.intensities.per.peptide.per.raw.file[[x]] < 0, na.rm = TRUE)) {
			stop(sprintf("Negative %s", x))
		}
	})

	#### Choose a final normalization
	normalized.sum.intensities.per.peptide.per.raw.file$norm.I <- normalized.sum.intensities.per.peptide.per.raw.file$normlm.I

	#### Aggregate the replicates
	# summarized.normalized.spikein.mq.data <- normalized.spikein.sum.intensities %>%
	# 	aggregate.replica.intensities() %>% print

	summarized.normalized.ups.mq.data <- normalized.sum.intensities.per.peptide.per.raw.file %>%
		aggregate.replica.intensities() %>% ungroup %>% print

	summarized.normalized.ups.control <- summarized.normalized.ups.mq.data %>%
		filter(Experiment == reference.experiment) %>%
		rename(
			n.evidence.control = n.evidence,
			n.replicate.control = n.replicate,
			Retention.length.control = Retention.length,
			Number.of.data.points.control = Number.of.data.points,
			I.tot.control = I.tot,
			I.mean.control = I.mean,
			I.sd.control = I.sd,
			#normmedian.I.mean.control = normmedian.I.mean,
			#normmedian.I.sd.control = normmedian.I.sd,
			#normgmean.I.mean.control = normgmean.I.mean,
			#normgmean.I.sd.control = normgmean.I.sd,
			normlm.I.mean.control = normlm.I.mean,
			normlm.I.sd.control = normlm.I.sd
		) %>%
		select(-Experiment) %>% print
	ups.ratios <- left_join(summarized.normalized.ups.mq.data, summarized.normalized.ups.control, by = c("Modified.sequence", "run")) %>%
		mutate(
			# log ratio
			ratio = log(I.mean.control / I.mean),
			#normmedian.ratio = log(normmedian.I.mean / normmedian.I.mean.control),
			#normgmean.ratio = log(normgmean.I.mean / normgmean.I.mean.control),
			normlm.ratio = log(normlm.I.mean / normlm.I.mean.control),
			# sd
			ratio.sd = sqrt((I.sd / I.mean)^2 + (I.sd.control / I.mean.control)^2),
			#normmedian.ratio.sd = sqrt((normmedian.I.sd / normmedian.I.mean)^2 + (normmedian.I.sd.control / normmedian.I.mean.control)^2),
			#normgmean.ratio.sd = sqrt((normgmean.I.sd / normgmean.I.mean)^2 + (normgmean.I.sd.control / normgmean.I.mean.control)^2),
			normlm.ratio.sd = sqrt((normlm.I.sd / normlm.I.mean)^2 + (normlm.I.sd.control / normlm.I.mean.control)^2),
			# n
			n.replicate.total = n.replicate.control + n.replicate,
			n.eff = pmin(n.replicate.control, n.replicate),
			# q
			ratio.q = (n.replicate.total-1) * normlm.ratio.sd^2,
			#normmedian.ratio.q = (n.replicate.total-1) * normmedian.ratio.sd^2,
			#normgmean.ratio.q = (n.replicate.total-1) * normgmean.ratio.sd^2,
			normlm.ratio.q = (n.replicate.total-1) * normlm.ratio.sd^2,

			# Fix peptide.ID
			Peptide.ID = ifelse(is.na(Peptide.ID.x), Peptide.ID.y, Peptide.ID.x)
		) %>%
		select(-Peptide.ID.x, -Peptide.ID.y) %>% print


	# Map peptides to the original protein
	peptides <- read.peptides(dir)
	#mapped.ratios <- ups.ratios %>% safe.mapping(peptides, by = "Peptide.ID")

	#ups.ratios$Peptide <- str_replace_all(ups.ratios$Modified.sequence, "(_|\\([a-z]+\\))", "")
	#mapping <- read.mapper.file(longest.isoform.map.file)
	#if (anyDuplicated(mapping$Peptide)) {stop("Duplicated peptides in mapper file!")}
	n.eff <- ups.ratios %>% safe.mapping(peptides, by = "Peptide.ID") %>%
		filter(!is.na(n.eff), # remove if n is missing
			   !is.na(normlm.ratio),
			   Experiment != control.experiment) %>%
		mutate(
			reference = control.experiment,
			# Calculate modifications
			modifications = xavamess::constructModifiedPeptide(Modified.sequence, Start.position)#,
			#length = nchar(Peptide),
			#end = Position + length - 1
		) %>%
		select(Leading.razor.protein, Sequence, modifications, Experiment, reference,
			   Start.position, End.position, Length,
			   normlm.ratio, normlm.ratio.q, n.eff) %>%
		rename(
			protein = Leading.razor.protein,
			sequence = Sequence,
			start = Start.position,
			end = End.position,
			length = Length,
			sample = Experiment,
			ratio = normlm.ratio,
			q = normlm.ratio.q,
			n = n.eff
		)

	# Ensure sequence matches length...
	stopifnot(identical(n.eff$length, nchar(n.eff$sequence)))
	stopifnot(identical(n.eff$length + n.eff$start - 1L, n.eff$end))
}


#' Normalize the sum intensities data and a linear model
#' @param data the MS data to normalize
#' @param lm.model an optional \code{\link{lm}} model of the form lm(log.Intensity ~ Raw.file). If provided, a column named norm.I will be calculated. Alternative models can be used if they provide methods for \code{\link{coef}} and \code{\link{predict}}
#'
normalize.labelfree <- function(data, lm.model) {
	sum.intensities <- data %>%
		# Integrate all Modified.sequences
		group_by(Raw.file, replicate, run, Experiment) %>%
		summarize(I.tot = sum(peptide.Intensity, na.rm = TRUE),
				  I.tot.log = sum(log(peptide.Intensity), na.rm = TRUE)#,
				  #I.gmean = exp(mean(log(peptide.Intensity), na.rm = TRUE)),
				  #I.median = median(peptide.Intensity, na.rm = TRUE),
				  #I.gsd = exp(sd(log(peptide.Intensity), na.rm = TRUE)),
				  #I.IQR = IQR(peptide.Intensity, na.rm = TRUE)
				  )#,
	#n = n(), n.evidence = sum(n.evidence))

	# Get the grand mean
	#grand.mean <- mean(sum.intensities$I.gmean) # If the runs behave normally
	#geom.grand.mean <- exp(mean(log(sum.intensities$I.gmean))) # If the runs behave log-normally, geometric grand mean
	# Grand median
	#grand.median <- median(sum.intensities$I.median)

	# Apply the normalization
	data <- data %>%
		left_join(sum.intensities %>% ungroup %>% select(-replicate, -run, -Experiment), by = "Raw.file")# %>%
		#mutate(
		#	#norm.I = (peptide.Intensity - I.median) / I.IQR  + grand.median,
		#	normmedian.I = exp((log(peptide.Intensity) - log(I.median)) / log(I.IQR) + log(grand.median)),
		#	#norm2.I = peptide.Intensity / I.tot,
		#	#norm3.I = exp(log.peptide.Intensity - I.tot.log),
		#	#normgmean.I = exp(log.peptide.Intensity - (log(I.gmean) - log(grand.mean))),
		#	normgmean.I = exp(log.peptide.Intensity - log(I.gmean) + log(geom.grand.mean))
		#	#normlm.I = exp(log.Intensity - predict(lm.raw.file, data))
		#)

	# Check if we got an lm model
	if (is.null(lm.model)) {
		lm.model <- lm(log.peptide.Intensity ~ Raw.file, data = sum.intensities.per.peptide.per.raw.file)
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
			# Then just warn
			warning(sprintf("Missing %s Raw file in coefficients: %s. Using the intercept term", length(missing), paste(missing, collapse=", ")))
			normalization.coefs[missing] <- 0
		}
		else {
			stop(sprintf("Missing %s Raw file in coefficients: %s", length(missing), paste(missing, collapse=", ")))
		}

	}
	data$norm.I = exp(data[[norm.colname]] - normalization.coefs[data$Raw.file])
	return(data)
}
