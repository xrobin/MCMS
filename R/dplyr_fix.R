## quiets concerns of R CMD check re non-standard evaluation
# Messages like:
# make_unique_modifications: no visible binding for global variable
# 'start'
# Undefined global functions or variables:
#   modifications pair ratio start
# See https://stackoverflow.com/questions/9439256/ for possibly updated fixes:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("modifications", "pair", "ratio", "start",
".", "Acetyl..Protein.N.term.", "Calibrated.retention.time", "End.position",
"Evidence.ID", "Experiment", "I.mean", "I.mean.control", "I.sd", "I.sd.control", "I.tot",
"Intensity", "Leading.Razor.Protein", "Leading.razor.protein", "Length", "Mass",
"Match.score", "Modifications", "Modified.sequence", "Number.of.data.points",
"Oxidation..M.", "Oxidation..M..Probabilities", "PEP", "Peptide.ID", "Peptide.ID.x",
"Peptide.ID.y", "Phospho..STY.", "Phospho..STY..Probabilities",
"Potential.contaminant", "Raw.file", "Retention.length", "Reverse", "Scan.number",
"Score", "Sequence", "Start.position", "Type",
"broken", "lu", "modified.peptide.p", "n.evidence", "n.replicate", "n.replicate.control",
"n.replicate.total", "norm.I", "norm.I.mean",
"norm.I.mean.control", "norm.I.sd", "norm.I.sd.control", "norm.ratio",
"norm.ratio.q", "norm.ratio.sd", "peptide.Intensity", "reference", "run", "variance"))