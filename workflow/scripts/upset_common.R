# Shared by R_Upset_plot_{peaks,TSS,clusters}.R.

# Parse the "label=type,label=type" map from define_samples_for_upset into a
# named vector (names = labels). Split at the last '=': only the label side can
# carry one, since it is built from Levels.
parse_label_types <- function(arg) {
	pairs <- unlist(strsplit(arg, ","))
	setNames(sub("^.*=", "", pairs), sub("=[^=]*$", "", pairs))
}

# Sample columns of each type, looked up rather than pattern-matched: one type
# can be a substring of another (H3K9me / H3K9me2), and the plotting matrix
# also carries annotation columns. Ordered as in sampleslist.
upset_type_cols <- function(types, label_types, sampleslist) {
	cols <- lapply(types, function(t) intersect(sampleslist, names(label_types)[label_types == t]))
	setNames(cols, types)
}
