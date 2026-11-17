# Shared by R_Upset_plot_{peaks,TSS,clusters}.R.

# Parse a "name=value,name=value" map from define_samples_for_upset into a
# named vector (names = the left side): label=type, and sample id=label. Split
# at the last '=': only the left side can carry one, since it is built from
# Levels.
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

# The merged peak/cluster/TSS file names samples by sample id, so the merge never
# reruns over a label change; swap the ids in its comma-separated Samples column
# for the plot labels. Clusters carry the size class as '{id}_{class}', so with
# `suffixes` (the size classes) such a token maps to '{label}_{class}'. Any id
# not in the map is an error: it means the merged file and the sample sheet
# disagree, and passing it through would plot a stale or unlabelled set.
recode_sample_ids <- function(samples, id_labels, suffixes = NULL) {
	recode <- function(tok) {
		if (tok %in% names(id_labels)) return(id_labels[[tok]])
		for (sfx in suffixes) {
			id <- sub(paste0("_", sfx, "$"), "", tok)
			if (id != tok && id %in% names(id_labels)) return(paste0(id_labels[[id]], "_", sfx))
		}
		NA_character_
	}
	tokens <- strsplit(as.character(samples), ",", fixed = TRUE)
	recoded <- lapply(tokens, function(t) vapply(t, recode, ""))
	unknown <- unique(unlist(mapply(function(t, r) t[is.na(r)], tokens, recoded)))
	if (length(unknown) > 0) {
		stop("Sample id(s) in the merged file with no plot label: ",
		     paste(unknown, collapse = ", "),
		     ". The merged file predates the current sample sheet; delete it to rebuild.")
	}
	vapply(recoded, paste, "", collapse = ",")
}
