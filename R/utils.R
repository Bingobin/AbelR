# Internal helpers shared across AbelR analysis modules belong in this file.

.abel_normalize_species <- function(species) {
  species <- tolower(species[1])
  species <- switch(
    species,
    hsa = "human",
    "homo sapiens" = "human",
    mmu = "mouse",
    "mus musculus" = "mouse",
    species
  )
  match.arg(species, c("human", "mouse"))
}

.abel_gene_annotation <- function(species, gene_anno_file = NULL) {
  annotation_files <- c(
    human = "gene_len.v43.new.txt.gz",
    mouse = "gene_len.vM38.txt.gz"
  )
  option_name <- paste0("AbelR.", species, "_gene_anno_file")

  if (is.null(gene_anno_file)) {
    gene_anno_file <- getOption(option_name)
  }
  if (is.null(gene_anno_file)) {
    gene_anno_file <- system.file(
      "extdata",
      annotation_files[[species]],
      package = "AbelR"
    )
  }
  if (!file.exists(gene_anno_file)) {
    stop(
      "A valid gene annotation file is required. Supply gene_anno_file or ",
      "set options(", option_name, " = '/path/to/annotation.txt.gz')."
    )
  }

  annotation_source <- if (grepl("[.]gz$", gene_anno_file, ignore.case = TRUE)) {
    gzfile(gene_anno_file, open = "rt")
  } else {
    gene_anno_file
  }
  if (inherits(annotation_source, "connection")) {
    on.exit(close(annotation_source), add = TRUE)
  }

  annotation <- utils::read.table(
    annotation_source,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = 1
  )
  if (species == "human") {
    annotation <- annotation[
      !grepl("_PAR_", rownames(annotation)),
      ,
      drop = FALSE
    ]
  }
  rownames(annotation) <- sub("\\..*$", "", rownames(annotation))
  annotation
}
