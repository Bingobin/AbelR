#!/usr/bin/env Rscript

print_help <- function() {
  cat(
    "Usage:\n",
    "  Rscript run_cytotrace2_by_sample.R \\\n",
    "    --input_seurat input.rds \\\n",
    "    --out_rds cytotrace2_metadata.rds [options]\n\n",
    "Required:\n",
    "  --input_seurat       Input Seurat RDS file\n",
    "  --out_rds            Output combined CytoTRACE2 metadata RDS file\n\n",
    "Options:\n",
    "  --sample_col         Sample metadata column (default: SampleID)\n",
    "  --assay              Expression assay (default: RNA)\n",
    "  --slot_type          CytoTRACE2 slot/layer (default: counts)\n",
    "  --species            human or mouse (default: human)\n",
    "  --ncores             Cores used by CytoTRACE2 (default: 4)\n",
    "  --subset_col         Optional metadata column used for subsetting\n",
    "  --subset_values      Comma-separated values retained from subset_col\n",
    "  --save_result_list   Save complete per-sample results (default: FALSE)\n",
    "  --out_result_list    Output path for complete per-sample results\n",
    "  --verbose            Show analysis progress (default: TRUE)\n",
    sep = ""
  )
}

parse_args <- function(args) {
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    print_help()
    quit(save = "no", status = 0)
  }
  if (length(args) %% 2 != 0) {
    stop("Every --argument must be followed by a value.")
  }
  keys <- args[seq.int(1, length(args), by = 2)]
  values <- args[seq.int(2, length(args), by = 2)]
  if (any(!grepl("^--", keys))) {
    stop("Argument names must begin with '--'.")
  }
  stats::setNames(as.list(values), sub("^--", "", keys))
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (!is.null(value)) {
    return(value)
  }
  if (required) {
    stop("Missing required argument: --", name)
  }
  default
}

as_flag <- function(value) {
  toupper(value) %in% c("TRUE", "T", "YES", "Y", "1")
}

ensure_parent <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE)
  }
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  input_seurat <- get_arg(args, "input_seurat", required = TRUE)
  out_rds <- get_arg(args, "out_rds", required = TRUE)
  save_result_list <- as_flag(
    get_arg(args, "save_result_list", default = "FALSE")
  )
  out_result_list <- get_arg(
    args,
    "out_result_list",
    default = sub("\\.rds$", ".CytoTRACE2_result_list.rds", out_rds)
  )
  subset_values <- get_arg(args, "subset_values")
  if (!is.null(subset_values)) {
    subset_values <- trimws(strsplit(subset_values, ",", fixed = TRUE)[[1]])
  }

  if (!requireNamespace("AbelR", quietly = TRUE)) {
    stop("Install AbelR before running this script.")
  }
  message("Reading Seurat object: ", input_seurat)
  seu <- readRDS(input_seurat)
  result <- AbelR::run_cytotrace2_by_sample(
    seu = seu,
    sample_col = get_arg(args, "sample_col", "SampleID"),
    assay = get_arg(args, "assay", "RNA"),
    slot_type = get_arg(args, "slot_type", "counts"),
    species = get_arg(args, "species", "human"),
    ncores = as.integer(get_arg(args, "ncores", "4")),
    subset_col = get_arg(args, "subset_col"),
    subset_values = subset_values,
    keep_result_list = save_result_list,
    verbose = as_flag(get_arg(args, "verbose", "TRUE"))
  )

  ensure_parent(out_rds)
  saveRDS(result$metadata, out_rds)
  message("Saved CytoTRACE2 metadata: ", out_rds)
  if (save_result_list) {
    ensure_parent(out_result_list)
    saveRDS(result$results, out_result_list)
    message("Saved complete CytoTRACE2 results: ", out_result_list)
  }
}

main()
