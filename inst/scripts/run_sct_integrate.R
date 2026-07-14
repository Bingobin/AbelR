#!/usr/bin/env Rscript

print_help <- function() {
  cat(
    "Usage:\n",
    "  Rscript run_sct_integrate.R \\\n",
    "    --input_seurat input.rds \\\n",
    "    --out_seurat integrated.rds [options]\n\n",
    "Required:\n",
    "  --input_seurat        Input Seurat RDS file\n",
    "  --out_seurat          Output integrated Seurat RDS file\n\n",
    "QC options:\n",
    "  --min_nFeature_RNA    Default: 300\n",
    "  --max_nFeature_RNA    Default: 9000\n",
    "  --max_percent_mt      Default: 12.5\n",
    "  --min_nCount_RNA      Default: 800\n",
    "  --max_nCount_RNA      Default: 60000\n",
    "  --doublet_col         Default: scDblFinder.class\n",
    "  --doublet_keep        Default: singlet\n",
    "  --use_doublet_filter  Default: TRUE\n\n",
    "SCT/integration options:\n",
    "  --sample_col          Default: SampleID\n",
    "  --var_features_n      Default: 3000\n",
    "  --sct_method          Default: glmGamPoi\n",
    "  --npcs                Default: 30\n",
    "  --dims                Range or list, for example 1:30 or 1,2,3\n",
    "  --k_anchor            Default: 20\n",
    "  --k_weight            Default: 20\n",
    "  --resolution          Default: 0.5\n",
    "  --umap_n_neighbors    Default: 30\n",
    "  --umap_min_dist       Default: 0.3\n",
    "  --umap_spread         Default: 1.0\n",
    "  --join_layers         Default: TRUE\n",
    "  --maxSize_GB          Default: 35\n",
    "  --workers             future workers; 1 uses sequential (default: 1)\n",
    "  --verbose             Default: FALSE\n",
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

parse_dims <- function(value) {
  value <- gsub(" ", "", value, fixed = TRUE)
  if (grepl(":", value, fixed = TRUE)) {
    bounds <- strsplit(value, ":", fixed = TRUE)[[1]]
    if (length(bounds) != 2) {
      stop("--dims must be a range such as 1:30 or a comma-separated list.")
    }
    return(seq.int(as.integer(bounds[1]), as.integer(bounds[2])))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
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
  out_seurat <- get_arg(args, "out_seurat", required = TRUE)
  workers <- as.integer(get_arg(args, "workers", "1"))
  max_size_gb <- as.numeric(get_arg(args, "maxSize_GB", "35"))
  if (is.na(workers) || workers < 1) {
    stop("--workers must be a positive integer.")
  }

  if (!requireNamespace("AbelR", quietly = TRUE)) {
    stop("Install AbelR before running this script.")
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for this script.")
  }
  old_plan <- future::plan()
  old_max_size <- getOption("future.globals.maxSize")
  on.exit({
    future::plan(old_plan)
    options(future.globals.maxSize = old_max_size)
  }, add = TRUE)
  if (workers == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  options(future.globals.maxSize = max_size_gb * 1024^3)

  message("Reading Seurat object: ", input_seurat)
  seu <- readRDS(input_seurat)
  seu <- AbelR::filter_seurat_qc(
    seu = seu,
    min_nFeature_RNA = as.numeric(get_arg(args, "min_nFeature_RNA", "300")),
    max_nFeature_RNA = as.numeric(get_arg(args, "max_nFeature_RNA", "9000")),
    max_percent_mt = as.numeric(get_arg(args, "max_percent_mt", "12.5")),
    min_nCount_RNA = as.numeric(get_arg(args, "min_nCount_RNA", "800")),
    max_nCount_RNA = as.numeric(get_arg(args, "max_nCount_RNA", "60000")),
    doublet_col = get_arg(args, "doublet_col", "scDblFinder.class"),
    doublet_keep = get_arg(args, "doublet_keep", "singlet"),
    use_doublet_filter = as_flag(
      get_arg(args, "use_doublet_filter", "TRUE")
    )
  )
  invisible(gc())

  message("Running SCT and RPCA integration with ", workers, " worker(s).")
  seu <- AbelR::SCT_METHOD_V3(
    seu = seu,
    sample_col = get_arg(args, "sample_col", "SampleID"),
    var_features_n = as.integer(get_arg(args, "var_features_n", "3000")),
    sct_method = get_arg(args, "sct_method", "glmGamPoi"),
    npcs = as.integer(get_arg(args, "npcs", "30")),
    dims = parse_dims(get_arg(args, "dims", "1:30")),
    k.anchor = as.integer(get_arg(args, "k_anchor", "20")),
    k.weight = as.integer(get_arg(args, "k_weight", "20")),
    resolution = as.numeric(get_arg(args, "resolution", "0.5")),
    umap_n_neighbors = as.integer(
      get_arg(args, "umap_n_neighbors", "30")
    ),
    umap_min_dist = as.numeric(get_arg(args, "umap_min_dist", "0.3")),
    umap_spread = as.numeric(get_arg(args, "umap_spread", "1.0")),
    join_layers = as_flag(get_arg(args, "join_layers", "TRUE")),
    maxSize_GB = max_size_gb,
    verbose = as_flag(get_arg(args, "verbose", "FALSE"))
  )

  ensure_parent(out_seurat)
  saveRDS(seu, out_seurat)
  message("Saved integrated Seurat object: ", out_seurat)
}

main()
