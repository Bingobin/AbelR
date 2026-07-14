#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat("
Usage:
  Rscript run_sct_integrate.R \\
    --input_seurat BMO_scRNA_26609.seu_obj.rds \\
    --out_seurat BMO_scRNA_26609.seu_sct.rds \\
    --sample_col SampleID \\
    --maxSize_GB 120

Required:
  --input_seurat        Input Seurat RDS file
  --out_seurat          Output integrated Seurat RDS file

Optional QC:
  --min_nFeature_RNA    Default: 300
  --max_nFeature_RNA    Default: 9000
  --max_percent_mt      Default: 12.5
  --min_nCount_RNA      Default: 800
  --max_nCount_RNA      Default: 60000
  --doublet_col         Default: scDblFinder.class
  --doublet_keep        Default: singlet
  --use_doublet_filter  Default: TRUE

Optional SCT/integration:
  --sample_col          Default: SampleID
  --var_features_n      Default: 3000
  --sct_method          Default: glmGamPoi
  --npcs                Default: 30
  --dims                Default: 1:30
  --k_anchor            Default: 20
  --k_weight            Default: 20
  --resolution          Default: 0.5
  --umap_n_neighbors    Default: 30
  --umap_min_dist       Default: 0.3
  --umap_spread         Default: 1.0
  --join_layers         Default: TRUE
  --maxSize_GB          Default: 35
  --sequential          Default: TRUE
  --verbose             Default: FALSE

Optional files:
  --sample_design       Optional sample design txt file. Read but not required for running.
")
    quit(save = "no", status = 0)
  }
  
  res <- list()
  i <- 1
  while (i <= length(args)) {
    key <- gsub("^--", "", args[i])
    value <- args[i + 1]
    if (is.na(value) || grepl("^--", value)) {
      stop("Argument ", args[i], " has no value.")
    }
    res[[key]] <- value
    i <- i + 2
  }
  
  return(res)
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  if (!is.null(args[[name]])) {
    return(args[[name]])
  }
  if (required) {
    stop("Missing required argument: --", name)
  }
  return(default)
}

as_logical_arg <- function(x) {
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "Y", "1")
}

parse_dims <- function(x) {
  x <- gsub(" ", "", x)
  
  if (grepl(":", x)) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    return(seq(as.integer(parts[1]), as.integer(parts[2])))
  }
  
  as.integer(strsplit(x, ",")[[1]])
}

args <- parse_args()

input_seurat <- get_arg(args, "input_seurat", required = TRUE)
out_seurat <- get_arg(args, "out_seurat", required = TRUE)

sample_design <- get_arg(args, "sample_design", default = NULL)

min_nFeature_RNA <- as.numeric(get_arg(args, "min_nFeature_RNA", default = "300"))
max_nFeature_RNA <- as.numeric(get_arg(args, "max_nFeature_RNA", default = "9000"))
max_percent_mt <- as.numeric(get_arg(args, "max_percent_mt", default = "12.5"))
min_nCount_RNA <- as.numeric(get_arg(args, "min_nCount_RNA", default = "800"))
max_nCount_RNA <- as.numeric(get_arg(args, "max_nCount_RNA", default = "60000"))

doublet_col <- get_arg(args, "doublet_col", default = "scDblFinder.class")
doublet_keep <- get_arg(args, "doublet_keep", default = "singlet")
use_doublet_filter <- as_logical_arg(get_arg(args, "use_doublet_filter", default = "TRUE"))

sample_col <- get_arg(args, "sample_col", default = "SampleID")
var_features_n <- as.integer(get_arg(args, "var_features_n", default = "3000"))
sct_method <- get_arg(args, "sct_method", default = "glmGamPoi")
npcs <- as.integer(get_arg(args, "npcs", default = "30"))
dims <- parse_dims(get_arg(args, "dims", default = "1:30"))

k_anchor <- as.integer(get_arg(args, "k_anchor", default = "20"))
k_weight <- as.integer(get_arg(args, "k_weight", default = "20"))
resolution <- as.numeric(get_arg(args, "resolution", default = "0.5"))

umap_n_neighbors <- as.integer(get_arg(args, "umap_n_neighbors", default = "30"))
umap_min_dist <- as.numeric(get_arg(args, "umap_min_dist", default = "0.3"))
umap_spread <- as.numeric(get_arg(args, "umap_spread", default = "1.0"))

join_layers <- as_logical_arg(get_arg(args, "join_layers", default = "TRUE"))
maxSize_GB <- as.numeric(get_arg(args, "maxSize_GB", default = "35"))
sequential <- as_logical_arg(get_arg(args, "sequential", default = "TRUE"))
verbose <- as_logical_arg(get_arg(args, "verbose", default = "FALSE"))

message("========== SCT integration pipeline ==========")
message("Input Seurat: ", input_seurat)
message("Output Seurat: ", out_seurat)
message("Sample column: ", sample_col)
message("QC: nFeature_RNA > ", min_nFeature_RNA,
        ", nFeature_RNA < ", max_nFeature_RNA,
        ", percent.mt < ", max_percent_mt,
        ", nCount_RNA > ", min_nCount_RNA,
        ", nCount_RNA < ", max_nCount_RNA)
message("Use doublet filter: ", use_doublet_filter)
if (use_doublet_filter) {
  message("Doublet column: ", doublet_col)
  message("Doublet keep value: ", doublet_keep)
}
message("SCT variable features: ", var_features_n)
message("SCT method: ", sct_method)
message("npcs: ", npcs)
message("dims: ", paste(dims, collapse = ","))
message("k.anchor: ", k_anchor)
message("k.weight: ", k_weight)
message("resolution: ", resolution)
message("UMAP n.neighbors: ", umap_n_neighbors)
message("UMAP min.dist: ", umap_min_dist)
message("UMAP spread: ", umap_spread)
message("JoinLayers: ", join_layers)
message("future.globals.maxSize GB: ", maxSize_GB)
message("Sequential plan: ", sequential)
message("Verbose: ", verbose)

if (sequential) {
  plan(sequential)
}

options(future.globals.maxSize = maxSize_GB * 1024^3)

if (!is.null(sample_design)) {
  message("Reading sample design: ", sample_design)
  sample.design.df <- read.table(sample_design, header = TRUE, sep = "\t", check.names = FALSE)
  message("Sample design dimension: ", paste(dim(sample.design.df), collapse = " x "))
}

message("Reading Seurat object...")
seu_obj <- readRDS(input_seurat)

required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt", sample_col)
missing_cols <- setdiff(required_cols, colnames(seu_obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

if (use_doublet_filter && !doublet_col %in% colnames(seu_obj@meta.data)) {
  stop("Doublet metadata column not found: ", doublet_col)
}

message("Original object:")
print(seu_obj)

message("Cells before QC: ", ncol(seu_obj))

qc_keep <- seu_obj$nFeature_RNA > min_nFeature_RNA &
  seu_obj$nFeature_RNA < max_nFeature_RNA &
  seu_obj$percent.mt < max_percent_mt &
  seu_obj$nCount_RNA > min_nCount_RNA &
  seu_obj$nCount_RNA < max_nCount_RNA

if (use_doublet_filter) {
  qc_keep <- qc_keep & seu_obj@meta.data[[doublet_col]] == doublet_keep
}

cells_keep <- colnames(seu_obj)[qc_keep]

message("Cells after QC: ", length(cells_keep))

if (length(cells_keep) == 0) {
  stop("No cells retained after QC filtering.")
}

seu_obj_f <- subset(seu_obj, cells = cells_keep)

rm(seu_obj)
gc()

SCT_METHOD_V3 <- function(
  seu,
  sample_col = "SampleID",
  var_features_n = 3000,
  sct_method = "glmGamPoi",
  npcs = 30,
  dims = 1:30,
  k.anchor = 20,
  k.weight = 20,
  resolution = 0.5,
  umap_n_neighbors = 30,
  umap_min_dist = 0.3,
  umap_spread = 1.0,
  join_layers = TRUE,
  maxSize_GB = 35,
  verbose = FALSE
) {
  options(future.globals.maxSize = maxSize_GB * 1024^3)
  
  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("sample_col not found in metadata: ", sample_col)
  }
  
  seu[[sample_col]] <- factor(seu[[sample_col]][, 1])
  
  if (join_layers) {
    message("Joining RNA layers...")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }
  
  message("Splitting RNA assay by ", sample_col, "...")
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu[[sample_col]][, 1])
  
  message("Running SCTransform...")
  seu <- SCTransform(
    seu,
    variable.features.n = var_features_n,
    method = sct_method,
    assay = "RNA",
    verbose = verbose
  )
  
  DefaultAssay(seu) <- "SCT"
  
  message("Running PCA...")
  seu <- RunPCA(
    seu,
    npcs = max(npcs, max(dims)),
    verbose = verbose,
    assay = "SCT"
  )
  
  message("Running RPCA integration...")
  seu <- IntegrateLayers(
    object = seu,
    method = RPCAIntegration,
    normalization.method = "SCT",
    orig.reduction = "pca",
    new.reduction = "integrated.dr",
    verbose = verbose,
    k.anchor = k.anchor,
    k.weight = k.weight,
    dims = dims
  )
  
  message("Finding neighbors...")
  seu <- FindNeighbors(
    seu,
    reduction = "integrated.dr",
    dims = dims,
    graph.name = "integrated_snn"
  )
  
  message("Finding clusters...")
  seu <- FindClusters(
    seu,
    graph.name = "integrated_snn",
    resolution = resolution
  )
  
  message("Running UMAP...")
  seu <- RunUMAP(
    seu,
    reduction = "integrated.dr",
    dims = dims,
    n.neighbors = umap_n_neighbors,
    min.dist = umap_min_dist,
    spread = umap_spread,
    verbose = verbose
  )
  
  return(seu)
}

message("Running SCT integration...")
seu_sct <- SCT_METHOD_V3(
  seu = seu_obj_f,
  sample_col = sample_col,
  var_features_n = var_features_n,
  sct_method = sct_method,
  npcs = npcs,
  dims = dims,
  k.anchor = k_anchor,
  k.weight = k_weight,
  resolution = resolution,
  umap_n_neighbors = umap_n_neighbors,
  umap_min_dist = umap_min_dist,
  umap_spread = umap_spread,
  join_layers = join_layers,
  maxSize_GB = maxSize_GB,
  verbose = verbose
)

message("Saving Seurat object...")
saveRDS(seu_sct, file = out_seurat)

message("Final object:")
print(seu_sct)

message("Done.")
