#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat("
Usage:
  Rscript run_cytotrace2_by_sample.R \\
    --input_seurat Nevus_scRNA.combind.seu_sct.rds \\
    --out_rds Nevus_scRNA.combind.CytoTRACE2_df.rds \\
    --sample_col SampleID \\
    --assay RNA \\
    --slot_type counts \\
    --species human \\
    --ncores 4

Optional:
    --subset_col CellType
    --subset_values Melanocyte
    --save_result_list TRUE
    --out_result_list Nevus_scRNA.combind.CytoTRACE2_result_list.rds
")
    quit(save = "no", status = 0)
  }
  
  res <- list()
  i <- 1
  while (i <= length(args)) {
    key <- gsub("^--", "", args[i])
    value <- args[i + 1]
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

args <- parse_args()

input_seurat <- get_arg(args, "input_seurat", required = TRUE)
out_rds <- get_arg(args, "out_rds", required = TRUE)

sample_col <- get_arg(args, "sample_col", default = "SampleID")
assay_use <- get_arg(args, "assay", default = "RNA")
slot_type <- get_arg(args, "slot_type", default = "counts")
species_use <- get_arg(args, "species", default = "human")
ncores_use <- as.integer(get_arg(args, "ncores", default = "4"))

subset_col <- get_arg(args, "subset_col", default = NULL)
subset_values <- get_arg(args, "subset_values", default = NULL)

save_result_list <- get_arg(args, "save_result_list", default = "FALSE")
save_result_list <- toupper(save_result_list) %in% c("TRUE", "T", "YES", "Y", "1")

out_result_list <- get_arg(
  args,
  "out_result_list",
  default = sub("\\.rds$", ".CytoTRACE2_result_list.rds", out_rds)
)

message("========== CytoTRACE2 by sample ==========")
message("Input Seurat RDS: ", input_seurat)
message("Output cytotrace2_df RDS: ", out_rds)
message("Sample column: ", sample_col)
message("Assay: ", assay_use)
message("Slot type: ", slot_type)
message("Species: ", species_use)
message("ncores: ", ncores_use)

seu <- readRDS(input_seurat)

if (!sample_col %in% colnames(seu@meta.data)) {
  stop("sample_col not found in seu@meta.data: ", sample_col)
}

if (!assay_use %in% names(seu@assays)) {
  stop("assay not found in Seurat object: ", assay_use)
}

DefaultAssay(seu) <- assay_use

if (!is.null(subset_col) && !is.null(subset_values)) {
  if (!subset_col %in% colnames(seu@meta.data)) {
    stop("subset_col not found in seu@meta.data: ", subset_col)
  }
  
  subset_values_vec <- unlist(strsplit(subset_values, ","))
  subset_values_vec <- trimws(subset_values_vec)
  
  message("Subset column: ", subset_col)
  message("Subset values: ", paste(subset_values_vec, collapse = ", "))
  
  cells_use <- rownames(seu@meta.data)[seu@meta.data[[subset_col]] %in% subset_values_vec]
  
  if (length(cells_use) == 0) {
    stop("No cells found after subsetting by ", subset_col, " in ", paste(subset_values_vec, collapse = ", "))
  }
  
  seu <- subset(seu, cells = cells_use)
  message("Cells after subset: ", ncol(seu))
}

unique_samples <- unique(as.character(seu@meta.data[[sample_col]]))
unique_samples <- unique_samples[!is.na(unique_samples)]

message("Number of samples: ", length(unique_samples))
message("Samples: ", paste(unique_samples, collapse = ", "))

cytotrace2_results_list <- list()

for (x in unique_samples) {
  message("========== Running CytoTRACE2 for sample: ", x, " ==========")
  
  cells_x <- rownames(seu@meta.data)[as.character(seu@meta.data[[sample_col]]) == x]
  
  if (length(cells_x) == 0) {
    warning("No cells found for sample: ", x)
    next
  }
  
  seu_sub <- subset(seu, cells = cells_x)
  DefaultAssay(seu_sub) <- assay_use
  
  message("Cells in sample ", x, ": ", ncol(seu_sub))
  
  result <- cytotrace2(
    seu_sub,
    is_seurat = TRUE,
    slot_type = slot_type,
    species = species_use,
    ncores = ncores_use
  )
  
  cytotrace2_results_list[[x]] <- result
  
  rm(seu_sub, result)
  gc()
}

message("========== Extracting CytoTRACE2 metadata ==========")

cytotrace2_df_list <- lapply(names(cytotrace2_results_list), function(x) {
  obj <- cytotrace2_results_list[[x]]
  
  cols_use <- grep("CytoTRACE2", colnames(obj@meta.data), value = TRUE)
  
  if (length(cols_use) == 0) {
    warning("No CytoTRACE2 columns found in sample: ", x)
    return(NULL)
  }
  
  df <- obj@meta.data[, cols_use, drop = FALSE]
  df[[sample_col]] <- x
  df$cell_barcode <- rownames(df)
  
  return(df)
})

cytotrace2_df_list <- cytotrace2_df_list[!sapply(cytotrace2_df_list, is.null)]

if (length(cytotrace2_df_list) == 0) {
  stop("No CytoTRACE2 metadata extracted.")
}

cytotrace2_df <- do.call(rbind, cytotrace2_df_list)

message("Final cytotrace2_df dimension:")
print(dim(cytotrace2_df))

saveRDS(cytotrace2_df, file = out_rds)
message("Saved cytotrace2_df to: ", out_rds)

if (save_result_list) {
  saveRDS(cytotrace2_results_list, file = out_result_list)
  message("Saved CytoTRACE2 result list to: ", out_result_list)
}

message("Done.")
