# Functions for CytoTRACE2 analyses.

#' Run CytoTRACE2 separately for each sample
#'
#' Splits a Seurat object by sample, runs CytoTRACE2 independently for each
#' sample, and combines the resulting CytoTRACE2 metadata. Full per-sample
#' Seurat result objects are retained only when requested to reduce memory use
#' during long server jobs.
#'
#' @param seu A Seurat object containing the expression assay and sample
#'   metadata.
#' @param sample_col Metadata column identifying samples.
#' @param assay Seurat assay used by CytoTRACE2.
#' @param slot_type Expression slot or layer requested by CytoTRACE2, commonly
#'   `"counts"`.
#' @param species Species accepted by CytoTRACE2; `"human"` or `"mouse"`.
#' @param ncores Number of cores passed to CytoTRACE2 for each sample.
#' @param subset_col Optional metadata column used to subset cells before the
#'   per-sample analysis.
#' @param subset_values Optional character vector of values retained from
#'   `subset_col`.
#' @param keep_result_list Logical; retain complete CytoTRACE2 Seurat objects in
#'   the returned `results` component.
#' @param verbose Logical; report samples, cell counts, and progress.
#'
#' @return A list with `metadata`, the combined cell-level CytoTRACE2 metadata,
#'   and `results`, the optional named list of complete per-sample results.
#' @export
run_cytotrace2_by_sample <- function(
  seu,
  sample_col = "SampleID",
  assay = "RNA",
  slot_type = "counts",
  species = c("human", "mouse"),
  ncores = 4,
  subset_col = NULL,
  subset_values = NULL,
  keep_result_list = FALSE,
  verbose = TRUE
) {
  species <- match.arg(species)
  if (!requireNamespace("CytoTRACE2", quietly = TRUE)) {
    stop("Package 'CytoTRACE2' is required for this analysis.")
  }
  if (!inherits(seu, "Seurat")) {
    stop("seu must be a Seurat object.")
  }
  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("sample_col not found in seu@meta.data: ", sample_col)
  }
  if (!assay %in% names(seu@assays)) {
    stop("assay not found in Seurat object: ", assay)
  }
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores < 1) {
    stop("ncores must be a positive number.")
  }
  if (xor(is.null(subset_col), is.null(subset_values))) {
    stop("subset_col and subset_values must be supplied together.")
  }

  SeuratObject::DefaultAssay(seu) <- assay
  if (!is.null(subset_col)) {
    if (!subset_col %in% colnames(seu@meta.data)) {
      stop("subset_col not found in seu@meta.data: ", subset_col)
    }
    cells_use <- rownames(seu@meta.data)[
      seu@meta.data[[subset_col]] %in% subset_values
    ]
    if (!length(cells_use)) {
      stop(
        "No cells retained for ",
        subset_col,
        " in: ",
        paste(subset_values, collapse = ", ")
      )
    }
    seu <- subset(seu, cells = cells_use)
  }

  samples <- unique(as.character(seu@meta.data[[sample_col]]))
  samples <- samples[!is.na(samples) & nzchar(samples)]
  if (!length(samples)) {
    stop("No valid sample identifiers found in metadata column: ", sample_col)
  }
  if (verbose) {
    message("CytoTRACE2 samples (n = ", length(samples), "): ")
    message(paste(samples, collapse = ", "))
  }

  metadata_list <- vector("list", length(samples))
  names(metadata_list) <- samples
  result_list <- if (isTRUE(keep_result_list)) {
    vector("list", length(samples))
  } else {
    NULL
  }
  if (!is.null(result_list)) {
    names(result_list) <- samples
  }

  for (sample_id in samples) {
    cells <- rownames(seu@meta.data)[
      as.character(seu@meta.data[[sample_col]]) == sample_id
    ]
    if (verbose) {
      message(
        "Running CytoTRACE2: ",
        sample_id,
        " (",
        length(cells),
        " cells)"
      )
    }
    seu_sub <- subset(seu, cells = cells)
    SeuratObject::DefaultAssay(seu_sub) <- assay
    result <- CytoTRACE2::cytotrace2(
      seu_sub,
      is_seurat = TRUE,
      slot_type = slot_type,
      species = species,
      ncores = as.integer(ncores)
    )
    cols_use <- grep("CytoTRACE2", colnames(result@meta.data), value = TRUE)
    if (!length(cols_use)) {
      stop("No CytoTRACE2 metadata columns produced for sample: ", sample_id)
    }
    sample_metadata <- result@meta.data[, cols_use, drop = FALSE]
    sample_metadata[[sample_col]] <- sample_id
    sample_metadata$cell_barcode <- rownames(sample_metadata)
    metadata_list[[sample_id]] <- sample_metadata

    if (isTRUE(keep_result_list)) {
      result_list[[sample_id]] <- result
    }
    rm(seu_sub, result, sample_metadata)
    invisible(gc())
  }

  metadata <- do.call(rbind, metadata_list)
  list(metadata = metadata, results = result_list)
}
