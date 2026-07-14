# Functions for monocle3 analyses.

#' Build a Monocle3 trajectory from a Seurat object
#'
#' Creates a Monocle3 cell-data set from a Seurat assay, preprocesses and
#' optionally aligns it, reuses an existing Seurat UMAP embedding, learns a
#' principal graph, and optionally orders cells from selected root clusters.
#'
#' @param seu A Seurat object providing counts and cell metadata.
#' @param ref_seu Optional Seurat object providing the reference embedding. If
#'   `NULL`, `seu` is used.
#' @param assay Assay from which expression values are extracted.
#' @param layer Assay layer containing the count matrix.
#' @param reduction Name of the Seurat dimensional reduction to reuse.
#' @param alignment_group Metadata column used by `monocle3::align_cds()`. Set
#'   to `NULL` to skip alignment.
#' @param seurat_cluster_col Metadata column containing cluster labels.
#' @param root_clusters Optional cluster labels used to select trajectory root
#'   cells. If `NULL`, cells are not ordered.
#' @param num_dim Number of dimensions used during Monocle3 preprocessing.
#' @param use_partition Logical passed to `monocle3::learn_graph()`.
#' @param cluster_cells_first Logical; run `monocle3::cluster_cells()` before
#'   graph learning.
#'
#' @return A Monocle3 `cell_data_set` object.
#' @export
run_monocle3_from_seurat_umap <- function(
    seu,
    ref_seu = NULL,
    assay = "SCT",
    layer = "counts",
    reduction = "umap",
    alignment_group = "SampleID",
    seurat_cluster_col = "seurat_clusters",
    root_clusters = NULL,
    num_dim = 50,
    use_partition = FALSE,
    cluster_cells_first = TRUE
) {
  for (pkg in c("monocle3", "SingleCellExperiment", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for this analysis.")
    }
  }
  # 1. prepare counts
  counts <- Seurat::GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )
  gene_anno <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts)
  )
  cell_meta <- seu@meta.data
  # make sure rownames of metadata match cells
  cell_meta <- cell_meta[colnames(counts), , drop = FALSE]
  cds <- monocle3::new_cell_data_set(
    counts,
    cell_metadata = cell_meta,
    gene_metadata = gene_anno
  )
  # 2. preprocess and align
  cds <- monocle3::preprocess_cds(
    cds,
    num_dim = num_dim
  )
  if (!is.null(alignment_group)) {
    if (!alignment_group %in% colnames(SummarizedExperiment::colData(cds))) {
      stop("alignment_group not found in cell metadata: ", alignment_group)
    }
    cds <- monocle3::align_cds(
      cds,
      alignment_group = alignment_group
    )
  }
  # 3. use external Seurat UMAP
  if (is.null(ref_seu)) {
    ref_seu <- seu
  }
  if (!reduction %in% Seurat::Reductions(ref_seu)) {
    stop("Reduction not found in ref_seu: ", reduction)
  }
  umap_coord <- Seurat::Embeddings(ref_seu, reduction = reduction)
  missing_cells <- setdiff(colnames(cds), rownames(umap_coord))
  if (length(missing_cells) > 0) {
    stop(
      "Some cells in cds are missing from reference UMAP coordinates. n = ",
      length(missing_cells)
    )
  }
  umap_coord <- umap_coord[colnames(cds), , drop = FALSE]
  SingleCellExperiment::reducedDims(cds)$UMAP <- umap_coord
  # 4. cluster and learn graph
  if (cluster_cells_first) {
    cds <- monocle3::cluster_cells(
      cds,
      reduction_method = "UMAP"
    )
  }
  cds <- monocle3::learn_graph(
    cds,
    use_partition = use_partition
  )
  # 5. order cells
  if (!is.null(root_clusters)) {
    if (!seurat_cluster_col %in% colnames(SummarizedExperiment::colData(cds))) {
      stop("seurat_cluster_col not found in metadata: ", seurat_cluster_col)
    }
    root_cells <- colnames(cds)[
      as.character(SummarizedExperiment::colData(cds)[[seurat_cluster_col]]) %in%
        as.character(root_clusters)
    ]
    if (length(root_cells) == 0) {
      stop("No root cells found. Check root_clusters and seurat_cluster_col.")
    }
    cds <- monocle3::order_cells(
      cds,
      root_cells = root_cells
    )
  }
  return(cds)
}
