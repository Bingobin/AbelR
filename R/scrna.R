# Functions for scrna analyses.

#' Filter a Seurat object using common RNA quality-control metrics
#'
#' Retains cells within configurable feature, count, and mitochondrial-content
#' thresholds and optionally restricts a doublet-classification column to one
#' requested value.
#'
#' @param seu A Seurat object containing RNA quality-control metadata.
#' @param min_nFeature_RNA,max_nFeature_RNA Exclusive lower and upper bounds for
#'   detected RNA features.
#' @param max_percent_mt Exclusive upper bound for mitochondrial percentage.
#' @param min_nCount_RNA,max_nCount_RNA Exclusive lower and upper RNA-count
#'   bounds.
#' @param doublet_col Metadata column containing doublet classifications.
#' @param doublet_keep Value retained from `doublet_col`.
#' @param use_doublet_filter Logical; apply the doublet classification filter.
#' @param verbose Logical; report cell counts before and after filtering.
#'
#' @return A Seurat object containing only cells that pass all enabled filters.
#' @export
filter_seurat_qc <- function(
  seu,
  min_nFeature_RNA = 300,
  max_nFeature_RNA = 9000,
  max_percent_mt = 12.5,
  min_nCount_RNA = 800,
  max_nCount_RNA = 60000,
  doublet_col = "scDblFinder.class",
  doublet_keep = "singlet",
  use_doublet_filter = TRUE,
  verbose = TRUE
) {
  if (!inherits(seu, "Seurat")) {
    stop("seu must be a Seurat object.")
  }
  required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))
  if (length(missing_cols)) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (isTRUE(use_doublet_filter) && !doublet_col %in% colnames(seu@meta.data)) {
    stop("Doublet metadata column not found: ", doublet_col)
  }

  meta <- seu@meta.data
  keep <- meta$nFeature_RNA > min_nFeature_RNA &
    meta$nFeature_RNA < max_nFeature_RNA &
    meta$percent.mt < max_percent_mt &
    meta$nCount_RNA > min_nCount_RNA &
    meta$nCount_RNA < max_nCount_RNA
  keep[is.na(keep)] <- FALSE
  if (isTRUE(use_doublet_filter)) {
    keep <- keep & !is.na(meta[[doublet_col]]) &
      meta[[doublet_col]] == doublet_keep
  }
  cells_keep <- rownames(meta)[keep]
  if (!length(cells_keep)) {
    stop("No cells retained after QC filtering.")
  }
  if (verbose) {
    message(
      "Seurat QC retained ",
      length(cells_keep),
      " of ",
      ncol(seu),
      " cells."
    )
  }
  subset(seu, cells = cells_keep)
}


#' Run Seurat v5 SCT and RPCA integration
#'
#' Normalizes a layered Seurat object with SCTransform, performs PCA and RPCA
#' layer integration, constructs a shared-nearest-neighbour graph, clusters the
#' cells, and calculates a UMAP embedding.
#'
#' @param seu A Seurat object containing a sample-identity metadata column.
#' @param sample_col Metadata column used to split RNA layers by sample.
#' @param var_features_n Number of variable features retained by SCTransform.
#' @param sct_method SCTransform fitting method, such as `"glmGamPoi"`.
#' @param npcs Number of principal components to calculate.
#' @param dims Principal-component dimensions used for integration, neighbours,
#'   and UMAP. `max(dims)` must not exceed `npcs`.
#' @param k.anchor Number of anchors used by RPCA integration.
#' @param k.weight Number of neighbours used when weighting anchors.
#' @param resolution Clustering resolution.
#' @param umap_n_neighbors Number of UMAP neighbours.
#' @param umap_min_dist Minimum UMAP distance.
#' @param umap_spread UMAP spread parameter.
#' @param join_layers Logical; join existing RNA layers before splitting them by
#'   `sample_col`.
#' @param maxSize_GB Maximum future global size in gigabytes.
#' @param verbose Logical; show progress from Seurat functions.
#'
#' @return The processed Seurat object with an `integrated.dr` reduction,
#'   clusters, and UMAP coordinates.
#' @export
scRNA_SCT_norm <- function(
    seu,
    sample_col = "SampleID",
    # SCTransform / PCA
    var_features_n = 3000,
    sct_method = "glmGamPoi",
    npcs = 30,
    
    # Integration / graph
    dims = 1:30,
    k.anchor = 20,
    k.weight = 20,
    resolution = 0.5,
    
    # UMAP (影响“分得太开”的视觉表现很大)
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    umap_spread = 1.0,
    
    # housekeeping
    join_layers = TRUE,
    maxSize_GB = 35,
    verbose = FALSE
){
  old_max_size <- getOption("future.globals.maxSize")
  on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
  options(future.globals.maxSize = maxSize_GB * 1024^3)

  if (max(dims) > npcs) {
    stop("max(dims) must be less than or equal to npcs.")
  }

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("sample_col not found in metadata: ", sample_col)
  }
  if (identical(sct_method, "glmGamPoi") &&
      !requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop("Package 'glmGamPoi' is required when sct_method = 'glmGamPoi'.")
  }
  sample_values <- factor(seu@meta.data[[sample_col]])
  if (anyNA(sample_values)) {
    stop("sample_col contains missing values: ", sample_col)
  }
  seu[[sample_col]] <- sample_values
  
  # RNA layers: join then split by sample
  if (join_layers) {
    seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
  }
  seu[["RNA"]] <- split(seu[["RNA"]], f = sample_values)
  
  # SCT
  seu <- SCTransform(
    seu,
    variable.features.n = var_features_n,
    method = sct_method,
    assay = "RNA",
    verbose = verbose
  )
  
  DefaultAssay(seu) <- "SCT"
  
  # PCA
  seu <- RunPCA(seu, npcs = npcs, verbose = verbose, assay = "SCT")
  
  # Integration
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
  
  # Neighbors / clusters
  seu <- FindNeighbors(seu, reduction = "integrated.dr", dims = dims, graph.name = "integrated_snn")
  seu <- FindClusters(seu, graph.name = "integrated_snn", resolution = resolution)
  
  # UMAP
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


#' Annotate Seurat clusters with scType
#'
#' Loads the scType marker database and scoring scripts, scores positive and
#' negative marker sets, selects the highest-scoring cell type for each Seurat
#' cluster, and stores the annotation in `customclassif` metadata.
#'
#' @param scRNA_object A Seurat object containing scaled assay data and a
#'   `seurat_clusters` metadata column.
#' @param assay Assay whose `scale.data` layer is used for scoring.
#' @param tissue Tissue category present in the scType marker database.
#' @param sctype_dir Directory containing `gene_sets_prepare.R`,
#'   `sctype_score_.R`, and, by default, `ScTypeDB_full.xlsx`. It may also be set
#'   with `options(AbelR.sctype_dir = "/path/to/scTYPE")`.
#' @param db_file Optional explicit path to the scType marker database file.
#'
#' @return A list containing the annotated `scRNA_object`, the colour vector
#'   `ccolss`, per-cluster `sctype_scores`, and the complete `CL_results` table.
#' @export
scTYPE_annotation <- function(
    scRNA_object,
    assay,
    tissue = "Immune system",
    sctype_dir = getOption("AbelR.sctype_dir"),
    db_file = NULL
) {
  #  scRNA_object <- pbmc.combined.sct
  #  assay <- "integrated"
  # tissue = "Immune system"

  ccolss <- c(
    "#5f75ae",
    "#92bbb8",
    "#64a841",
    "#e5486e",
    "#de8e06",
    "#eccf5a",
    "#b5aa0f",
    "#e4b680",
    "#7ba39d",
    "#b15928",
    "#ffff99",
    "#6a3d9a",
    "#cab2d6",
    "#ff7f00",
    "#fdbf6f",
    "#e31a1c",
    "#fb9a99",
    "#33a02c",
    "#b2df8a",
    "#1f78b4",
    "#a6cee3"
  )
  if (is.null(sctype_dir) || !nzchar(sctype_dir)) {
    stop(
      "Set sctype_dir or options(AbelR.sctype_dir = '/path/to/scTYPE')."
    )
  }
  prepare_file <- file.path(sctype_dir, "gene_sets_prepare.R")
  score_file <- file.path(sctype_dir, "sctype_score_.R")
  if (is.null(db_file)) {
    db_file <- file.path(sctype_dir, "ScTypeDB_full.xlsx")
  }
  required_files <- c(prepare_file, score_file, db_file)
  if (any(!file.exists(required_files))) {
    stop(
      "Missing scType file(s): ",
      paste(required_files[!file.exists(required_files)], collapse = ", ")
    )
  }
  source(prepare_file, local = environment())
  source(score_file, local = environment())

  if (!assay %in% names(scRNA_object@assays)) {
    stop("Assay not found in scRNA_object: ", assay)
  }
  gs_list <- gene_sets_prepare(db_file, tissue)
  scale_data <- Seurat::GetAssayData(
    scRNA_object,
    assay = assay,
    layer = "scale.data"
  )
  ex.max <- sctype_score(
    scRNAseqData = scale_data,
    scale = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  CL_results <- do.call(
    "rbind",
    lapply(unique(scRNA_object@meta.data$seurat_clusters), function(cl) {
      es.max.cl <- sort(
        rowSums(ex.max[, rownames(scRNA_object@meta.data[
          scRNA_object@meta.data$seurat_clusters == cl,
        ])]),
        decreasing = !0
      )
      head(
        data.frame(
          cluster = cl,
          type = names(es.max.cl),
          scores = es.max.cl,
          ncells = sum(scRNA_object@meta.data$seurat_clusters == cl)
        ),
        10
      )
    })
  )
  sctype_scores <- CL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)
  sctype_scores <- sctype_scores[!duplicated(sctype_scores$cluster), ]
  sctype_scores$type[
    as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4
  ] <- "Unknown"
  scRNA_object@meta.data$customclassif <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    scRNA_object@meta.data$customclassif[
      scRNA_object@meta.data$seurat_clusters == j
    ] <-
      paste0(as.character(cl_type$type[1]))
  }

  # lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
  #
  # #prepare edges
  # cL_resutls <- CL_results
  # cL_resutls=cL_resutls[order(cL_resutls$cluster),]
  # edges = cL_resutls
  # edges$type = paste0(edges$type,"_",edges$cluster)
  # edges$cluster = paste0("cluster ", edges$cluster)
  # edges = edges[,c("cluster", "type")]
  # colnames(edges) = c("from", "to")
  # rownames(edges) <- NULL
  #
  # #prepare nodes
  # nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]
  # nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster)
  # nodes_lvl1$Colour = "#f1f1ef"
  # nodes_lvl1$ord = 1
  # nodes_lvl1$realname = nodes_lvl1$cluster
  # nodes_lvl1 = as.data.frame(nodes_lvl1)
  # nodes_lvl2 = c()
  # for (i in 1:length(unique(cL_resutls$cluster))){
  #   dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]
  #   nodes_lvl2 = rbind(nodes_lvl2,
  #                      data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster),
  #                                 ncells = dt_tmp$scores,
  #                                 Colour = ccolss[i], ord = 2,
  #                                 realname = dt_tmp$type))
  # }
  # nodes = rbind(nodes_lvl1, nodes_lvl2)
  # nodes$ncells[nodes$ncells<1] = 1
  # files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]
  # files_db = unique(files_db)
  # nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
  # nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]
  # nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  # mygraph <- graph_from_data_frame(edges, vertices=nodes)
  # ##plot
  # gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
  #   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) +
  #   geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  #   theme_void() +
  #   geom_node_text(aes(filter=ord==2,
  #                      label=shortName, colour=I("#ffffff"),
  #                      fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1))) +
  #   geom_node_label(aes(filter=ord==1,  label=shortName,
  #                       colour=I("#000000"), size = I(3), fill="white", parse = T),
  #                   repel = !0, segment.linetype="dotted")
  result <- list()
  result[["scRNA_object"]] <- scRNA_object
  #  result[["gggr_plot"]] <- gggr
  result[["ccolss"]] <- ccolss
  result[["sctype_scores"]] <- sctype_scores
  result[["CL_results"]] <- CL_results

  return(result)
}


# CytoTRACE2 analysis ------------------------------------------------------

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


# Monocle3 trajectory analysis --------------------------------------------

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


# CellChat analysis -------------------------------------------------------

#' Build and run a CellChat analysis
#'
#' Extracts a normalized expression layer and metadata from a Seurat object,
#' selects the human or mouse CellChat database, identifies overexpressed
#' ligand-receptor genes and interactions, and computes communication networks
#' and pathway centrality.
#'
#' @param seu A Seurat object with normalized expression and cell metadata.
#' @param species Species matching the expression matrix and CellChat database;
#'   either `"human"` or `"mouse"`.
#' @param group.by Metadata column defining the cell groups used by CellChat.
#' @param sample.by Metadata column identifying samples or replicates.
#' @param cluster.by Optional metadata column copied to `meta$clusters`.
#' @param assay Seurat assay from which expression data are extracted.
#' @param layer Normalized expression layer passed to CellChat.
#' @param min.cells Minimum number of cells required for retained communication.
#' @param workers_overexpress Number of future workers used during
#'   overexpression analysis.
#' @param workers_prob Number of future workers used for communication
#'   probability calculations.
#' @param maxSize Maximum allowed future global size in bytes.
#' @param type Averaging method passed to `CellChat::computeCommunProb()`.
#' @param use_parallel Logical; use multisession parallel processing.
#' @param verbose Logical; print input summaries and progress messages.
#'
#' @return A processed CellChat object.
#' @export
Build_CellChat_object <- function(
  seu,
  species = c("human", "mouse"),
  group.by = "CellType",
  sample.by = "SampleID",
  cluster.by = NULL,
  assay = "RNA",
  layer = "data",
  min.cells = 10,
  workers_overexpress = 6,
  workers_prob = 2,
  maxSize = 40 * 1024^3,
  type = "triMean",
  use_parallel = TRUE,
  verbose = TRUE
) {
  species <- match.arg(species)
  for (pkg in c("Seurat", "CellChat", "future")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for Build_CellChat_object().")
    }
  }
  if (!group.by %in% colnames(seu@meta.data)) {
    stop("group.by column not found in seu@meta.data: ", group.by)
  }
  if (!sample.by %in% colnames(seu@meta.data)) {
    stop("sample.by column not found in seu@meta.data: ", sample.by)
  }
  if (!is.null(cluster.by) && !cluster.by %in% colnames(seu@meta.data)) {
    stop("cluster.by column not found in seu@meta.data: ", cluster.by)
  }
  data.input <- GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )
  meta <- data.frame(
    labels = seu@meta.data[[group.by]],
    samples = seu@meta.data[[sample.by]],
    row.names = colnames(seu)
  )
  if (!is.null(cluster.by)) {
    meta$clusters <- seu@meta.data[[cluster.by]]
  }
  meta$labels <- as.factor(meta$labels)
  meta$samples <- as.factor(meta$samples)
  if (verbose) {
    message("CellChat grouping column: ", group.by)
    message("Sample column: ", sample.by)
    message("Number of cells: ", ncol(data.input))
    message("Number of groups: ", length(unique(meta$labels)))
    print(table(meta$labels))
  }
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels"
  )
  if (species == "human") {
    cellchat@DB <- CellChatDB.human
  } else {
    cellchat@DB <- CellChatDB.mouse
  }
  message("1. subsetData")
  print(system.time({
    cellchat <- subsetData(cellchat)
  }))

  options(future.globals.maxSize = maxSize)
  if (use_parallel) {
    future::plan(future::multisession, workers = workers_overexpress)
  } else {
    future::plan(future::sequential)
  }
  message("2. identifyOverExpressedGenes")
  print(system.time({
    cellchat <- identifyOverExpressedGenes(cellchat)
  }))

  message("3. identifyOverExpressedInteractions")
  print(system.time({
    cellchat <- identifyOverExpressedInteractions(cellchat)
  }))

  if (use_parallel) {
    future::plan(future::multisession, workers = workers_prob)
  } else {
    future::plan(future::sequential)
  }
  message("4. computeCommunProb")
  print(system.time({
    cellchat <- computeCommunProb(
      cellchat,
      type = type
    )
  }))

  message("5. filterCommunication")
  print(system.time({
    cellchat <- filterCommunication(
      cellchat,
      min.cells = min.cells
    )
  }))

  message("6. computeCommunProbPathway")
  print(system.time({
    cellchat <- computeCommunProbPathway(cellchat)
  }))

  message("7. aggregateNet")
  print(system.time({
    cellchat <- aggregateNet(cellchat)
  }))

  message("8. centrality")
  print(system.time({
    cellchat <- netAnalysis_computeCentrality(
      cellchat,
      slot.name = "netP"
    )
  }))

  future::plan(future::sequential)
  return(cellchat)
}
