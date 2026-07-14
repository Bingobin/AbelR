# Functions for scrna analyses.

#' Integrate Seurat objects with the SCT workflow
#'
#' Applies SCTransform to a list of Seurat objects, identifies integration
#' features and anchors, integrates the data, and runs scaling, PCA, UMAP,
#' neighbour finding, and clustering.
#'
#' @param object_list A list of Seurat objects to integrate.
#'
#' @return An integrated Seurat object with PCA, UMAP, neighbour graph, and
#'   clustering results.
#' @export
scRNA_SCT_norm <- function(object_list) {
  mut.list <- lapply(X = object_list, FUN = SCTransform)
  features <- SelectIntegrationFeatures(
    object.list = mut.list,
    nfeatures = 3000
  )
  mut.list <- PrepSCTIntegration(
    object.list = mut.list,
    anchor.features = features
  )
  mut.anchors <- FindIntegrationAnchors(
    object.list = mut.list,
    normalization.method = "SCT",
    anchor.features = features
  )
  combined.sct <- IntegrateData(
    anchorset = mut.anchors,
    normalization.method = "SCT"
  )
  DefaultAssay(combined.sct) <- "integrated"
  combined.sct <- ScaleData(combined.sct, verbose = FALSE)
  combined.sct <- RunPCA(combined.sct, npcs = 50, verbose = FALSE)
  combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:40)
  combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:40)
  #  combined.sct <- FindClusters(combined.sct, resolution = c(0.2,0.5,0.8,1,1.5,2))
  combined.sct <- FindClusters(combined.sct, resolution = c(0.5, 1))
  return(combined.sct)
}


#' Run Seurat v5 SCT and RPCA integration
#'
#' Normalizes a layered Seurat object with SCTransform, performs PCA and RPCA
#' layer integration, constructs a shared-nearest-neighbour graph, clusters the
#' cells, and calculates a UMAP embedding.
#'
#' @param seu A Seurat object containing a `SampleID` metadata column.
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
#'   `SampleID`.
#' @param maxSize_GB Maximum future global size in gigabytes.
#' @param verbose Logical; show progress from Seurat functions.
#'
#' @return The processed Seurat object with an `integrated.dr` reduction,
#'   clusters, and UMAP coordinates.
#' @export
SCT_METHOD_V3 <- function(
    seu,
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
  
  # ensure SampleID
  seu$SampleID <- factor(seu$SampleID)
  
  # RNA layers: join then split by SampleID
  if (join_layers) {
    seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
  }
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$SampleID)
  
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
