# Functions for cellchat analyses.

Build_CellChat_object <- function(
    seu,
    species = c("human", "mouse"),
    group.by = "CellType",
    sample.by = "SampleID",
    cluster.by = NULL,
    assay = "RNA",
    layer = "data",
    min.cells = 10,
    workers = 4,
    maxSize = 40 * 1024^3,
    type = "triMean",
    use_parallel = TRUE,
    verbose = TRUE
) {
  species <- match.arg(species)
  suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(future)
  })
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
  cellchat <- subsetData(cellchat)
  if (use_parallel) {
    options(future.globals.maxSize = maxSize)
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = type
  )
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min.cells
  )
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )
  future::plan(future::sequential)
  return(cellchat)
}


