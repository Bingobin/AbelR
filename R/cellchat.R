# Functions for cellchat analyses.

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
