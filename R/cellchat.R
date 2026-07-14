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
  data.input <- Seurat::GetAssayData(
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
  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels"
  )
  if (species == "human") {
    cellchat@DB <- getExportedValue("CellChat", "CellChatDB.human")
  } else {
    cellchat@DB <- getExportedValue("CellChat", "CellChatDB.mouse")
  }
  cellchat <- CellChat::subsetData(cellchat)
  old_max_size <- getOption("future.globals.maxSize")
  old_plan <- future::plan()
  on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
  on.exit(future::plan(old_plan), add = TRUE)
  if (use_parallel) {
    options(future.globals.maxSize = maxSize)
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = type
  )
  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = min.cells
  )
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat <- CellChat::netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )
  return(cellchat)
}

