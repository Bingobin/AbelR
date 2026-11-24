test_that("Seurat UMAP labels cells by barcode when the reduction has fewer cells", {
  cells <- paste0("cell", 1:5)
  counts <- matrix(1L, nrow = 2, ncol = 5,
                   dimnames = list(c("gene1", "gene2"), cells))
  seu <- SeuratObject::CreateSeuratObject(counts = SeuratObject::as.sparse(counts))
  seu$cell_annotation <- c("unplotted", "T", NA, "T", "B")
  xy <- cbind(UMAP_1 = c(0, 1, 2, 5), UMAP_2 = c(2, 7, 4, 8))
  rownames(xy) <- cells[-1L]
  seu[["UMAP_dim30"]] <- SeuratObject::CreateDimReducObject(
    embeddings = xy, key = "UMAP_", assay = "RNA"
  )

  plot <- plot_seurat_umap(seu, "UMAP_dim30", "cell_annotation")

  expect_s3_class(plot, "ggplot")
  expect_identical(rownames(plot$data), cells[-1L])
  expect_identical(as.character(plot$data$cell_type), c("T", NA, "T", "B"))
  centers <- plot$layers[[2L]]$data
  expect_equal(centers$x[centers$cell_type == "T"], 1)
  expect_equal(centers$y[centers$cell_type == "T"], 3)
  expect_setequal(as.character(centers$cell_type), c("T", "B"))
})
