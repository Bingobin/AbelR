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

  withr::local_seed(42L)
  rng_before <- .Random.seed
  shuffled <- plot_seurat_umap(
    seu, "UMAP_dim30", "cell_annotation",
    pt.size = 0.5, alpha = 0.6, shuffle = TRUE, seed = 7L
  )
  shuffled_again <- plot_seurat_umap(
    seu, "UMAP_dim30", "cell_annotation", shuffle = TRUE, seed = 7L
  )
  expect_identical(.Random.seed, rng_before)
  expect_identical(shuffled$data, shuffled_again$data)
  expect_false(identical(rownames(shuffled$data), rownames(plot$data)))
  expect_identical(shuffled$data[rownames(plot$data), ], plot$data)
  expect_identical(shuffled$layers[[2L]]$data, centers)
  drawn_points <- ggplot2::ggplot_build(shuffled)$data[[1L]]
  expect_equal(drawn_points$x, shuffled$data$x)
  expect_true(all(drawn_points$size == 0.5))
  expect_true(all(drawn_points$alpha == 0.6))
})

test_that("Seurat feature points use continuous values and highlight independently of labels", {
  cells <- paste0("cell", 1:5)
  counts <- rbind(CD44 = c(0, 2, 4, 8, 5), MYC = rep(1, 5))
  colnames(counts) <- cells
  seu <- SeuratObject::CreateSeuratObject(SeuratObject::as.sparse(counts))
  seu$score <- c(-1, 0, 2, 10, NA)
  seu$seurat_clusters <- c(0, 1, 0, 1, 0)
  seu$celltype <- c("T", "B", "T", "B", "T")
  xy <- cbind(UMAP_1 = 0:4, UMAP_2 = c(2, 1, 0, 1, 3))
  rownames(xy) <- cells
  seu[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = xy, key = "UMAP_", assay = "RNA"
  )

  gene_plot <- plot_seurat_feature(seu, "CD44", layer = "counts")
  expect_equal(gene_plot$data$feature_value, c(0, 2, 4, 8, 5))
  expect_equal(nrow(gene_plot$layers[[2L]]$data), 4L)
  expect_identical(rownames(gene_plot$layers[[1L]]$data), cells[1])
  expect_identical(gene_plot$labels$title, "CD44 \u00b7 n = 4 / 5")
  expect_setequal(gene_plot$layers[[3L]]$data$cell_type, c("C0", "C1"))

  layers <- ggplot2::ggplot_build(gene_plot)$data
  expect_true(all(layers[[1L]]$colour == "grey85"))
  expect_equal(layers[[2L]]$x, unname(xy[c(2, 3, 5, 4), 1L]))

  all_scores <- plot_seurat_feature(seu, "score", min_cutoff = 2)
  expect_identical(all_scores$labels$title, "score \u00b7 n = 1 / 5")
  expect_identical(rownames(all_scores$layers[[2L]]$data), cells[4])

  score_plot <- plot_seurat_feature(
    seu, "score", group_by = "celltype", highlight_by = "seurat_clusters",
    highlight = 0
  )
  expect_setequal(rownames(score_plot$layers[[1L]]$data), cells[c(1, 2, 4, 5)])
  expect_identical(score_plot$labels$title, "score \u00b7 C0 \u00b7 n = 1 / 3")
  scores <- score_plot$layers[[2L]]$data
  expect_identical(rownames(scores), cells[3])
  expect_equal(scores$feature_value, 2)
  expect_setequal(as.character(score_plot$layers[[3L]]$data$cell_type), c("T", "B"))
  built <- ggplot2::ggplot_build(score_plot)
  expect_true(all(built$data[[1L]]$colour == "grey85"))
  expect_equal(built$plot$scales$get_scales("colour")$get_limits(), c(2, 2))

  none <- plot_seurat_feature(seu, "score", highlight = "1", min_cutoff = 10)
  expect_identical(none$labels$title, "score \u00b7 C1 \u00b7 n = 0 / 2")
  expect_equal(nrow(none$layers[[2L]]$data), 0L)
  withr::with_pdf(tempfile(fileext = ".pdf"), {
    expect_s3_class(ggplot2::ggplotGrob(none), "gtable")
  })
})
