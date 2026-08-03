test_that("DESeq2 volcano labels break P-value ties by fold change", {
  skip_if_not_installed("ggrepel")

  result <- data.frame(
    Symbol = c(
      paste0("UP", 1:4),
      paste0("DOWN", 1:4),
      "BACKGROUND"
    ),
    log2FoldChange = c(1, 2, 3, 4, -1, -2, -3, -4, 0),
    pvalue = c(rep(0, 8), 0.8),
    padj = c(rep(0, 8), 0.9),
    stringsAsFactors = FALSE
  )

  plot <- volcano_plot_Deseq2(
    result,
    gene.list = character(),
    n = 0,
    pv = 0.05,
    fc = 1.5,
    top = 2,
    adjust = TRUE
  )

  labelled <- plot$data$Symbol[plot$data$label != ""]
  expect_setequal(labelled, c("UP3", "UP4", "DOWN3", "DOWN4"))
  expect_equal(length(labelled), 4)
})


test_that("DESeq2 volcano prioritizes P value before fold change", {
  skip_if_not_installed("ggrepel")

  result <- data.frame(
    Symbol = c("LOW_P", "HIGH_FC", "DOWN_LOW_P", "DOWN_HIGH_FC"),
    log2FoldChange = c(1, 10, -1, -10),
    pvalue = c(1e-20, 1e-10, 1e-20, 1e-10),
    padj = c(1e-18, 1e-8, 1e-18, 1e-8),
    stringsAsFactors = FALSE
  )

  plot <- volcano_plot_Deseq2(
    result,
    gene.list = character(),
    n = 0,
    pv = 0.05,
    fc = 1.5,
    top = 1,
    adjust = TRUE
  )

  labelled <- plot$data$Symbol[plot$data$label != ""]
  expect_setequal(labelled, c("LOW_P", "DOWN_LOW_P"))
})


test_that("DESeq2 volcano accepts custom up and down colors", {
  skip_if_not_installed("ggrepel")

  result <- data.frame(
    Symbol = c("UP", "DOWN", "BACKGROUND"),
    log2FoldChange = c(2, -2, 0),
    pvalue = c(1e-5, 1e-5, 0.8),
    padj = c(1e-4, 1e-4, 0.9),
    stringsAsFactors = FALSE
  )

  plot <- volcano_plot_Deseq2(
    result,
    gene.list = character(),
    n = 1,
    top = 1,
    adjust = TRUE,
    up_color = "#D73027",
    down_color = "#4575B4"
  )

  color_scale <- plot$scales$get_scales("colour")
  expect_equal(
    unname(color_scale$palette(4)),
    c("grey", "#D73027", "#4575B4", "black")
  )
})
