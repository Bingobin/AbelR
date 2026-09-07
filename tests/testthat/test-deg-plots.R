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


test_that("DESeq2 volcano replaces zeros using the full selected P-value column", {
  skip_if_not_installed("ggrepel")
  result <- data.frame(
    Symbol = c("ZERO", "NONZERO", "UNPLOTTED_MIN", "MISSING"),
    log2FoldChange = c(2, -2, 0, 0),
    pvalue = c(0, 1e-5, 1e-15, NA),
    padj = c(0, 1e-3, 1e-10, NA)
  )
  original <- result
  for (adjust in c(TRUE, FALSE)) {
    plot <- volcano_plot_Deseq2(
      result, gene.list = character(), n = 0, top = 1,
      adjust = adjust, max_y = 100
    )
    expect_false("UNPLOTTED_MIN" %in% plot$data$Symbol)
    expect_equal(
      plot$data$padj[plot$data$Symbol == "ZERO"],
      if (adjust) 1e-10 else 1e-15
    )
    expect_equal(plot$data$group[plot$data$Symbol == "ZERO"], "up")
    expect_equal(plot$data$label[plot$data$Symbol == "ZERO"], "ZERO")
  }
  expect_identical(result, original)
})


test_that("DESeq2 volcano warns when all available P values are zero", {
  skip_if_not_installed("ggrepel")
  result <- data.frame(
    Symbol = c("UP", "DOWN"), log2FoldChange = c(2, -2),
    pvalue = c(0, 0), padj = c(0, 0)
  )
  expect_warning(
    plot <- volcano_plot_Deseq2(
      result, gene.list = character(), n = 0, top = 0,
      adjust = TRUE, max_y = 20
    ),
    "No finite positive P value.*padj.*max_y"
  )
  expect_equal(-log10(plot$data$padj), c(20, 20))
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


test_that("DEG comparison accepts independent axis columns and settings", {
  skip_if_not_installed("ggrepel")

  x_deg_result <- data.frame(
    GID = paste0("ENSG", 1:6),
    Symbol = paste0("GENE", 1:6),
    statistic = c(3, -3, 3, -3, 0.2, 0.4),
    pvalue = c(rep(0.005, 4), 0.5, 0.6),
    FDR = c(rep(0.01, 4), 0.7, 0.8)
  )
  y_deg_result <- data.frame(
    gene_id = paste0("ENSG", 1:6),
    score = c(4, -4, -4, 4, 0.1, 0.2),
    pvalue = c(rep(0.004, 4), 0.4, 0.5),
    qvalue = c(rep(0.02, 4), 0.6, 0.7)
  )

  plot <- plot_deg_comparison(
    x_deg_result = x_deg_result,
    y_deg_result = y_deg_result,
    x_label = "X statistic",
    y_label = "Y score",
    pc = FALSE,
    bg_num = 2,
    show_cor = FALSE,
    top = 0,
    x_col = "statistic",
    y_col = "score",
    x_pvalue_col = "FDR",
    y_pvalue_col = "qvalue",
    pvalue_cutoff = 0.05,
    x_threshold = 2,
    y_threshold = 3,
    x_limits = c(-2.5, 2.5),
    y_limits = c(-5, 5),
    x_merge_col = "GID",
    y_merge_col = "gene_id"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$x, "X statistic")
  expect_equal(plot$labels$y, "Y score")
  expect_equal(plot$coordinates$limits$x, c(-2.5, 2.5))
  expect_equal(plot$coordinates$limits$y, c(-5, 5))
  expect_equal(
    plot$data$Group[match(paste0("GENE", 1:4), plot$data$Symbol)],
    c("upup", "dwdw", "updw", "dwup")
  )
  expect_equal(range(plot$data$Plot_X), c(-2.5, 2.5))
  expect_equal(range(plot$data$Plot_Y), c(-4, 4))
  expect_equal(plot$data$Merge_ID, paste0("ENSG", 1:6))
})


test_that("DEG comparison always labels found goi and warns for missing goi", {
  skip_if_not_installed("ggrepel")

  x_deg_result <- data.frame(
    GID = paste0("ENSG", 1:6),
    Symbol = paste0("GENE", 1:6),
    statistic = c(3, -3, 3, -3, 0.2, 0.4),
    FDR = c(rep(0.01, 4), 0.7, 0.8)
  )
  y_deg_result <- data.frame(
    gene_id = paste0("ENSG", 1:6),
    score = c(4, -4, -4, 4, 0.1, 0.2),
    qvalue = c(rep(0.02, 4), 0.6, 0.7)
  )

  expect_warning(
    plot <- plot_deg_comparison(
      x_deg_result = x_deg_result,
      y_deg_result = y_deg_result,
      pc = FALSE,
      bg_num = 0,
      show_cor = FALSE,
      top = 0,
      goi = c("GENE5", "NOT_FOUND"),
      x_col = "statistic",
      y_col = "score",
      x_pvalue_col = "FDR",
      y_pvalue_col = "qvalue",
      x_threshold = 2,
      y_threshold = 3,
      x_merge_col = "GID",
      y_merge_col = "gene_id"
    ),
    "NOT_FOUND"
  )

  expect_equal(plot$data$GOI[plot$data$Symbol == "GENE5"], "GENE5")
  expect_equal(plot$data$Group[plot$data$Symbol == "GENE5"], "other")
  expect_false("GENE6" %in% plot$data$Symbol)
  label_layer <- plot$layers[[4]]
  expect_equal(label_layer$data$Symbol, "GENE5")
})


test_that("DEG comparison annotates missing Gene_Type by species", {
  skip_if_not_installed("ggrepel")

  gene_ids <- c(
    "ENSG00000186092",
    "ENSG00000284733",
    "ENSG00000284662",
    "ENSG00000187634",
    "ENSG00000290825"
  )
  symbols <- c("OR4F5", "OR4F29", "OR4F16", "SAMD11", "DDX11L2")
  x_deg_result <- data.frame(
    GID = gene_ids,
    Symbol = symbols,
    log2FoldChange = c(2, -2, 1, -1, 3),
    padj = rep(0.01, 5)
  )
  y_deg_result <- data.frame(
    gene_id = gene_ids,
    log2FoldChange = c(1, -1, 2, -2, 3),
    padj = rep(0.01, 5)
  )

  plot <- plot_deg_comparison(
    x_deg_result = x_deg_result,
    y_deg_result = y_deg_result,
    x_merge_col = "GID",
    y_merge_col = "gene_id",
    species = "human",
    bg_num = 10,
    show_cor = FALSE,
    top = 0
  )

  expect_s3_class(plot, "ggplot")
  expect_setequal(plot$data$Symbol, symbols[1:4])
  expect_true(all(plot$data$Gene_Type == "protein_coding"))
})


test_that("DEG comparison validates selected axis columns", {
  skip_if_not_installed("ggrepel")

  x_deg_result <- data.frame(
    Row.names = paste0("ENSG", 1:3),
    Symbol = paste0("GENE", 1:3),
    log2FoldChange = c(-1, 0, 1),
    pvalue = c(0.1, 0.2, 0.3),
    padj = c(0.2, 0.3, 0.4)
  )
  y_deg_result <- data.frame(
    Row.names = paste0("ENSG", 1:3),
    log2FoldChange = c(-2, 0, 2),
    pvalue = c(0.1, 0.2, 0.3),
    padj = c(0.2, 0.3, 0.4)
  )

  expect_error(
    plot_deg_comparison(
      x_deg_result,
      y_deg_result,
      pc = FALSE,
      x_col = "missing_column"
    ),
    "x_col.*missing_column.*x_deg_result"
  )
  expect_error(
    plot_deg_comparison(
      x_deg_result,
      y_deg_result,
      pc = FALSE,
      x_limits = c(1, -1)
    ),
    "x_limits must be ordered"
  )
  expect_error(
    plot_deg_comparison(
      x_deg_result,
      y_deg_result,
      pc = FALSE,
      x_merge_col = "missing_id"
    ),
    "x_merge_col.*missing_id.*x_deg_result"
  )

  duplicated_x_result <- rbind(
    x_deg_result,
    x_deg_result[1, , drop = FALSE]
  )
  expect_error(
    plot_deg_comparison(
      duplicated_x_result,
      y_deg_result,
      pc = FALSE
    ),
    "x_merge_col.*duplicated"
  )
  expect_error(
    plot_deg_comparison(
      x_deg_result,
      y_deg_result,
      pc = FALSE,
      x_pvalue_col = "missing_pvalue"
    ),
    "x_pvalue_col.*missing_pvalue.*x_deg_result"
  )
})


test_that("DEG Manhattan labels top up and down genes separately", {
  skip_if_not_installed("ggrepel")

  annotation <- data.frame(
    GID = paste0("ENSG", 1:7),
    Symbol = c(paste0("UP", 1:3), paste0("DOWN", 1:3), "NS"),
    Chr = "chr1",
    Start = seq(100, 700, by = 100),
    End = seq(150, 750, by = 100),
    Gene_Type = "protein_coding"
  )
  annotation_file <- tempfile(fileext = ".txt")
  utils::write.table(
    annotation,
    annotation_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  deg <- data.frame(
    Symbol = annotation$Symbol,
    log2FoldChange = c(1, 2, 3, -1, -2, -3, 0),
    pvalue = c(1e-4, 1e-5, 1e-6, 1e-4, 1e-5, 1e-6, 0.5),
    padj = c(1e-3, 1e-4, 1e-5, 1e-3, 1e-4, 1e-5, 0.8)
  )

  result <- plot_deg_manhattan(
    deg_list = list(Treatment = deg),
    color_map = c(Treatment = "#333333"),
    species = "human",
    gene_anno_file = annotation_file,
    chromosome_lengths = c(chr1 = 1000),
    chr_keep = "chr1",
    top_n = 2
  )

  expect_equal(nrow(result$top_genes), 4)
  regulation_counts <- table(as.character(result$top_genes$Regulation))
  expect_equal(
    as.integer(regulation_counts[c("Down", "Up")]),
    c(2L, 2L)
  )
  expect_setequal(
    result$top_genes$Symbol,
    c("UP2", "UP3", "DOWN2", "DOWN3")
  )
})
