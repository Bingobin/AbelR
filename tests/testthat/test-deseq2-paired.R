test_that("paired DESeq2 fits pair effects and retains the treatment contrast", {
  skip_if_not_installed("DESeq2")
  set.seed(614)
  design <- data.frame(
    Group = rep(c("SH2", "SCR"), 3), Batch = "B1", Library = "RNA",
    Pair = rep(c("dish1", "dish2", "dish3"), each = 2),
    row.names = c("SH2_1", "SCR_1", "SH2_2", "SCR_2", "SH2_3", "SCR_3")
  )
  counts <- data.frame(GID = paste0("gene", 1:300))
  for (id in rownames(design)) counts[[id]] <- rnbinom(300, mu = 100, size = 5)
  annotation_file <- tempfile(fileext = ".txt")
  plot_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(annotation_file, plot_file)), add = TRUE)
  utils::write.table(
    data.frame(GID = counts$GID, Gene_Type = "protein_coding"),
    annotation_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  args <- list(
    count.matrix = counts, design.df = design,
    library = "RNA", batch = "B1", tr = "SH2", ctr = "SCR",
    gene_anno_file = annotation_file
  )
  for (paired in c(FALSE, TRUE)) {
    args$pair_col <- if (paired) "Pair" else NULL
    expect_warning(res <- do.call(DESeq2_DEG_analysis_batch, args), "no Symbol column")
    expect_equal(
      sum(grepl("^PairID_", res$inter)),
      if (paired) 2L else 0L
    )
    expect_true("ContrastGroup_SH2_vs_SCR" %in% res$inter)
    expect_equal(nrow(res$result), 300L)
    expect_equal(rownames(res$design), colnames(res$vsd))
  }
  expect_equal(as.character(res$design$PairID), design$Pair)

  args$pair_col <- "missing"
  expect_error(do.call(DESeq2_DEG_analysis_batch, args), "pair_col")
  args$pair_col <- "Pair"
  args$design.df$Pair[1] <- NA
  expect_error(do.call(DESeq2_DEG_analysis_batch, args), "missing or empty pair IDs")
  args$design.df <- design[-1, ]
  expect_error(do.call(DESeq2_DEG_analysis_batch, args), "Check pair IDs: dish1")
  args$design.df <- design[1:2, ]
  expect_error(do.call(DESeq2_DEG_analysis_batch, args), "at least two complete")
  args$design.df <- design
  args$design.df$Pair[3:4] <- "dish1"
  expect_error(do.call(DESeq2_DEG_analysis_batch, args), "Check pair IDs: dish1")
})
