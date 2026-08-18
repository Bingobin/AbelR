test_that("EnrichMSigDB validates species and database names", {
  expect_error(
    EnrichMSigDB("TP53", species = "rat"),
    "human.*mouse"
  )
  expect_error(
    EnrichMSigDB("TP53", database = "UNKNOWN"),
    "Unsupported human database"
  )
  expect_error(
    EnrichMSigDB("TP53", database = "MH", species = "human"),
    "Unsupported human database"
  )
  expect_error(
    EnrichMSigDB("Trp53", database = "H", species = "mouse"),
    "Unsupported mouse database"
  )
})

test_that("EnrichMSigDB validates gene and enrichment parameters", {
  result <- "not evaluated"
  expect_warning(
    result <- EnrichMSigDB(character()),
    "No valid genes.*returning NULL"
  )
  expect_null(result)
  expect_error(EnrichMSigDB("TP53", pvalueCutoff = 2), "between 0 and 1")
  expect_error(EnrichMSigDB("TP53", minGSSize = 0), "positive integer")
  expect_error(
    EnrichMSigDB("TP53", minGSSize = 20, maxGSSize = 10),
    "greater than or equal"
  )
})

test_that("MSigDB database maps cover human and native mouse collections", {
  human_map <- .abel_msigdb_database_map("human")
  mouse_map <- .abel_msigdb_database_map("mouse")

  expect_true(all(c(
    "H", "C1", "C2", "C2:CP:REACTOME", "C3:TFT:GTRD",
    "C4", "C5:GO:BP", "C6", "C7", "C8", "C9"
  ) %in% names(human_map)))
  expect_true(all(c(
    "MH", "M1", "M2", "M2:CP:REACTOME", "M3:MIRDB",
    "M3:GTRD", "M5:GO:BP", "M7", "M8"
  ) %in% names(mouse_map)))
  expect_false("H" %in% names(mouse_map))
  expect_false("MH" %in% names(human_map))
})

test_that("MSigDB database names accept colon, underscore, and legacy forms", {
  expect_equal(
    .abel_match_msigdb_database("c2:cp:reactome", "human"),
    "C2:CP:REACTOME"
  )
  expect_equal(
    .abel_match_msigdb_database("C2_CP_REACTOME", "human"),
    "C2:CP:REACTOME"
  )
  expect_equal(
    .abel_match_msigdb_database("M3:MIR:MIRDB", "mouse"),
    "M3:MIRDB"
  )
  expect_true(is.na(.abel_match_msigdb_database("H", "mouse")))
})

make_gsea_result <- function(prefix, p_adjust) {
  data.frame(
    Description = paste0(prefix, seq_along(p_adjust)),
    NES = seq_along(p_adjust) * c(1, -1, 1)[seq_along(p_adjust)],
    pvalue = p_adjust / 2,
    p.adjust = p_adjust,
    stringsAsFactors = FALSE
  )
}

test_that("plot_gsea_dotplot preserves single-result behaviour", {
  result <- make_gsea_result("term", c(0.03, 0.001, 0.02))
  plot <- plot_gsea_dotplot(result, top_n = 2)

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 2L)
  expect_setequal(plot$data$Description, c("term2", "term3"))
})

test_that("plot_gsea_dotplot compares group results on a shared pathway axis", {
  first <- data.frame(
    Description = c("P1", "P2", "P3"),
    NES = c(1, -1, 0.5),
    pvalue = c(0.01, 0.02, 0.03),
    p.adjust = c(0.01, 0.02, 0.03)
  )
  second <- data.frame(
    Description = c("P1", "P2", "P3"),
    NES = c(-1, 1, 0.8),
    pvalue = c(0.03, 0.001, 0.02),
    p.adjust = c(0.03, 0.001, 0.02)
  )
  results <- list(
    Control = first,
    Treatment = second
  )

  plot <- plot_gsea_dotplot(
    results,
    top_n = c(Control = 1, Treatment = 2),
    group_order = c("Treatment", "Control")
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(levels(plot$data$PlotGroup), c("Treatment", "Control"))
  expect_equal(nrow(plot$data), 6L)
  expect_setequal(as.character(plot$data$Description), c("P1", "P2", "P3"))
  expect_equal(rev(levels(plot$data$Description)), c("P2", "P3", "P1"))
})

test_that("plot_gsea_dotplot keeps requested pathway order from top to bottom", {
  results <- list(
    Control = make_gsea_result("P", c(0.03, 0.001, 0.02)),
    Treatment = make_gsea_result("P", c(0.04, 0.002, 0.01))
  )

  plot <- plot_gsea_dotplot(
    results,
    top_n = 1,
    label_pathways = c("P3", "P1")
  )

  displayed_order <- rev(levels(plot$data$Description))
  expect_equal(displayed_order[seq_len(2)], c("P3", "P1"))
  expect_equal(displayed_order, c("P3", "P1", "P2"))
})

test_that("plot_gsea_dotplot labels unnamed group results", {
  results <- unname(list(
    make_gsea_result("P", c(0.03, 0.001, 0.02)),
    make_gsea_result("P", c(0.04, 0.002, 0.01))
  ))

  plot <- NULL
  expect_warning(
    plot <- plot_gsea_dotplot(results, top_n = 1),
    "assigning Group1"
  )
  expect_equal(levels(plot$data$PlotGroup), c("Group1", "Group2"))
})

test_that("plot_gsea_dotplot moves a common MSigDB prefix to the y-axis", {
  hallmark <- data.frame(
    Description = c(
      "HALLMARK_APOPTOSIS",
      "HALLMARK_MYC_TARGETS_V1",
      "HALLMARK_E2F_TARGETS"
    ),
    NES = c(1.2, -1.5, 1.8),
    pvalue = c(0.01, 0.02, 0.03),
    p.adjust = c(0.01, 0.02, 0.03)
  )

  plot <- plot_gsea_dotplot(hallmark, top_n = 3)

  expect_equal(plot$labels$x, "GSEA")
  expect_equal(plot$labels$y, "HALLMARK pathway")
  expect_equal(
    rev(levels(plot$data$DisplayTerm)),
    c("APOPTOSIS", "MYC_TARGETS_V1", "E2F_TARGETS")
  )
  expect_true(all(grepl("^HALLMARK_", plot$data$Description)))
})
