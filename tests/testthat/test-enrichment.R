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
    "H", "C1", "C2", "C2_CP_REACTOME", "C3_TFT_GTRD",
    "C4", "C5_GO_BP", "C6", "C7", "C8", "C9"
  ) %in% names(human_map)))
  expect_true(all(c(
    "MH", "M1", "M2", "M2_CP_REACTOME", "M3_MIRDB",
    "M3_GTRD", "M5_GO_BP", "M7", "M8"
  ) %in% names(mouse_map)))
  expect_false("H" %in% names(mouse_map))
  expect_false("MH" %in% names(human_map))
})

test_that("MSigDB database names accept colon, underscore, and legacy forms", {
  expect_equal(
    .abel_match_msigdb_database("c2:cp:reactome", "human"),
    "C2_CP_REACTOME"
  )
  expect_equal(
    .abel_match_msigdb_database("C2_CP_REACTOME", "human"),
    "C2_CP_REACTOME"
  )
  expect_equal(
    .abel_match_msigdb_database("M3:MIR:MIRDB", "mouse"),
    "M3_MIRDB"
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

test_that("plot_gsea_summary handles a single result", {
  result <- make_gsea_result("term", c(0.03, 0.001, 0.02))
  summary <- plot_gsea_summary(result, top_n = 2)

  expect_named(summary, c("dotplot", "heatmap", "table"))
  expect_s3_class(summary$dotplot, "ggplot")
  expect_s3_class(summary$heatmap, "ggplot")
  expect_equal(nrow(summary$table), 2L)
  expect_setequal(summary$table$Description, c("term2", "term3"))
})

test_that("plot_gsea_summary compares groups on a shared pathway axis", {
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

  summary <- plot_gsea_summary(
    results,
    top_n = c(Control = 1, Treatment = 2),
    group_order = c("Treatment", "Control")
  )

  expect_equal(levels(summary$table$Group), c("Treatment", "Control"))
  expect_equal(nrow(summary$table), 6L)
  expect_setequal(
    as.character(summary$table$Description),
    c("P1", "P2", "P3")
  )
  expect_equal(
    rev(levels(summary$table$Description)),
    c("P2", "P3", "P1")
  )
})

test_that("plot_gsea_summary keeps requested pathway order", {
  results <- list(
    Control = make_gsea_result("P", c(0.03, 0.001, 0.02)),
    Treatment = make_gsea_result("P", c(0.04, 0.002, 0.01))
  )

  summary <- plot_gsea_summary(
    results,
    top_n = 1,
    label_pathways = c("P3", "P1")
  )

  displayed_order <- rev(levels(summary$table$Description))
  expect_equal(displayed_order[seq_len(2)], c("P3", "P1"))
  expect_equal(displayed_order, c("P3", "P1", "P2"))
})

test_that("plot_gsea_summary labels unnamed group results", {
  results <- unname(list(
    make_gsea_result("P", c(0.03, 0.001, 0.02)),
    make_gsea_result("P", c(0.04, 0.002, 0.01))
  ))

  summary <- NULL
  expect_warning(
    summary <- plot_gsea_summary(results, top_n = 1),
    "assigning Group1"
  )
  expect_equal(levels(summary$table$Group), c("Group1", "Group2"))
})

test_that("plot_gsea_summary moves a common MSigDB prefix to the y-axis", {
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

  summary <- plot_gsea_summary(hallmark, top_n = 3)

  expect_equal(summary$dotplot$labels$x, "GSEA")
  expect_equal(summary$dotplot$labels$y, "HALLMARK pathway")
  expect_equal(summary$heatmap$labels$y, "HALLMARK pathway")
  expect_equal(
    rev(levels(summary$table$DisplayTerm)),
    c("APOPTOSIS", "MYC_TARGETS_V1", "E2F_TARGETS")
  )
  expect_true(all(grepl("^HALLMARK_", summary$table$Description)))
})

test_that("plot_gsea_summary formats pathway labels in sentence case", {
  hallmark <- data.frame(
    Description = c(
      "HALLMARK_MYC-TARGETS_V1",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
    ),
    NES = c(1.2, -1.5),
    pvalue = c(0.001, 0.02),
    p.adjust = c(0.002, 0.03)
  )

  summary <- plot_gsea_summary(
    hallmark,
    top_n = 2,
    sentence_case = TRUE
  )

  expect_equal(
    rev(levels(summary$table$DisplayTerm)),
    c("Myc targets v1", "Tnfa signaling via nfkb")
  )
  expect_true(all(grepl("^HALLMARK_", summary$table$Description)))
})

test_that("plot_gsea_summary completes missing groups and adds plot layers", {
  results <- list(
    GroupA = data.frame(
      Description = c("HALLMARK_A", "HALLMARK_B"),
      NES = c(1.5, -1.2),
      pvalue = c(0.001, 0.2),
      p.adjust = c(0.004, 0.3)
    ),
    GroupB = data.frame(
      Description = "HALLMARK_A",
      NES = -1.1,
      pvalue = 0.02,
      p.adjust = 0.04
    )
  )

  summary <- plot_gsea_summary(
    results,
    top_n = 0,
    label_pathways = c("HALLMARK_A", "HALLMARK_B")
  )

  expect_equal(nrow(summary$table), 4L)
  missing_row <- summary$table[
    summary$table$Description == "HALLMARK_B" &
      summary$table$Group == "GroupB",
  ]
  expect_true(is.na(missing_row$NES))
  expect_equal(nrow(missing_row), 1L)
  expect_equal(length(summary$dotplot$layers), 2L)
  expect_equal(length(summary$heatmap$layers), 2L)
  expect_equal(
    as.character(summary$table$Significance[
      summary$table$Description == "HALLMARK_A" &
        summary$table$Group == "GroupA"
    ]),
    "**"
  )
})
