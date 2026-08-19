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
  pathway_id <- paste0(prefix, seq_along(p_adjust))
  data.frame(
    ID = pathway_id,
    Description = pathway_id,
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
    ID = c("P1", "P2", "P3"),
    Description = c("P1", "P2", "P3"),
    NES = c(1, -1, 0.5),
    pvalue = c(0.01, 0.02, 0.03),
    p.adjust = c(0.01, 0.02, 0.03)
  )
  second <- data.frame(
    ID = c("P1", "P2", "P3"),
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
    rev(levels(summary$table$ID)),
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

  displayed_order <- rev(levels(summary$table$ID))
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
    ID = c(
      "HALLMARK_APOPTOSIS",
      "HALLMARK_MYC_TARGETS_V1",
      "HALLMARK_E2F_TARGETS"
    ),
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
  expect_true(all(grepl("^HALLMARK_", summary$table$ID)))
})

test_that("plot_gsea_summary formats pathway labels in sentence case", {
  hallmark <- data.frame(
    ID = c(
      "HALLMARK_MYC-TARGETS_V1",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
    ),
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
  expect_true(all(grepl("^HALLMARK_", summary$table$ID)))
})

test_that("plot_gsea_summary completes missing groups and adds plot layers", {
  results <- list(
    GroupA = data.frame(
      ID = c("HALLMARK_A", "HALLMARK_B"),
      Description = c("HALLMARK_A", "HALLMARK_B"),
      NES = c(1.5, -1.2),
      pvalue = c(0.001, 0.2),
      p.adjust = c(0.004, 0.3)
    ),
    GroupB = data.frame(
      ID = "HALLMARK_A",
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
    summary$table$SelectedPathway == "HALLMARK_B" &
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

test_that("plot_gsea_summary selects P-value and pathway columns", {
  result <- data.frame(
    ID = c("ID_A", "ID_B"),
    Description = c("Description A", "Description B"),
    NES = c(1.5, -1.2),
    pvalue = c(0.001, 0.02),
    p.adjust = c(0.5, 0.03)
  )

  by_adjusted <- plot_gsea_summary(result, top_n = 1)
  by_raw_description <- plot_gsea_summary(
    result,
    p_column = "pvalue",
    pathway_column = "Description",
    top_n = 1
  )

  expect_equal(as.character(by_adjusted$table$ID), "ID_B")
  expect_equal(
    as.character(by_raw_description$table$Description),
    "Description A"
  )
  expect_equal(by_raw_description$table$PlotP, 0.001)
  expect_equal(
    by_raw_description$dotplot$labels$size,
    "-log10(pvalue)"
  )
  expect_equal(by_raw_description$table$Significance, "**")
})

test_that("plot_gsea_summary selects one database from nested groups", {
  make_database_result <- function(database, id, p_adjust, nes) {
    data.frame(
      Database = database,
      ID = id,
      Description = paste(id, "description"),
      NES = nes,
      pvalue = p_adjust / 2,
      p.adjust = p_adjust
    )
  }
  results <- list(
    GroupA = list(
      H = make_database_result("H", "HALLMARK_A", 0.01, 1.2),
      C5_GO_BP = make_database_result("C5_GO_BP", "GOBP_A", 0.02, -1.1)
    ),
    GroupB = list(
      H = make_database_result("H", "HALLMARK_B", 0.03, -1.3),
      C5_GO_BP = make_database_result("C5_GO_BP", "GOBP_B", 0.04, 1.4)
    )
  )

  summary <- plot_gsea_summary(
    results,
    database = "C5:GO:BP",
    top_n = 1
  )

  expect_equal(levels(summary$table$Group), c("GroupA", "GroupB"))
  expect_setequal(as.character(summary$table$ID), c("GOBP_A", "GOBP_B"))
  expect_true(all(
    stats::na.omit(summary$table$Database) == "C5_GO_BP"
  ))
  expect_equal(summary$dotplot$labels$y, "GOBP pathway")
})

make_enrichment_result <- function(database, descriptions, p_adjust, count) {
  data.frame(
    Database = database,
    Species = "human",
    ID = paste0("ID", seq_along(descriptions)),
    Description = descriptions,
    GeneRatio = paste0(count, "/100"),
    BgRatio = "10/1000",
    RichFactor = count / 10,
    FoldEnrichment = count / 2,
    zScore = count / 3,
    pvalue = p_adjust / 2,
    p.adjust = p_adjust,
    qvalue = p_adjust,
    geneID = "A/B",
    Count = count,
    stringsAsFactors = FALSE
  )
}

test_that("plot_enrichment_summary compares pathway counts across groups", {
  results <- list(
    GroupA = make_enrichment_result(
      "C2_CP_KEGG_LEGACY",
      c("Pathway A", "Pathway B"),
      c(0.001, 0.2),
      c(10, 4)
    ),
    GroupB = make_enrichment_result(
      "C2_CP_KEGG_LEGACY",
      c("Pathway A", "Pathway C"),
      c(0.03, 0.002),
      c(6, 12)
    )
  )

  summary <- plot_enrichment_summary(
    results,
    pathway_column = "Description",
    top_n = 1
  )

  expect_named(summary, c("dotplot", "table"))
  expect_s3_class(summary$dotplot, "ggplot")
  expect_equal(nrow(summary$table), 4L)
  expect_equal(
    rev(levels(summary$table$DisplayTerm)),
    c("Pathway A", "Pathway C")
  )
  expect_equal(
    levels(summary$table$Group),
    c("GroupA", "GroupB")
  )
  expect_equal(summary$dotplot$labels$x, "Group")
  expect_equal(summary$dotplot$labels$y, "Pathway")
  expect_equal(summary$dotplot$labels$size, "Count")
})

test_that("plot_enrichment_summary selects one database from nested groups", {
  results <- list(
    GroupA = list(
      H = make_enrichment_result("H", "Hallmark A", 0.01, 5),
      C2_CP_KEGG_LEGACY = make_enrichment_result(
        "C2_CP_KEGG_LEGACY", "KEGG A", 0.02, 8
      )
    ),
    GroupB = list(
      H = make_enrichment_result("H", "Hallmark B", 0.03, 6),
      C2_CP_KEGG_LEGACY = make_enrichment_result(
        "C2_CP_KEGG_LEGACY", "KEGG B", 0.04, 7
      )
    )
  )

  summary <- plot_enrichment_summary(
    results,
    database = "C2:CP:KEGG_LEGACY",
    pathway_column = "Description",
    top_n = 1,
    group_order = c("GroupB", "GroupA"),
    low_color = "white",
    high_color = "darkgreen"
  )

  expect_setequal(
    as.character(summary$table$Description),
    c("KEGG A", "KEGG B")
  )
  expect_true(all(
    stats::na.omit(as.character(summary$table$Database)) ==
      "C2_CP_KEGG_LEGACY"
  ))
  expect_equal(
    levels(summary$table$Group),
    c("GroupB", "GroupA")
  )
  expect_equal(summary$dotplot$scales$scales[[1]]$palette(c(0, 1)),
    c("#FFFFFF", "#006400")
  )
})

test_that("plot_enrichment_summary keeps labelled pathways first", {
  results <- list(
    GroupA = make_enrichment_result(
      "H",
      c("Pathway A", "Pathway B", "Pathway C"),
      c(0.03, 0.001, 0.02),
      c(4, 8, 6)
    ),
    GroupB = make_enrichment_result(
      "H",
      c("Pathway A", "Pathway B", "Pathway C"),
      c(0.04, 0.002, 0.01),
      c(5, 7, 9)
    )
  )

  summary <- plot_enrichment_summary(
    results,
    pathway_column = "Description",
    top_n = 1,
    label_pathways = c("Pathway C", "Pathway A")
  )

  expect_equal(
    rev(levels(summary$table$Description)),
    c("Pathway C", "Pathway A", "Pathway B")
  )
})

test_that("plot_enrichment_summary selects the P-value and pathway columns", {
  result <- make_enrichment_result(
    "H",
    c("Description A", "Description B"),
    c(0.5, 0.03),
    c(5, 8)
  )
  result$pvalue <- c(0.001, 0.02)

  by_adjusted <- plot_enrichment_summary(result, top_n = 1)
  by_raw_id <- plot_enrichment_summary(
    result,
    p_column = "pvalue",
    pathway_column = "ID",
    top_n = 1
  )

  expect_equal(as.character(by_adjusted$table$ID), "ID2")
  expect_equal(as.character(by_raw_id$table$ID), "ID1")
  expect_equal(by_raw_id$table$PlotP, 0.001)
  expect_equal(
    by_raw_id$dotplot$scales$scales[[1]]$name,
    "-log10(pvalue)"
  )
})

test_that("plot_enrichment_summary requires both pathway label columns", {
  result <- make_enrichment_result("H", "Pathway A", 0.01, 5)
  result$ID <- NULL

  expect_warning(
    expect_error(
      plot_enrichment_summary(result),
      "no valid enrichment results"
    ),
    "missing: ID"
  )
})

test_that("plot_enrichment_summary caps plotted negative log10 P values", {
  result <- make_enrichment_result(
    "H",
    c("Pathway A", "Pathway B"),
    c(1e-20, 0.01),
    c(5, 8)
  )

  uncapped <- plot_enrichment_summary(result, top_n = 2)
  capped <- plot_enrichment_summary(result, top_n = 2, max_logp = 10)

  expect_equal(max(uncapped$table$LogP), 20)
  expect_equal(max(capped$table$LogP), 10)
  expect_equal(
    capped$dotplot$scales$scales[[1]]$name,
    "-log10(p.adjust)\n(max 10)"
  )
  expect_error(
    plot_enrichment_summary(result, max_logp = 0),
    "positive finite"
  )
})

test_that("plot_enrichment_summary sentence-cases pathway IDs", {
  result <- make_enrichment_result(
    "C5_GO_BP",
    c("Long pathway definition", "Another definition"),
    c(0.01, 0.02),
    c(5, 8)
  )
  result$ID <- c("GOBP_PIGMENTATION", "GOBP_PROTEIN-FOLDING")

  summary <- plot_enrichment_summary(
    result,
    top_n = 2,
    sentence_case = TRUE
  )

  expect_equal(
    rev(levels(summary$table$DisplayTerm)),
    c("Gobp pigmentation", "Gobp protein folding")
  )
  expect_equal(
    as.character(summary$table$ID),
    c("GOBP_PIGMENTATION", "GOBP_PROTEIN-FOLDING")
  )
})
