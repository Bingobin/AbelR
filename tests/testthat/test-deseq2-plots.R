test_that("gene TPM plot uses build_bulkRNA_batches-style input", {
  design <- data.frame(
    CellLine = c("M13", "M13", "M13", "M13"),
    Batch = c("B1", "B1", "B2", "B2"),
    Group = c("SCR", "SH9", "SCR", "SH9"),
    Library = c("mRNA", "mRNA", "mRNA", "mRNA"),
    row.names = paste0("sample", 1:4),
    stringsAsFactors = FALSE
  )
  tpm <- data.frame(
    GID = c("ENSG000001.1", "ENSG000002.3"),
    Symbol = c("GENE1", "GENE2"),
    sample1 = c(0, 1),
    sample2 = c(3, 2),
    sample3 = c(7, 3),
    sample4 = c(15, 4),
    check.names = FALSE
  )
  bulk <- list(tpm = tpm, design = design)

  plot <- plot_gene_tpm_median(
    bulk,
    gene = "GENE1",
    group_colors = c(SCR = "grey70", SH9 = "firebrick")
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 4)
  expect_equal(levels(plot$data$Group), c("SCR", "SH9"))
  expect_equal(
    sort(plot$data$Median),
    sort(log2(c(0, 3, 7, 15) + 1))
  )
})


test_that("gene TPM plot matches versionless Ensembl IDs", {
  design <- data.frame(
    CellLine = "M13",
    Batch = "B1",
    Group = "SCR",
    Library = "mRNA",
    row.names = "sample1"
  )
  tpm <- data.frame(
    GID = "ENSG000001.12",
    sample1 = 4,
    check.names = FALSE
  )

  plot <- plot_gene_tpm_median(
    list(tpm = tpm, design = design),
    gene = "ENSG000001",
    group_colors = "steelblue"
  )

  expect_equal(plot$data$Median, log2(5))
})


test_that("gene TPM plot rejects design samples missing from TPM", {
  design <- data.frame(
    CellLine = "M13",
    Batch = "B1",
    Group = "SCR",
    Library = "mRNA",
    row.names = "missing_sample"
  )
  tpm <- data.frame(GID = "ENSG000001", sample1 = 1)

  expect_error(
    plot_gene_tpm_median(
      list(tpm = tpm, design = design),
      gene = "ENSG000001"
    ),
    "Design samples missing from the TPM table"
  )
})
