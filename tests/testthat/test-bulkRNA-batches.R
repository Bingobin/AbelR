test_that("bulk RNA design derives Pair from the original replicate suffix", {
  count_dir <- tempfile("bulk-pair-")
  dir.create(count_dir)
  on.exit(unlink(count_dir, recursive = TRUE), add = TRUE)
  sample_ids <- c("M13_KAP1_sh9_1", "M13_KAP1_SCR_1",
                  "M13_KAP1_sh9_02", "M13_KAP1_SCR_02")
  counts <- data.frame(GID = c("gene1", "gene2"))
  for (id in sample_ids) counts[[id]] <- c(10, 20)
  for (batch in c("B1", "B2")) {
    for (type in c("count", "tpm")) {
      utils::write.table(
        counts, file.path(count_dir, paste0("KAP1_mRNA_", batch, ".merge.", type, ".txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
  bulk <- build_bulkRNA_batches(
    batches = c("B1", "B2"), target = "KAP1", library_name = "mRNA",
    count_dir = count_dir
  )
  expect_identical(bulk$design$Pair, rep(c("1", "1", "02", "02"), 2))
  expect_identical(bulk$design$Group, rep(c("SH9", "SCR"), 4))
  expect_identical(bulk$design$Batch, rep(c("B1", "B2"), each = 4))
  expect_identical(rownames(bulk$design), colnames(bulk$count)[-1])
  expect_identical(rownames(bulk$design), colnames(bulk$tpm)[-1])
})
