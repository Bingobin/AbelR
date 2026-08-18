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
