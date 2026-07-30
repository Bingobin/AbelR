test_that("UKB ICD-10 records are aligned with diagnosis-date arrays", {
  phenotype <- data.frame(
    participant.eid = c(1, 2),
    diagnoses = c(
      "I26 Pulmonary embolism|I63 Cerebral infarction",
      "N18 Chronic kidney disease"
    ),
    assessment = as.Date(c("2020-01-01", "2020-01-01"))
  )
  diagnosis_dates <- data.frame(
    participant.eid = c(1, 2),
    p41280_a0 = c("2021-01-01", "2022-01-01"),
    p41280_a1 = c("2023-01-01", NA)
  )

  result <- prepare_ukb_icd10_events(
    phenotype,
    diagnosis_dates,
    diagnosis_col = "diagnoses",
    assessment_date_col = "assessment"
  )

  expect_equal(result$ICD10, c("I26", "I63", "N18"))
  expect_equal(result$array, c(0L, 1L, 0L))
  expect_equal(result$interval_days[1], 366)
})


test_that("thrombosis codes use the requested classification groups", {
  result <- classify_ukb_thrombosis(
    c("I26", "I63", "N18", NA),
    include_groups = c("ATE", "VTE")
  )

  expect_equal(
    result$thromb_group,
    c("VTE", "ATE", NA_character_, NA_character_)
  )
})


test_that("CHIP call sets and carrier classes preserve VAF information", {
  inhouse <- data.frame(
    Hugo_Symbol = c("TET2", "DNMT3A"),
    Chromosome = "1",
    Start_Position = 1:2,
    End_Position = 1:2,
    Reference_Allele = "A",
    Tumor_Seq_Allele2 = "T",
    Tumor_Sample_Barcode = c("p1", "p2"),
    Mean_VAF = c(0.12, 0.05)
  )
  released <- data.frame(
    Hugo_Symbol = c("TET2", "JAK2"),
    Chromosome = "1",
    Start_Position = c(1, 3),
    End_Position = c(1, 3),
    Reference_Allele = "A",
    Tumor_Seq_Allele2 = "T",
    Tumor_Sample_Barcode = c("p1", "p3"),
    VAF = c(0.11, 0.25)
  )

  merged <- merge_ukb_chip_callsets(
    inhouse,
    released,
    inhouse_vaf_col = "Mean_VAF"
  )
  carriers <- classify_chip_carriers(
    merged,
    cohort_ids = c("p1", "p2", "p3")
  )

  expect_equal(nrow(merged), 3)
  expect_equal(sum(merged$CHIP_source == "Shared"), 1)
  expect_equal(carriers$CHIP.lar[carriers$eid == "p1"], 1)
  expect_equal(carriers$DNMT3A.sma[carriers$eid == "p2"], 1)
  expect_equal(carriers$JAK2.lar[carriers$eid == "p3"], 1)
})


test_that("competing-risk data exclude prevalent events and code death", {
  phenotype <- data.frame(
    eid = c("p1", "p2", "p3"),
    Date = as.Date("2020-01-01"),
    Date_Death = as.Date(c(NA, "2022-01-01", NA))
  )
  icd10 <- data.frame(
    eid = c("p1", "p3"),
    ICD10 = c("I26", "I26"),
    Diagnosis = "Pulmonary embolism",
    date = as.Date(c("2021-01-01", "2020-03-01"))
  )

  result <- make_ukb_competing_risk_data(
    icd10,
    phenotype,
    icd_pattern = "^I26",
    lag_days = 180,
    censor_date = as.Date("2023-01-01")
  )

  expect_setequal(result$eid, c("p1", "p2"))
  expect_equal(result$status_num[result$eid == "p1"], 1)
  expect_equal(result$status_num[result$eid == "p2"], 2)
})


test_that("Fine-Gray analysis returns model and cumulative-incidence outputs", {
  skip_if_not_installed("cmprsk")
  set.seed(12)
  n <- 80
  data <- data.frame(
    eid = paste0("p", seq_len(n)),
    time_years = stats::rexp(n, rate = 0.15),
    status_num = sample(
      c(0, 1, 2),
      n,
      replace = TRUE,
      prob = c(0.55, 0.30, 0.15)
    ),
    CHIP = factor(rep(c("CHIP-", "CHIP+"), n / 2)),
    Age = stats::rnorm(n, 60, 6)
  )

  result <- fit_ukb_competing_risk(
    data,
    exposure = "CHIP",
    covariates = "Age"
  )

  expect_s3_class(result$fit_crude, "coxph")
  expect_s3_class(result$fit_adjusted, "coxph")
  expect_true(nrow(result$estimates) > 0)
  expect_s3_class(result$cif_plot, "ggplot")
  expect_s3_class(result$forest_plot, "ggplot")
})


test_that("UKB omics limma returns one result per variable trait", {
  skip_if_not_installed("limma")
  set.seed(10)
  data <- data.frame(
    eid = paste0("p", seq_len(60)),
    CHIP.all = rep(0:1, each = 30),
    Age = stats::rnorm(60, 60, 5),
    Sex = rep(rep(0:1, each = 15), 2),
    feature_a = stats::rnorm(60),
    feature_b = stats::rnorm(60)
  )

  result <- run_ukb_omics_limma(
    data,
    features = c("feature_a", "feature_b"),
    covariates = c("Age", "Sex")
  )

  expect_equal(nrow(result$result), 2)
  expect_equal(result$coefficient, "CHIP.all")
  expect_setequal(result$result$Feature, c("feature_a", "feature_b"))
})


test_that("multi-omics scoring returns module and combined scores", {
  set.seed(11)
  data <- data.frame(
    feature_a = stats::rnorm(30),
    feature_b = stats::rnorm(30),
    feature_c = stats::rnorm(30)
  )
  feature_sets <- list(
    Olink = c("feature_a", "feature_b"),
    CBC = "feature_c"
  )
  effects <- list(
    Olink = c(feature_a = 0.2, feature_b = -0.3),
    CBC = c(feature_c = 0.4)
  )

  result <- build_ukb_multiomic_score(
    data,
    feature_sets = feature_sets,
    effect_tables = effects
  )

  expect_equal(nrow(result$scores), 30)
  expect_true(
    all(
      c("score_olink", "score_cbc", "score_multiomic_z") %in%
        colnames(result$data)
    )
  )
  expect_equal(
    unname(stats::sd(result$data$score_multiomic_z)),
    1,
    tolerance = 1e-10
  )
})
