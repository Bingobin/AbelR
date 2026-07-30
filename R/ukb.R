# Functions for UK Biobank phenotype, CHIP, survival, and multi-omics analyses.


.ukb_assert_columns <- function(data, columns, data_name = "data") {
  missing_columns <- setdiff(columns, colnames(data))
  if (length(missing_columns) > 0L) {
    stop(
      "The following columns are missing from ",
      data_name,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  invisible(TRUE)
}


.ukb_default_thrombosis_codes <- function() {
  data.frame(
    group = c(
      "VTE",
      "VTE",
      "VTE",
      "VTE (pregnancy)",
      "ATE",
      "ATE",
      "ATE",
      "ATE",
      "VTE",
      "Other embolism"
    ),
    subgroup = c(
      "Pulmonary embolism",
      "DVT/phlebitis",
      "Other venous thrombosis",
      "Pregnancy/puerperium VTE",
      "Myocardial infarction/acute coronary thrombosis",
      "Cerebral infarction",
      "Arterial embolism/thrombosis",
      "Visceral/mesenteric infarction",
      "Cerebral venous thrombosis",
      "Fat/air/amniotic embolism"
    ),
    regex = c(
      "^I26",
      "^I80",
      "^(I81|I82)",
      "^O(22\\.3|22\\.5|87\\.1|87\\.3|88\\.2)$",
      "^I2(1|2)|^I24\\.0",
      "^I63",
      "^I74",
      "^(K55\\.0|K76\\.3|N28\\.0|D73\\.5)$",
      "^(I67\\.6|G08)$",
      "^(T79\\.(0|1)|O88\\.(0|1))$"
    ),
    stringsAsFactors = FALSE
  )
}


.ukb_as_date <- function(x, format = NULL) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (is.null(format)) {
    as.Date(as.character(x))
  } else {
    as.Date(as.character(x), format = format)
  }
}


#' Prepare UK Biobank ICD-10 event records
#'
#' Converts a pipe-separated UK Biobank diagnosis field and its corresponding
#' wide diagnosis-date fields into one row per participant and diagnosis. The
#' diagnosis array index is used to align each ICD-10 record with its date.
#'
#' @param phenotype_df Data frame containing participant identifiers, the
#'   pipe-separated ICD-10 diagnosis field, and the assessment date.
#' @param diagnosis_date_df Data frame containing participant identifiers and
#'   one diagnosis-date column per array index.
#' @param diagnosis_col Name of the pipe-separated diagnosis column.
#' @param assessment_date_col Name of the baseline assessment-date column.
#' @param id_col Participant identifier column in `phenotype_df`.
#' @param date_id_col Participant identifier column in `diagnosis_date_df`.
#' @param date_cols Optional character vector of diagnosis-date columns. If
#'   `NULL`, columns ending in an integer array index are used.
#' @param diagnosis_sep Regular expression separating diagnosis entries.
#' @param date_format Optional format passed to [base::as.Date()] for diagnosis
#'   and assessment dates.
#' @param date_join Whether to retain only diagnoses with matched dates
#'   (`"inner"`) or all diagnosis records (`"left"`).
#'
#' @return A data frame with standardized columns `eid`, `ICD10`, `Diagnosis`,
#'   `array`, `date`, `assess_date`, and `interval_days`.
#' @export
prepare_ukb_icd10_events <- function(
  phenotype_df,
  diagnosis_date_df,
  diagnosis_col,
  assessment_date_col,
  id_col = "participant.eid",
  date_id_col = id_col,
  date_cols = NULL,
  diagnosis_sep = "\\|",
  date_format = NULL,
  date_join = c("inner", "left")
) {
  date_join <- match.arg(date_join)
  if (!is.data.frame(phenotype_df) || !is.data.frame(diagnosis_date_df)) {
    stop("phenotype_df and diagnosis_date_df must be data frames.")
  }
  .ukb_assert_columns(
    phenotype_df,
    c(id_col, diagnosis_col, assessment_date_col),
    "phenotype_df"
  )
  .ukb_assert_columns(diagnosis_date_df, date_id_col, "diagnosis_date_df")

  if (anyDuplicated(phenotype_df[[id_col]])) {
    stop("phenotype_df must contain one row per participant.")
  }
  if (anyDuplicated(diagnosis_date_df[[date_id_col]])) {
    stop("diagnosis_date_df must contain one row per participant.")
  }

  if (is.null(date_cols)) {
    candidates <- setdiff(colnames(diagnosis_date_df), date_id_col)
    date_cols <- candidates[grepl("[0-9]+$", candidates)]
  }
  if (length(date_cols) == 0L) {
    stop("No diagnosis-date columns were selected.")
  }
  .ukb_assert_columns(diagnosis_date_df, date_cols, "diagnosis_date_df")

  date_index_text <- regmatches(date_cols, regexpr("[0-9]+$", date_cols))
  date_indices <- suppressWarnings(as.integer(date_index_text))
  if (anyNA(date_indices)) {
    stop("Every diagnosis-date column must end in an integer array index.")
  }
  if (anyDuplicated(date_indices)) {
    stop("Diagnosis-date columns contain duplicated array indices.")
  }

  diagnosis_rows <- lapply(seq_len(nrow(phenotype_df)), function(i) {
    raw_value <- phenotype_df[[diagnosis_col]][i]
    if (is.na(raw_value) || !nzchar(trimws(as.character(raw_value)))) {
      return(NULL)
    }
    records <- trimws(
      strsplit(
        as.character(raw_value),
        split = diagnosis_sep,
        perl = TRUE
      )[[1]]
    )
    records <- records[
      nzchar(records) &
        !toupper(records) %in% c("NA", "N/A")
    ]
    if (length(records) == 0L) {
      return(NULL)
    }
    data.frame(
      .ukb_id = rep(as.character(phenotype_df[[id_col]][i]), length(records)),
      array = seq_along(records) - 1L,
      diagnosis_record = records,
      stringsAsFactors = FALSE
    )
  })
  diagnosis_long <- do.call(rbind, diagnosis_rows)
  if (is.null(diagnosis_long) || nrow(diagnosis_long) == 0L) {
    stop("No non-missing ICD-10 diagnosis records were found.")
  }

  diagnosis_long$ICD10 <- toupper(
    sub("[[:space:]].*$", "", diagnosis_long$diagnosis_record)
  )
  diagnosis_long$Diagnosis <- trimws(
    sub("^[^[:space:]]+[[:space:]]*", "", diagnosis_long$diagnosis_record)
  )
  no_description <- diagnosis_long$Diagnosis == diagnosis_long$diagnosis_record
  diagnosis_long$Diagnosis[no_description] <- NA_character_
  diagnosis_long$diagnosis_record <- NULL

  date_rows <- lapply(seq_along(date_cols), function(i) {
    values <- diagnosis_date_df[[date_cols[i]]]
    keep <- !is.na(values) & nzchar(trimws(as.character(values)))
    if (!any(keep)) {
      return(NULL)
    }
    data.frame(
      .ukb_id = as.character(diagnosis_date_df[[date_id_col]][keep]),
      array = rep(date_indices[i], sum(keep)),
      date = .ukb_as_date(values[keep], format = date_format),
      stringsAsFactors = FALSE
    )
  })
  date_long <- do.call(rbind, date_rows)
  if (is.null(date_long) || nrow(date_long) == 0L) {
    stop("No non-missing diagnosis dates were found.")
  }
  if (anyNA(date_long$date)) {
    stop("Some diagnosis dates could not be parsed as Date values.")
  }
  if (anyDuplicated(date_long[c(".ukb_id", "array")])) {
    stop("Diagnosis dates are not unique for participant and array index.")
  }

  events <- merge(
    diagnosis_long,
    date_long,
    by = c(".ukb_id", "array"),
    all.x = identical(date_join, "left"),
    all.y = FALSE,
    sort = FALSE
  )

  assessment <- data.frame(
    .ukb_id = as.character(phenotype_df[[id_col]]),
    assess_date = .ukb_as_date(
      phenotype_df[[assessment_date_col]],
      format = date_format
    ),
    stringsAsFactors = FALSE
  )
  if (anyNA(assessment$assess_date)) {
    stop("Some assessment dates could not be parsed as Date values.")
  }
  events <- merge(events, assessment, by = ".ukb_id", all.x = TRUE, sort = FALSE)
  events$interval_days <- as.numeric(events$date - events$assess_date)
  events$eid <- events$.ukb_id
  events$.ukb_id <- NULL
  events <- events[
    order(events$eid, events$array),
    c(
      "eid",
      "ICD10",
      "Diagnosis",
      "array",
      "date",
      "assess_date",
      "interval_days"
    ),
    drop = FALSE
  ]
  rownames(events) <- NULL
  events
}


#' Classify thrombosis-related ICD-10 codes
#'
#' Maps ICD-10 codes to arterial thromboembolism, venous thromboembolism, and
#' related subgroups using an ordered regular-expression table.
#'
#' @param icd10 Character vector of ICD-10 codes.
#' @param classification_table Optional data frame containing `group`,
#'   `subgroup`, and `regex` columns. If `NULL`, AbelR's default thrombosis
#'   classification is used.
#' @param include_groups Optional character vector of groups to retain. Other
#'   matches are returned as missing.
#'
#' @return A data frame with `thromb_group`, `thromb_subgroup`, and
#'   `matched_regex` columns in the same order as `icd10`.
#' @export
classify_ukb_thrombosis <- function(
  icd10,
  classification_table = NULL,
  include_groups = NULL
) {
  if (!is.atomic(icd10)) {
    stop("icd10 must be an atomic vector.")
  }
  if (is.null(classification_table)) {
    classification_table <- .ukb_default_thrombosis_codes()
  }
  if (!is.data.frame(classification_table)) {
    stop("classification_table must be a data frame.")
  }
  .ukb_assert_columns(
    classification_table,
    c("group", "subgroup", "regex"),
    "classification_table"
  )
  if (nrow(classification_table) == 0L) {
    stop("classification_table must contain at least one rule.")
  }

  codes <- toupper(trimws(as.character(icd10)))
  codes[is.na(icd10) | !nzchar(codes)] <- NA_character_
  thromb_group <- rep(NA_character_, length(codes))
  thromb_subgroup <- rep(NA_character_, length(codes))
  matched_regex <- rep(NA_character_, length(codes))

  for (i in seq_len(nrow(classification_table))) {
    hit <- !is.na(codes) &
      is.na(thromb_group) &
      grepl(classification_table$regex[i], codes, perl = TRUE)
    thromb_group[hit] <- as.character(classification_table$group[i])
    thromb_subgroup[hit] <- as.character(classification_table$subgroup[i])
    matched_regex[hit] <- as.character(classification_table$regex[i])
  }

  if (!is.null(include_groups)) {
    remove <- !is.na(thromb_group) & !thromb_group %in% include_groups
    thromb_group[remove] <- NA_character_
    thromb_subgroup[remove] <- NA_character_
    matched_regex[remove] <- NA_character_
  }

  data.frame(
    thromb_group = thromb_group,
    thromb_subgroup = thromb_subgroup,
    matched_regex = matched_regex,
    stringsAsFactors = FALSE
  )
}


#' Merge two UKB CHIP variant call sets
#'
#' Performs a full join of an in-house and a UKB-released CHIP call set using
#' explicit variant keys. Both source-specific VAF values are retained, a
#' preferred combined VAF is calculated, and each variant is labelled by call
#' source.
#'
#' @param inhouse_maf In-house CHIP MAF-like data frame.
#' @param ukb_maf UKB-released CHIP MAF-like data frame.
#' @param by Character vector of columns that uniquely identify a variant.
#' @param inhouse_vaf_col VAF column in `inhouse_maf`.
#' @param ukb_vaf_col VAF column in `ukb_maf`.
#' @param prefer Which VAF to use when both call sets contain a variant:
#'   `"inhouse"`, `"ukb"`, or the `"max"` value.
#'
#' @return A merged data frame containing `VAF_inhouse`, `VAF_ukb`, combined
#'   `VAF`, and `CHIP_source`.
#' @export
merge_ukb_chip_callsets <- function(
  inhouse_maf,
  ukb_maf,
  by = c(
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode"
  ),
  inhouse_vaf_col = "VAF",
  ukb_vaf_col = "VAF",
  prefer = c("inhouse", "ukb", "max")
) {
  prefer <- match.arg(prefer)
  if (!is.data.frame(inhouse_maf) || !is.data.frame(ukb_maf)) {
    stop("inhouse_maf and ukb_maf must be data frames.")
  }
  if (!is.character(by) || length(by) == 0L || anyDuplicated(by)) {
    stop("by must be a non-empty vector of unique column names.")
  }
  .ukb_assert_columns(
    inhouse_maf,
    c(by, inhouse_vaf_col),
    "inhouse_maf"
  )
  .ukb_assert_columns(ukb_maf, c(by, ukb_vaf_col), "ukb_maf")
  if (anyDuplicated(inhouse_maf[by])) {
    stop("inhouse_maf contains duplicated variant keys.")
  }
  if (anyDuplicated(ukb_maf[by])) {
    stop("ukb_maf contains duplicated variant keys.")
  }

  inhouse <- inhouse_maf
  ukb <- ukb_maf
  colnames(inhouse)[match(inhouse_vaf_col, colnames(inhouse))] <-
    "VAF_inhouse"
  colnames(ukb)[match(ukb_vaf_col, colnames(ukb))] <- "VAF_ukb"

  merged <- dplyr::full_join(
    inhouse,
    ukb,
    by = by,
    suffix = c(".inhouse", ".ukb")
  )
  merged$VAF_inhouse <- suppressWarnings(as.numeric(merged$VAF_inhouse))
  merged$VAF_ukb <- suppressWarnings(as.numeric(merged$VAF_ukb))

  merged$VAF <- switch(
    prefer,
    inhouse = ifelse(
      is.na(merged$VAF_inhouse),
      merged$VAF_ukb,
      merged$VAF_inhouse
    ),
    ukb = ifelse(
      is.na(merged$VAF_ukb),
      merged$VAF_inhouse,
      merged$VAF_ukb
    ),
    max = {
      result <- pmax(merged$VAF_inhouse, merged$VAF_ukb, na.rm = TRUE)
      result[is.infinite(result)] <- NA_real_
      result
    }
  )
  merged$CHIP_source <- ifelse(
    !is.na(merged$VAF_inhouse) & !is.na(merged$VAF_ukb),
    "Shared",
    ifelse(
      !is.na(merged$VAF_inhouse),
      "In-house only",
      ifelse(!is.na(merged$VAF_ukb), "UKB only", NA_character_)
    )
  )
  merged
}


#' Classify CHIP carriers by gene and clone size
#'
#' Converts a participant-level mutation table into binary carrier indicators
#' for overall CHIP and selected driver genes. Each category contains `.all`,
#' `.sma`, and `.lar` columns for any, small, and large clones.
#'
#' @param maf CHIP mutation data frame.
#' @param cohort_ids Vector of all cohort participant identifiers, including
#'   participants without detected CHIP mutations.
#' @param sample_col Participant identifier column in `maf`.
#' @param gene_col Gene-symbol column in `maf`.
#' @param vaf_col Variant allele fraction column in `maf`.
#' @param genes Driver genes for gene-specific indicators.
#' @param large_clone_cutoff Overall CHIP large-clone VAF cutoff.
#' @param gene_large_cutoffs Optional named numeric vector overriding the
#'   large-clone cutoff for selected genes.
#' @param id_col Name of the participant identifier column in the output.
#' @param na_vaf How mutations with missing VAF are handled: classify carriers
#'   as small (`"small"`), exclude those mutations (`"exclude"`), or stop with
#'   an error (`"error"`).
#'
#' @return A participant-level data frame of binary CHIP indicators.
#' @export
classify_chip_carriers <- function(
  maf,
  cohort_ids,
  sample_col = "Tumor_Sample_Barcode",
  gene_col = "Hugo_Symbol",
  vaf_col = "VAF",
  genes = c(
    "DNMT3A",
    "TET2",
    "ASXL1",
    "PPM1D",
    "TP53",
    "SRSF2",
    "JAK2",
    "SF3B1"
  ),
  large_clone_cutoff = 0.1,
  gene_large_cutoffs = c(JAK2 = 0.2),
  id_col = "eid",
  na_vaf = c("small", "exclude", "error")
) {
  na_vaf <- match.arg(na_vaf)
  if (!is.data.frame(maf)) {
    stop("maf must be a data frame.")
  }
  .ukb_assert_columns(maf, c(sample_col, gene_col, vaf_col), "maf")
  if (length(cohort_ids) == 0L || anyNA(cohort_ids)) {
    stop("cohort_ids must contain non-missing participant identifiers.")
  }
  cohort_ids <- as.character(cohort_ids)
  if (anyDuplicated(cohort_ids)) {
    stop("cohort_ids must be unique.")
  }
  if (!is.numeric(large_clone_cutoff) ||
      length(large_clone_cutoff) != 1L ||
      is.na(large_clone_cutoff) ||
      large_clone_cutoff < 0 ||
      large_clone_cutoff > 1) {
    stop("large_clone_cutoff must be a number between 0 and 1.")
  }
  if (!is.numeric(gene_large_cutoffs) ||
      (length(gene_large_cutoffs) > 0L &&
       (is.null(names(gene_large_cutoffs)) ||
        any(names(gene_large_cutoffs) == "")))) {
    stop("gene_large_cutoffs must be a named numeric vector.")
  }

  maf_work <- maf
  maf_work[[sample_col]] <- as.character(maf_work[[sample_col]])
  maf_work[[gene_col]] <- as.character(maf_work[[gene_col]])
  maf_work[[vaf_col]] <- suppressWarnings(as.numeric(maf_work[[vaf_col]]))
  maf_work <- maf_work[
    !is.na(maf_work[[sample_col]]) &
      maf_work[[sample_col]] %in% cohort_ids &
      !is.na(maf_work[[gene_col]]),
    ,
    drop = FALSE
  ]

  if (na_vaf == "error" && anyNA(maf_work[[vaf_col]])) {
    stop("maf contains missing or non-numeric VAF values.")
  }
  if (na_vaf == "exclude") {
    maf_work <- maf_work[!is.na(maf_work[[vaf_col]]), , drop = FALSE]
  }

  carrier_table <- data.frame(
    .ukb_id = cohort_ids,
    stringsAsFactors = FALSE
  )
  make_indicator <- function(ids) {
    as.integer(carrier_table$.ukb_id %in% unique(ids))
  }

  all_ids <- unique(maf_work[[sample_col]])
  large_ids <- unique(
    maf_work[[sample_col]][
      !is.na(maf_work[[vaf_col]]) &
        maf_work[[vaf_col]] >= large_clone_cutoff
    ]
  )
  small_ids <- setdiff(all_ids, large_ids)
  carrier_table$CHIP.all <- make_indicator(all_ids)
  carrier_table$CHIP.sma <- make_indicator(small_ids)
  carrier_table$CHIP.lar <- make_indicator(large_ids)

  for (gene in unique(as.character(genes))) {
    gene_maf <- maf_work[maf_work[[gene_col]] == gene, , drop = FALSE]
    gene_all_ids <- unique(gene_maf[[sample_col]])
    cutoff <- if (gene %in% names(gene_large_cutoffs)) {
      gene_large_cutoffs[[gene]]
    } else {
      large_clone_cutoff
    }
    if (is.na(cutoff) || cutoff < 0 || cutoff > 1) {
      stop("Invalid large-clone cutoff for gene: ", gene)
    }
    gene_large_ids <- unique(
      gene_maf[[sample_col]][
        !is.na(gene_maf[[vaf_col]]) &
          gene_maf[[vaf_col]] >= cutoff
      ]
    )
    gene_small_ids <- setdiff(gene_all_ids, gene_large_ids)
    carrier_table[[paste0(gene, ".all")]] <- make_indicator(gene_all_ids)
    carrier_table[[paste0(gene, ".sma")]] <- make_indicator(gene_small_ids)
    carrier_table[[paste0(gene, ".lar")]] <- make_indicator(gene_large_ids)
  }

  colnames(carrier_table)[1] <- id_col
  carrier_table
}


#' Build a UK Biobank competing-risk outcome
#'
#' Identifies each participant's earliest matching ICD-10 event, excludes
#' prevalent and lag-period events, and creates follow-up time with death as a
#' competing event.
#'
#' @param icd10_df Long ICD-10 event data, such as the output of
#'   [prepare_ukb_icd10_events()].
#' @param phenotype_df Participant phenotype data.
#' @param icd_pattern Regular expression defining the event of interest.
#' @param lag_days Number of days after assessment during which matching events
#'   are treated as prevalent and excluded.
#' @param censor_date Administrative censoring date. This must be supplied
#'   explicitly.
#' @param id_col Participant identifier column shared by both data frames.
#' @param code_col ICD-10 code column in `icd10_df`.
#' @param diagnosis_col Optional diagnosis-description column.
#' @param event_date_col Diagnosis-date column in `icd10_df`.
#' @param assessment_date_col Baseline assessment-date column in
#'   `phenotype_df`.
#' @param death_date_col Death-date column in `phenotype_df`.
#'
#' @return The phenotype data with event information, `time`, `time_years`,
#'   numeric `status_num`, and multi-state factor `status`. Status codes are
#'   `0` for censoring, `1` for the event, and `2` for death.
#' @export
make_ukb_competing_risk_data <- function(
  icd10_df,
  phenotype_df,
  icd_pattern,
  lag_days = 180,
  censor_date,
  id_col = "eid",
  code_col = "ICD10",
  diagnosis_col = "Diagnosis",
  event_date_col = "date",
  assessment_date_col = "Date",
  death_date_col = "Date_Death"
) {
  if (!is.data.frame(icd10_df) || !is.data.frame(phenotype_df)) {
    stop("icd10_df and phenotype_df must be data frames.")
  }
  .ukb_assert_columns(
    icd10_df,
    c(id_col, code_col, event_date_col),
    "icd10_df"
  )
  .ukb_assert_columns(
    phenotype_df,
    c(id_col, assessment_date_col, death_date_col),
    "phenotype_df"
  )
  if (!is.character(icd_pattern) || length(icd_pattern) != 1L) {
    stop("icd_pattern must be a single regular expression.")
  }
  if (!is.numeric(lag_days) ||
      length(lag_days) != 1L ||
      is.na(lag_days) ||
      lag_days < 0) {
    stop("lag_days must be a single non-negative number.")
  }
  censor_date <- .ukb_as_date(censor_date)
  if (length(censor_date) != 1L || is.na(censor_date)) {
    stop("censor_date must be one valid Date value.")
  }
  if (anyDuplicated(phenotype_df[[id_col]])) {
    stop("phenotype_df must contain one row per participant.")
  }

  matched <- icd10_df[
    !is.na(icd10_df[[code_col]]) &
      grepl(icd_pattern, as.character(icd10_df[[code_col]]), perl = TRUE),
    ,
    drop = FALSE
  ]
  matched$.ukb_id <- as.character(matched[[id_col]])
  matched$.ukb_event_date <- .ukb_as_date(matched[[event_date_col]])
  if (nrow(matched) > 0L && anyNA(matched$.ukb_event_date)) {
    stop("Some matched event dates could not be parsed as Date values.")
  }

  if (nrow(matched) > 0L) {
    matched <- matched[order(matched$.ukb_id, matched$.ukb_event_date), ]
    matched <- matched[!duplicated(matched$.ukb_id), , drop = FALSE]
    event_data <- data.frame(
      .ukb_id = matched$.ukb_id,
      event_code = as.character(matched[[code_col]]),
      event_date = matched$.ukb_event_date,
      stringsAsFactors = FALSE
    )
    if (!is.null(diagnosis_col) && diagnosis_col %in% colnames(matched)) {
      event_data$event_diagnosis <- as.character(matched[[diagnosis_col]])
    }
  } else {
    event_data <- data.frame(
      .ukb_id = character(),
      event_code = character(),
      event_date = as.Date(character()),
      stringsAsFactors = FALSE
    )
    if (!is.null(diagnosis_col)) {
      event_data$event_diagnosis <- character()
    }
  }

  result <- phenotype_df
  result$.ukb_id <- as.character(result[[id_col]])
  result <- merge(result, event_data, by = ".ukb_id", all.x = TRUE, sort = FALSE)
  result$assess_date <- .ukb_as_date(result[[assessment_date_col]])
  result$death_date <- .ukb_as_date(result[[death_date_col]])
  if (anyNA(result$assess_date)) {
    stop("Some assessment dates could not be parsed as Date values.")
  }

  prevalent_or_lag <- !is.na(result$event_date) &
    result$event_date <= result$assess_date + lag_days
  result <- result[!prevalent_or_lag, , drop = FALSE]

  event_number <- as.numeric(result$event_date)
  death_number <- as.numeric(result$death_date)
  censor_number <- rep(as.numeric(censor_date), nrow(result))
  end_number <- pmin(event_number, death_number, censor_number, na.rm = TRUE)
  result$end_date <- as.Date(end_number, origin = "1970-01-01")
  result$time <- as.numeric(result$end_date - result$assess_date)
  result$time_years <- result$time / 365.25
  result$status_num <- ifelse(
    !is.na(result$event_date) & result$event_date == result$end_date,
    1L,
    ifelse(
      !is.na(result$death_date) & result$death_date == result$end_date,
      2L,
      0L
    )
  )
  result$status <- factor(
    result$status_num,
    levels = c(0, 1, 2),
    labels = c("0", "1", "2")
  )
  result <- result[!is.na(result$time) & result$time >= 0, , drop = FALSE]
  result$.ukb_id <- NULL
  rownames(result) <- NULL
  result
}


.ukb_tidy_cuminc <- function(cuminc_object, event_code, conf_level) {
  curve_names <- setdiff(names(cuminc_object), "Tests")
  suffix <- paste0(" ", event_code)
  curve_names <- curve_names[endsWith(curve_names, suffix)]
  if (length(curve_names) == 0L) {
    return(data.frame())
  }
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  do.call(rbind, lapply(curve_names, function(curve_name) {
    curve <- cuminc_object[[curve_name]]
    standard_error <- sqrt(curve$var)
    data.frame(
      group = substr(
        curve_name,
        1L,
        nchar(curve_name) - nchar(suffix)
      ),
      time = curve$time,
      estimate = curve$est,
      variance = curve$var,
      conf.low = pmax(0, curve$est - z_value * standard_error),
      conf.high = pmin(1, curve$est + z_value * standard_error),
      stringsAsFactors = FALSE
    )
  }))
}


.ukb_extract_finegray <- function(fit, exposure, model_label, conf_level) {
  fit_summary <- summary(fit)
  coefficients <- as.data.frame(fit_summary$coefficients)
  coefficient_names <- rownames(coefficients)
  selected <- startsWith(coefficient_names, ".ukb_exposure")
  if (!any(selected)) {
    stop("No exposure coefficient was found in the Fine-Gray model.")
  }
  standard_error_col <- if ("robust se" %in% colnames(coefficients)) {
    "robust se"
  } else {
    "se(coef)"
  }
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  beta <- coefficients[selected, "coef"]
  standard_error <- coefficients[selected, standard_error_col]
  raw_terms <- coefficient_names[selected]
  suffix <- sub("^\\.ukb_exposure", "", raw_terms)
  terms <- ifelse(nzchar(suffix), paste0(exposure, suffix), exposure)

  data.frame(
    Model = model_label,
    term = terms,
    sHR = exp(beta),
    conf.low = exp(beta - z_value * standard_error),
    conf.high = exp(beta + z_value * standard_error),
    p.value = coefficients[selected, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Fit a Fine-Gray competing-risk analysis
#'
#' Calculates cumulative incidence by exposure group and fits crude and
#' covariate-adjusted Fine-Gray subdistribution hazard models using expanded
#' data from [survival::finegray()].
#'
#' @param data Participant-level competing-risk data.
#' @param exposure Exposure column to compare.
#' @param covariates Character vector of adjustment covariates.
#' @param time_col Follow-up time column.
#' @param status_col Status column coded with censoring, event, and competing
#'   event values.
#' @param id_col Participant identifier column.
#' @param event_code Code for the event of interest.
#' @param censor_code Code for censoring.
#' @param conf_level Confidence level for estimates and cumulative incidence.
#' @param colors Optional named colors for exposure groups.
#' @param title Optional title used for returned plots.
#'
#' @return A list containing analysis data, cumulative-incidence results,
#'   Fine-Gray models, estimate tables, plot data, and `cif_plot` and
#'   `forest_plot` ggplot objects.
#' @export
fit_ukb_competing_risk <- function(
  data,
  exposure,
  covariates = character(),
  time_col = "time_years",
  status_col = "status_num",
  id_col = "eid",
  event_code = 1,
  censor_code = 0,
  conf_level = 0.95,
  colors = NULL,
  title = NULL
) {
  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    stop("Package 'cmprsk' is required for competing-risk analysis.")
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  .ukb_assert_columns(
    data,
    unique(c(exposure, covariates, time_col, status_col, id_col)),
    "data"
  )
  if (!is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop("conf_level must be a single number between 0 and 1.")
  }

  analysis_data <- data[
    stats::complete.cases(
      data[, unique(c(exposure, covariates, time_col, status_col, id_col)),
           drop = FALSE]
    ),
    ,
    drop = FALSE
  ]
  if (nrow(analysis_data) == 0L) {
    stop("No complete observations are available for the requested model.")
  }
  analysis_data$.ukb_time <- as.numeric(analysis_data[[time_col]])
  analysis_data$.ukb_status_num <- suppressWarnings(
    as.numeric(as.character(analysis_data[[status_col]]))
  )
  analysis_data$.ukb_exposure <- analysis_data[[exposure]]
  analysis_data$.ukb_id <- as.character(analysis_data[[id_col]])

  if (any(!is.finite(analysis_data$.ukb_time)) ||
      any(analysis_data$.ukb_time < 0)) {
    stop("time_col must contain finite, non-negative follow-up times.")
  }
  if (anyNA(analysis_data$.ukb_status_num)) {
    stop("status_col could not be converted to numeric event codes.")
  }
  if (!event_code %in% analysis_data$.ukb_status_num) {
    stop("event_code is not present in status_col.")
  }
  if (!censor_code %in% analysis_data$.ukb_status_num) {
    stop("censor_code is not present in status_col.")
  }
  if (length(unique(analysis_data$.ukb_exposure)) < 2L) {
    stop("exposure must contain at least two observed groups.")
  }
  if (anyDuplicated(analysis_data$.ukb_id)) {
    stop("id_col must uniquely identify rows in data.")
  }

  status_levels <- c(
    as.character(censor_code),
    setdiff(
      sort(unique(as.character(analysis_data$.ukb_status_num))),
      as.character(censor_code)
    )
  )
  analysis_data$.ukb_status <- factor(
    as.character(analysis_data$.ukb_status_num),
    levels = status_levels
  )

  cumulative_incidence <- cmprsk::cuminc(
    ftime = analysis_data$.ukb_time,
    fstatus = analysis_data$.ukb_status_num,
    group = analysis_data$.ukb_exposure,
    cencode = censor_code
  )
  cif_data <- .ukb_tidy_cuminc(
    cumulative_incidence,
    event_code = event_code,
    conf_level = conf_level
  )

  carry_variables <- unique(c(".ukb_exposure", covariates, ".ukb_id"))
  finegray_formula <- stats::reformulate(
    carry_variables,
    response = "survival::Surv(.ukb_time, .ukb_status)"
  )
  finegray_data <- survival::finegray(
    finegray_formula,
    data = analysis_data,
    etype = as.character(event_code)
  )

  crude_fit <- survival::coxph(
    survival::Surv(fgstart, fgstop, fgstatus) ~ .ukb_exposure,
    data = finegray_data,
    weights = finegray_data[["fgwt"]],
    ties = "efron",
    robust = TRUE,
    cluster = finegray_data$.ukb_id
  )
  adjusted_formula <- stats::reformulate(
    c(".ukb_exposure", covariates),
    response = "survival::Surv(fgstart, fgstop, fgstatus)"
  )
  adjusted_fit <- survival::coxph(
    adjusted_formula,
    data = finegray_data,
    weights = finegray_data[["fgwt"]],
    ties = "efron",
    robust = TRUE,
    cluster = finegray_data$.ukb_id
  )

  estimates <- rbind(
    .ukb_extract_finegray(
      crude_fit,
      exposure = exposure,
      model_label = "Unadjusted",
      conf_level = conf_level
    ),
    .ukb_extract_finegray(
      adjusted_fit,
      exposure = exposure,
      model_label = "Adjusted",
      conf_level = conf_level
    )
  )
  estimates$Model <- factor(
    estimates$Model,
    levels = c("Unadjusted", "Adjusted")
  )

  if (is.null(colors)) {
    groups <- unique(cif_data$group)
    colors <- stats::setNames(
      grDevices::hcl.colors(max(length(groups), 2L), "Dark 3")[seq_along(groups)],
      groups
    )
  }
  if (is.null(title)) {
    title <- paste0("Competing-risk analysis: ", exposure)
  }

  cif_plot <- ggplot(
    cif_data,
    aes(x = .data$time, y = .data$estimate, color = .data$group)
  ) +
    geom_step(linewidth = 1) +
    labs(
      title = title,
      x = "Follow-up time",
      y = "Cumulative incidence",
      color = NULL
    ) +
    scale_color_manual(values = colors) +
    theme_test()

  forest_plot <- ggplot(
    estimates,
    aes(
      x = .data$sHR,
      y = .data$term,
      color = .data$Model,
      shape = .data$Model
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbarh(
      aes(xmin = .data$conf.low, xmax = .data$conf.high),
      position = position_dodge(width = 0.5),
      height = 0.15
    ) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    scale_x_log10() +
    labs(
      title = title,
      x = paste0(
        "Subdistribution hazard ratio (",
        round(conf_level * 100),
        "% CI)"
      ),
      y = NULL,
      color = NULL,
      shape = NULL
    ) +
    theme_test()

  list(
    data = analysis_data,
    cuminc = cumulative_incidence,
    fit_crude = crude_fit,
    fit_adjusted = adjusted_fit,
    estimates = estimates,
    cif_data = cif_data,
    finegray_data = finegray_data,
    cif_plot = cif_plot,
    forest_plot = forest_plot
  )
}


.ukb_apply_omics_transform <- function(matrix_data, transform) {
  if (is.null(transform) || identical(transform, "none")) {
    return(matrix_data)
  }
  if (is.function(transform)) {
    transformed <- apply(matrix_data, 2, transform)
    transformed <- as.matrix(transformed)
    dimnames(transformed) <- dimnames(matrix_data)
    return(transformed)
  }
  if (is.character(transform) &&
      length(transform) == 1L &&
      identical(transform, "log1p")) {
    if (any(matrix_data < 0, na.rm = TRUE)) {
      stop("log1p transformation cannot be applied to negative values.")
    }
    return(log1p(matrix_data))
  }
  stop("transform must be NULL, 'none', 'log1p', or a function.")
}


#' Run covariate-adjusted limma analysis for UKB omics traits
#'
#' Fits one limma model across a selected Olink, blood-count, biochemistry, NMR,
#' or other continuous trait panel.
#'
#' @param data Participant-level data containing traits, exposure, covariates,
#'   and identifiers.
#' @param features Character vector of omics trait columns.
#' @param exposure Exposure column included in the model.
#' @param covariates Character vector of adjustment covariates.
#' @param id_col Participant identifier column.
#' @param coefficient Optional design-matrix coefficient to test. If `NULL`, an
#'   exact or unique exposure-related coefficient is selected.
#' @param transform Optional `"log1p"` transformation or a function applied to
#'   every feature column.
#' @param standardize Logical; standardize each trait before modelling.
#' @param robust Logical passed to [limma::eBayes()].
#' @param trend Logical passed to [limma::eBayes()].
#'
#' @return A list containing the limma `result`, fitted model, design matrix,
#'   analysed expression matrix, samples, features, and tested coefficient.
#' @export
run_ukb_omics_limma <- function(
  data,
  features,
  exposure = "CHIP.all",
  covariates = c(
    "Age",
    "Sex",
    "Smoking",
    "BMI",
    paste0("pc", 1:10),
    "statin",
    "antidiabetic",
    "antithrombotic"
  ),
  id_col = "eid",
  coefficient = NULL,
  transform = NULL,
  standardize = TRUE,
  robust = FALSE,
  trend = FALSE
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for UKB omics analysis.")
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  features <- unique(as.character(features))
  if (length(features) == 0L) {
    stop("features must contain at least one column name.")
  }
  .ukb_assert_columns(
    data,
    unique(c(features, exposure, covariates, id_col)),
    "data"
  )
  if (!is.logical(standardize) || length(standardize) != 1L) {
    stop("standardize must be TRUE or FALSE.")
  }

  design_columns <- unique(c(exposure, covariates))
  keep <- stats::complete.cases(data[, design_columns, drop = FALSE])
  model_data <- data[keep, , drop = FALSE]
  if (nrow(model_data) < 3L) {
    stop("Fewer than three complete observations remain for modelling.")
  }
  sample_ids <- as.character(model_data[[id_col]])
  if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
    stop("id_col must contain unique, non-missing identifiers.")
  }

  non_numeric <- features[
    !vapply(model_data[features], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "All omics features must be numeric. Non-numeric columns: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  expression <- as.matrix(model_data[, features, drop = FALSE])
  storage.mode(expression) <- "double"
  expression <- .ukb_apply_omics_transform(expression, transform)

  feature_sd <- apply(expression, 2, stats::sd, na.rm = TRUE)
  usable <- is.finite(feature_sd) & feature_sd > 0
  if (any(!usable)) {
    message(
      "Removed constant or all-missing features: ",
      paste(colnames(expression)[!usable], collapse = ", ")
    )
    expression <- expression[, usable, drop = FALSE]
  }
  if (ncol(expression) == 0L) {
    stop("No variable omics features remain for modelling.")
  }
  if (standardize) {
    expression <- scale(expression)
  }
  rownames(expression) <- sample_ids

  design <- stats::model.matrix(
    stats::reformulate(design_columns),
    data = model_data
  )
  if (is.null(coefficient)) {
    if (exposure %in% colnames(design)) {
      coefficient <- exposure
    } else {
      candidates <- colnames(design)[startsWith(colnames(design), exposure)]
      if (length(candidates) != 1L) {
        stop(
          "Could not select one exposure coefficient. Available candidates: ",
          paste(candidates, collapse = ", "),
          ". Supply coefficient explicitly."
        )
      }
      coefficient <- candidates
    }
  }
  if (!coefficient %in% colnames(design)) {
    stop(
      "coefficient is not present in the design matrix. Available: ",
      paste(colnames(design), collapse = ", ")
    )
  }

  fit <- limma::lmFit(t(expression), design)
  fit <- limma::eBayes(fit, robust = robust, trend = trend)
  result <- limma::topTable(
    fit,
    coef = coefficient,
    number = Inf,
    sort.by = "P"
  )
  result$Feature <- rownames(result)

  list(
    result = result,
    fit = fit,
    design = design,
    expression = expression,
    samples = sample_ids,
    features = colnames(expression),
    coefficient = coefficient
  )
}


.ukb_extract_beta <- function(
  effect_table,
  beta_col,
  feature_col,
  module
) {
  if (is.numeric(effect_table) && !is.null(names(effect_table))) {
    beta <- effect_table
  } else if (is.data.frame(effect_table)) {
    .ukb_assert_columns(effect_table, beta_col, paste0("effect_tables$", module))
    if (is.null(feature_col)) {
      feature_names <- rownames(effect_table)
      if (is.null(feature_names) ||
          identical(feature_names, as.character(seq_len(nrow(effect_table))))) {
        stop(
          "effect table for module '",
          module,
          "' needs informative row names or feature_col."
        )
      }
    } else {
      .ukb_assert_columns(
        effect_table,
        feature_col,
        paste0("effect_tables$", module)
      )
      feature_names <- as.character(effect_table[[feature_col]])
    }
    beta <- suppressWarnings(as.numeric(effect_table[[beta_col]]))
    names(beta) <- feature_names
  } else {
    stop(
      "Each effect table must be a named numeric vector or a data frame. ",
      "Invalid module: ",
      module
    )
  }
  beta <- beta[!is.na(beta) & !is.na(names(beta)) & nzchar(names(beta))]
  if (anyDuplicated(names(beta))) {
    stop("Duplicated effect features in module: ", module)
  }
  beta
}


.ukb_module_setting <- function(setting, module, default = NULL) {
  if (is.null(setting)) {
    return(default)
  }
  if (is.list(setting)) {
    if (!module %in% names(setting)) {
      return(default)
    }
    return(setting[[module]])
  }
  if (length(setting) > 1L && !is.null(names(setting))) {
    if (!module %in% names(setting)) {
      return(default)
    }
    return(setting[[module]])
  }
  setting
}


.ukb_winsorize <- function(x, probs) {
  if (is.null(probs) || all(is.na(x))) {
    return(x)
  }
  limits <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
}


.ukb_standardize_vector <- function(x) {
  standard_deviation <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(standard_deviation) || standard_deviation == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / standard_deviation
}


#' Build beta-weighted UKB multi-omics scores
#'
#' Preprocesses selected features for multiple omics modules, calculates
#' beta-weighted module scores, standardizes them, and optionally combines the
#' module scores into one multi-omics score.
#'
#' @param data Participant-level multi-omics data.
#' @param feature_sets Named list of feature vectors, one per omics module.
#' @param effect_tables Named list of limma result data frames or named numeric
#'   beta vectors. Names must match `feature_sets`.
#' @param beta_col Effect-size column used for data-frame effect tables.
#' @param feature_col Optional feature-name column used by every effect table.
#'   If `NULL`, data-frame row names are used.
#' @param transform Optional transformation applied to module matrices. Supply
#'   `NULL`, `"none"`, `"log1p"`, a function, or a named list/vector with one
#'   setting per module.
#' @param winsor_probs Two quantile probabilities used for winsorization, or
#'   `NULL` to disable winsorization.
#' @param impute Whether standardized missing and non-finite values are replaced
#'   with zero (`"zero"`) or retained (`"none"`).
#' @param module_weights Optional named non-negative weights used to combine
#'   standardized module scores. Equal weights are used by default.
#' @param combine Logical; calculate combined `score_multiomic` columns.
#'
#' @return A list containing augmented `data`, participant `scores`,
#'   standardized module `matrices`, aligned `betas`, and `features_used`.
#' @export
build_ukb_multiomic_score <- function(
  data,
  feature_sets,
  effect_tables,
  beta_col = "logFC",
  feature_col = NULL,
  transform = NULL,
  winsor_probs = c(0.005, 0.995),
  impute = c("zero", "none"),
  module_weights = NULL,
  combine = TRUE
) {
  impute <- match.arg(impute)
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.list(feature_sets) ||
      length(feature_sets) == 0L ||
      is.null(names(feature_sets)) ||
      any(names(feature_sets) == "")) {
    stop("feature_sets must be a non-empty named list.")
  }
  if (anyDuplicated(names(feature_sets))) {
    stop("feature_sets module names must be unique.")
  }
  if (!is.list(effect_tables) ||
      is.null(names(effect_tables)) ||
      !all(names(feature_sets) %in% names(effect_tables))) {
    stop("effect_tables must be a named list covering every feature module.")
  }
  if (!is.null(winsor_probs) &&
      (!is.numeric(winsor_probs) ||
       length(winsor_probs) != 2L ||
       anyNA(winsor_probs) ||
       winsor_probs[1] < 0 ||
       winsor_probs[2] > 1 ||
       winsor_probs[1] >= winsor_probs[2])) {
    stop("winsor_probs must contain two increasing probabilities in [0, 1].")
  }
  if (!is.logical(combine) || length(combine) != 1L || is.na(combine)) {
    stop("combine must be TRUE or FALSE.")
  }

  modules <- names(feature_sets)
  score_names <- paste0("score_", make.names(tolower(modules)))
  if (anyDuplicated(score_names)) {
    stop("Module names produce duplicated score column names.")
  }
  matrices <- vector("list", length(modules))
  betas <- vector("list", length(modules))
  features_used <- vector("list", length(modules))
  names(matrices) <- names(betas) <- names(features_used) <- modules
  scores <- data.frame(row.names = rownames(data))

  for (i in seq_along(modules)) {
    module <- modules[i]
    requested_features <- unique(as.character(feature_sets[[module]]))
    available_features <- intersect(requested_features, colnames(data))
    missing_features <- setdiff(requested_features, available_features)
    if (length(missing_features) > 0L) {
      message(
        "Missing features in module ",
        module,
        ": ",
        paste(missing_features, collapse = ", ")
      )
    }
    if (length(available_features) == 0L) {
      stop("No selected features were found for module: ", module)
    }
    non_numeric <- available_features[
      !vapply(data[available_features], is.numeric, logical(1))
    ]
    if (length(non_numeric) > 0L) {
      stop(
        "Non-numeric features in module ",
        module,
        ": ",
        paste(non_numeric, collapse = ", ")
      )
    }

    module_matrix <- as.matrix(data[, available_features, drop = FALSE])
    storage.mode(module_matrix) <- "double"
    module_transform <- .ukb_module_setting(
      transform,
      module,
      default = NULL
    )
    module_matrix <- .ukb_apply_omics_transform(
      module_matrix,
      module_transform
    )
    if (!is.null(winsor_probs)) {
      module_matrix <- apply(
        module_matrix,
        2,
        .ukb_winsorize,
        probs = winsor_probs
      )
      module_matrix <- as.matrix(module_matrix)
      colnames(module_matrix) <- available_features
    }

    z_matrix <- apply(module_matrix, 2, .ukb_standardize_vector)
    z_matrix <- as.matrix(z_matrix)
    colnames(z_matrix) <- available_features
    rownames(z_matrix) <- rownames(data)
    if (impute == "zero") {
      z_matrix[!is.finite(z_matrix)] <- 0
    }

    beta <- .ukb_extract_beta(
      effect_tables[[module]],
      beta_col = beta_col,
      feature_col = feature_col,
      module = module
    )
    common_features <- intersect(colnames(z_matrix), names(beta))
    if (length(common_features) == 0L) {
      stop(
        "No common features between data and effects for module: ",
        module
      )
    }
    beta <- beta[common_features]
    denominator <- sum(abs(beta))
    if (!is.finite(denominator) || denominator == 0) {
      stop("Module beta values sum to zero: ", module)
    }
    z_matrix <- z_matrix[, common_features, drop = FALSE]
    module_score <- as.numeric(z_matrix %*% beta) / denominator

    scores[[score_names[i]]] <- module_score
    scores[[paste0(score_names[i], "_z")]] <-
      .ukb_standardize_vector(module_score)
    matrices[[module]] <- z_matrix
    betas[[module]] <- beta
    features_used[[module]] <- common_features
  }

  if (combine) {
    module_z_columns <- paste0(score_names, "_z")
    if (is.null(module_weights)) {
      weights <- rep(1, length(modules))
      names(weights) <- modules
    } else {
      if (!is.numeric(module_weights) ||
          is.null(names(module_weights)) ||
          !all(modules %in% names(module_weights))) {
        stop("module_weights must be a named numeric vector for all modules.")
      }
      weights <- module_weights[modules]
      if (anyNA(weights) || any(weights < 0) || sum(weights) == 0) {
        stop("module_weights must be non-negative with a positive sum.")
      }
    }
    combined <- as.numeric(
      as.matrix(scores[, module_z_columns, drop = FALSE]) %*% weights
    ) / sum(weights)
    scores$score_multiomic <- combined
    scores$score_multiomic_z <- .ukb_standardize_vector(combined)
  }

  augmented_data <- cbind(data, scores)
  list(
    data = augmented_data,
    scores = scores,
    matrices = matrices,
    betas = betas,
    features_used = features_used
  )
}
