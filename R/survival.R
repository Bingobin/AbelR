# Functions for survival analyses.

#' Fit a multivariable Cox risk-score model
#'
#' Combines selected-gene expression with clinical features, fits or reuses a
#' multivariable Cox model, calculates standardized sample risk scores, and
#' divides samples into high- and low-risk groups at the median score.
#'
#' @param sample.list Character vector of sample identifiers to analyse.
#' @param gene.list Character vector of gene symbols used as expression
#'   predictors.
#' @param feature.list Character vector of clinical or molecular feature columns
#'   used as additional predictors.
#' @param tpm.matrix Expression data frame containing a `Symbol` column and
#'   sample columns.
#' @param clinic.df Clinical data frame containing sample identifiers,
#'   `overallSurvival`, `vitalStatus`, and the requested feature columns.
#' @param cox.df Optional coefficient table containing at least `Feature` and
#'   `Coefficient`. If `NULL`, coefficients are fitted from the input data.
#'
#' @return A list with `HR`, the Cox coefficient and hazard-ratio table, and
#'   `Clinc`, the merged sample risk-score and clinical table.
#' @export
RiskScore_multivar_cox <- function(
  sample.list,
  gene.list,
  feature.list,
  tpm.matrix,
  clinic.df,
  cox.df = NULL
) {
  #  sample.list <- aml_all_wt1.clinc.rm[aml_all_wt1.clinc.rm$Source.x == "BEAT",]$SampleID
  #  gene.list <- gene_sig.list
  #  tpm.matrix <- aml_all_tpm.matrix
  #  clinic.df <- aml_all_wt1.clinc
  #  feature.list <- c("WT1","TET2","IDH1","IDH2","ageAtDiagnosis")

  clinic.df$ageAtDiagnosis <- as.numeric(clinic.df$ageAtDiagnosis)
  sample.list <- gsub("-", "\\.", sample.list)
  col_index <- match(sample.list, colnames(tpm.matrix))
  row_index <- match(gene.list, tpm.matrix$Symbol)
  risk.tpm <- tpm.matrix[row_index, col_index, ]
  rownames(risk.tpm) <- tpm.matrix$Symbol[row_index]
  risk.tpm <- t(risk.tpm)
  rownames(risk.tpm) <- gsub("\\.", "-", rownames(risk.tpm))
  risk.cox <- merge(risk.tpm, clinic.df, by.x = 0, by.y = 2, all = FALSE)
  risk.cox[risk.cox$vitalStatus == "Death", ]$vitalStatus <- 1
  risk.cox[risk.cox$vitalStatus == "Alive", ]$vitalStatus <- 0
  risk.cox$vitalStatus <- as.numeric(risk.cox$vitalStatus)
  risk.cox$overallSurvival <- as.numeric(risk.cox$overallSurvival)
  risk.cox <- risk.cox %>% filter(overallSurvival > 0)
  #  survpre_mut <- data.table::as.data.table(risk.cox[,c("WT1","TET2","IDH1","IDH2","ageAtDiagnosis")])
  #  colnames(survpre_mut) <- c("WT1_mut","TET2_mut", "IDH1_mut", "IDH2_mut","Age")
  #  survpre_mut$Age <- as.numeric(survpre_mut$Age)
  survpre_mut <- data.table::as.data.table(risk.cox[,
    feature.list,
    drop = FALSE
  ])
  survpre_sig <- risk.cox[, gene.list, drop = FALSE]
  survpre_os <- risk.cox[, c("overallSurvival", "vitalStatus")]
  colnames(survpre_os) <- c("Time", "Status")
  surv_model.df <- cbind(survpre_sig, survpre_mut)
  surv_model.df <- as.data.frame(
    row.names = risk.cox$Row.names,
    cbind(surv_model.df, survpre_os)
  )
  surv_model.df <- na.omit(surv_model.df)
  survobj <- Surv(surv_model.df$Time, surv_model.df$Status)
  survpre <- as.matrix(surv_model.df[,
    c(gene.list, feature.list),
    drop = FALSE
  ])
  if (is.null(cox.df)) {
    model <- coxph(survobj ~ survpre)
    cox.df.sum <- summary(model)
    cox.df <- as.data.frame(cbind(cox.df.sum$coefficients, cox.df.sum$conf.int))
    cox.df$Feature <- gsub("survpre", "", rownames(cox.df))
    cox.df <- data.table::as.data.table(cox.df[, c(10, 1, 5, 6, 8, 9)])
    colnames(cox.df) <- c(
      "Feature",
      "Coefficient",
      "Pvalue",
      "Hazard.ratio",
      "HR.Lower.95",
      "HR.Upper.95"
    )
    cox.df <- na.omit(cox.df)
  } else {
    required_cox_cols <- c("Feature", "Coefficient")
    if (!all(required_cox_cols %in% colnames(cox.df))) {
      stop("cox.df must contain Feature and Coefficient columns.")
    }
  }
  ls_df <- lapply(1:nrow(surv_model.df), function(x) {
    sample <- rownames(surv_model.df)[x]
    score <- sum(sapply(1:nrow(cox.df), function(i) {
      coef <- cox.df$Coefficient[i]
      tpm <- surv_model.df[x, colnames(surv_model.df) == cox.df$Feature[i]]
      coef * tpm
    }))
    data.frame(Sample = sample, Risk_score = score)
  })
  rs.df <- do.call(rbind, ls_df)
  rs.df$RS_norm <- sapply(1:nrow(rs.df), function(i) {
    (rs.df$Risk_score[i] - mean(rs.df$Risk_score)) / sd(rs.df$Risk_score)
  })
  rs.df$Risk_Group <- "High_risk"
  rs.df[rs.df$RS_norm < median(rs.df$RS_norm), ]$Risk_Group <- "Low_risk"
  rs.df$Sample <- gsub("\\.", "-", rs.df$Sample)
  rownames(surv_model.df) <- gsub("\\.", "-", rownames(surv_model.df))
  rs.df.clinc <- merge(rs.df, surv_model.df, by.x = 1, by.y = 0, all = FALSE)
  rs.df.clinc <- merge(
    rs.df.clinc,
    clinic.df,
    by.x = 1,
    by.y = 2,
    all = FALSE
  )
  result.ls <- list()
  result.ls[["HR"]] <- cox.df
  result.ls[["Clinc"]] <- rs.df.clinc
  return(result.ls)
}


#' Calculate gene-based Cox risk groups
#'
#' Fits a univariable Cox model for each selected gene, combines the fitted
#' coefficients into a sample risk score, standardizes the score, and assigns
#' median-based risk groups.
#'
#' @param gene.list Character vector of gene-expression columns.
#' @param design.df Data frame containing `SampleID`, `Time`, `Status`, and the
#'   selected gene-expression columns.
#'
#' @return The input design data merged with `Risk_score`, `RS_norm`, and
#'   `Risk_Group` columns.
#' @export
Risk_grouping_cox <- function(gene.list, design.df) {
  #gene.list<-risk_gene.list
  #design.df<-Risk_model.beat.df

  ls_df <- lapply(match(gene.list, colnames(design.df)), function(i) {
    pretty_aml <- data.frame(
      time = design.df$Time,
      status = design.df$Status,
      TPM = design.df[, i]
    )
    Gene_label <- colnames(design.df)[i]
    model <- coxph(Surv(time, status) ~ TPM, data = pretty_aml)
    cox_result <- summary(model)
    data.frame(
      row.names = Gene_label,
      cox_result$conf.int,
      PV = cox_result$coefficients[5]
    )
  })
  hr.df <- do.call(rbind, ls_df)
  ls_df <- lapply(1:nrow(design.df), function(x) {
    sample <- design.df$SampleID[x]
    score <- sum(sapply(1:length(gene.list), function(i) {
      coef <- log(hr.df[rownames(hr.df) == gene.list[i], ]$exp.coef.)
      tpm <- design.df[x, colnames(design.df) == gene.list[i]]
      coef * tpm
    }))
    data.frame(Sample = sample, Risk_score = score)
  })
  rs.df <- do.call(rbind, ls_df)
  rs.df$RS_norm <- sapply(1:nrow(rs.df), function(i) {
    (rs.df$Risk_score[i] - mean(rs.df$Risk_score)) / sd(rs.df$Risk_score)
  })
  rs.df$Risk_Group <- "High_risk"
  rs.df[rs.df$RS_norm < median(rs.df$RS_norm), ]$Risk_Group <- "Low_risk"
  rs.df$Sample <- gsub("\\.", "-", rs.df$Sample)
  design.df$SampleID <- gsub("\\.", "-", design.df$SampleID)
  rs.df.clinc <- merge(rs.df, design.df, by.x = 1, by.y = 1, all = FALSE)
  return(rs.df.clinc)
}


#' Prepare expression and survival data for risk modelling
#'
#' Extracts selected genes from an expression table, joins their expression to
#' survival metadata, removes invalid survival records, and converts survival
#' status to numeric event coding.
#'
#' @param gene.list Character vector of gene symbols to extract.
#' @param design.df Survival metadata containing `SampleID`, `Time`, and
#'   `Status`.
#' @param tpm.df Expression data frame containing a `Symbol` column and sample
#'   columns.
#'
#' @return A data frame with one row per matched sample, selected-gene
#'   expression, numeric `Time`, and binary `Status`.
#' @export
Risk_model_df <- function(gene.list, design.df, tpm.df) {
  #gene.list <- risk_gene.list
  #design.df <- aml_design_TCGA
  #tpm.df <- aml_combind_tpm.df
  index_col <- match(design.df$SampleID, colnames(tpm.df))
  index_row <- match(gene.list, tpm.df$Symbol)
  Risk.df <- tpm.df[index_row, index_col]
  rownames(Risk.df) <- gene.list
  Risk.df <- t(Risk.df)
  Risk.df <- merge(Risk.df, design.df, by.x = 0, by.y = 1)
  Risk.df <- Risk.df %>%
    filter(!(is.na(Status) | Status == "Unknown" | Time == 0))
  colnames(Risk.df)[1] <- "SampleID"
  Risk.df$Time <- as.numeric(as.character(Risk.df$Time))
  Risk.df$Status <- as.character(Risk.df$Status)
  Risk.df[Risk.df$Status == "Death", ]$Status <- 1
  Risk.df[Risk.df$Status == "Alive", ]$Status <- 0
  Risk.df$Status <- as.numeric(Risk.df$Status)
  return(Risk.df)
}


#' Plot risk scores and survival status
#'
#' Draws vertically aligned panels showing ranked risk scores and survival time
#' or status for samples assigned to high- and low-risk groups.
#'
#' @param clinic_rs.df Data frame containing `Risk_score`, `Risk_Group`, `Time`,
#'   and `vitalStatus` columns.
#'
#' @return A combined plot produced by [cowplot::plot_grid()].
#' @export
plot_survial_risk <- function(clinic_rs.df) {
  #clinic_rs.df<-rs.df.clinc
  blank <- theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank()
  )
  clinic_rs.df$rank <- rank(clinic_rs.df$Risk_score)
  x1 <- ggplot(clinic_rs.df, aes(x = rank, y = Risk_score)) +
    geom_point(aes(color = Risk_Group)) +
    blank +
    theme(
      legend.justification = c("right", "top"),
      legend.position = c(1, 1),
      legend.title = element_blank()
    ) +
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    ) +
    #    scale_x_continuous(breaks=c(dim(clinic_rs.df)[1]*0.25, dim(clinic_rs.df)[1]*0.75),labels=c("Low_risk","High_risk")) +
    geom_vline(xintercept = dim(clinic_rs.df)[1] * 0.5, linetype = 5) +
    scale_color_manual(values = c(c("#FF7F00", "#1F78B4"))) +
    ylab("Risk score") +
    theme(
      panel.border = element_rect(linetype = 1, size = 1, fill = NA),
      axis.line = element_blank(),
      axis.text.x = element_blank()
    ) +
    theme(
      legend.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA)
    )

  x2 <- ggplot(clinic_rs.df, aes(x = rank, y = Time)) +
    geom_point(aes(color = vitalStatus, shape = vitalStatus)) +
    blank +
    theme(
      legend.justification = c("right", "top"),
      legend.position = c(1, 1),
      legend.title = element_blank()
    ) +
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    ) +
    scale_x_continuous(
      breaks = c(dim(clinic_rs.df)[1] * 0.25, dim(clinic_rs.df)[1] * 0.75),
      labels = c("Low_risk", "High_risk")
    ) +
    geom_vline(xintercept = dim(clinic_rs.df)[1] * 0.5, linetype = 5) +
    scale_color_manual(values = c("grey40", "black")) +
    scale_shape_manual(values = c(1, 2)) +
    ylab("Survial days") +
    theme(
      panel.border = element_rect(linetype = 1, size = 1, fill = NA),
      axis.line = element_blank()
    ) +
    theme(
      legend.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA)
    )

  return(plot_grid(x1, x2, ncol = 1))
}


#' Plot time-dependent survival ROC curves
#'
#' Calculates survival ROC curves at one, three, and five years and displays
#' their AUC values in a combined plot.
#'
#' @param risk.df Data frame containing survival `Time`, binary `Status`, and
#'   `Risk_score` columns.
#'
#' @return A [ggplot2::ggplot] object containing the three ROC curves.
#' @export
plot_roc_curve <- function(risk.df) {
  #  risk.df <- rs.df.clinc
  #  risk.df$lp <- a
  tmp.df.1 <- survivalROC(
    Stime = risk.df$Time,
    status = risk.df$Status,
    marker = risk.df$Risk_score,
    predict.time = 365 * 1,
    span = 0.1 * nrow(risk.df)^(-0.30)
    #                 method="KM"
  )
  tmp.df.3 <- survivalROC(
    Stime = risk.df$Time,
    status = risk.df$Status,
    marker = risk.df$Risk_score,
    predict.time = 365 * 3,
    span = 0.1 * nrow(risk.df)^(-0.30)
    #                 method="KM"
  )
  tmp.df.5 <- survivalROC(
    Stime = risk.df$Time,
    status = risk.df$Status,
    marker = risk.df$Risk_score,
    predict.time = 365 * 5,
    span = 0.1 * nrow(risk.df)^(-0.30)
    #                 method="KM"
  )
  #  str(tmp.df)
  auc1 <- paste0("1 years", "(AUC=", round(tmp.df.1$AUC, 4), ")")
  auc3 <- paste0("3 years", "(AUC=", round(tmp.df.3$AUC, 4), ")")
  auc5 <- paste0("5 years", "(AUC=", round(tmp.df.5$AUC, 4), ")")
  tmp.df.plot.1 <- data.frame(
    FP = sort(tmp.df.1$FP),
    TP = sort(tmp.df.1$TP),
    PT = rep(auc1, length(tmp.df.1$TP))
  )
  tmp.df.plot.3 <- data.frame(
    FP = sort(tmp.df.3$FP),
    TP = sort(tmp.df.3$TP),
    PT = rep(auc3, length(tmp.df.3$TP))
  )
  tmp.df.plot.5 <- data.frame(
    FP = sort(tmp.df.5$FP),
    TP = sort(tmp.df.5$TP),
    PT = rep(auc5, length(tmp.df.5$TP))
  )
  #  tmp.df.plot.1 <- data.frame(FP=tmp.df.1$FP,TP=tmp.df.1$TP, PT=rep(auc1,length(tmp.df.1$TP)))
  #  tmp.df.plot.3 <- data.frame(FP=tmp.df.3$FP,TP=tmp.df.3$TP, PT=rep(auc3,length(tmp.df.3$TP)))
  #  tmp.df.plot.5 <- data.frame(FP=tmp.df.5$FP,TP=tmp.df.5$TP, PT=rep(auc5,length(tmp.df.5$TP)))
  tmp.df.plot <- rbind(tmp.df.plot.1, tmp.df.plot.3)
  tmp.df.plot <- rbind(tmp.df.plot, tmp.df.plot.5)

  ggplot(tmp.df.plot, aes(x = FP, y = TP)) +
    geom_line(aes(color = PT), size = 0.5) +
    blank +
    geom_abline(slope = 1, linetype = 2, size = 0.5) +
    xlab("1-Specificity") +
    ylab("Sensitivity") +
    scale_color_manual(values = brewer.pal(9, "Set1")[c(2, 4, 5)]) +
    theme(
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = c(.95, .05),
      legend.justification = c("right", "bottom")
    ) +
    theme(
      panel.border = element_rect(linetype = 1, size = 0.8, fill = NA),
      axis.line = element_blank()
    )
}


#' Plot odds ratios from logistic regression
#'
#' Fits a multivariable binary logistic regression model, extracts odds ratios
#' with Wald confidence intervals, and draws a forest plot. Predictor names and
#' display labels are supplied by the caller, so the function can be used with
#' UK Biobank data or other cohort data.
#'
#' @param data A data frame containing the outcome and predictor variables.
#' @param outcome Character scalar naming a binary outcome column. Numeric
#'   outcomes must be coded as `0` and `1`; logical and two-level factor
#'   outcomes are also supported.
#' @param predictors Character vector naming predictor columns.
#' @param term_labels Optional named character vector used to relabel model
#'   coefficient terms. Names must be coefficient names, including factor
#'   levels where applicable, and values are the labels displayed in the plot.
#' @param conf_level Confidence level for the Wald confidence intervals.
#' @param digits Number of decimal places used in plot labels.
#' @param point_color Color used for points and confidence intervals.
#' @param title Optional plot title. By default, the outcome name is shown.
#' @param show_labels Logical; if `TRUE`, display odds ratios, confidence
#'   intervals, and significance symbols beside the estimates.
#'
#' @return A list containing `model`, the fitted [stats::glm()] object; `or`, a
#'   data frame of coefficient-level odds ratios and confidence intervals; and
#'   `plot`, the forest plot.
#'
#' @examples
#' fit <- plot_logistic_forest(
#'   data = mtcars,
#'   outcome = "am",
#'   predictors = c("wt", "hp"),
#'   term_labels = c(
#'     wt = "Weight (1000 lbs)",
#'     hp = "Horsepower"
#'   )
#' )
#' fit$or
#' fit$plot
#'
#' @export
plot_logistic_forest <- function(
  data,
  outcome,
  predictors,
  term_labels = NULL,
  conf_level = 0.95,
  digits = 2,
  point_color = "firebrick",
  title = NULL,
  show_labels = TRUE
) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.character(outcome) || length(outcome) != 1L || is.na(outcome)) {
    stop("outcome must be a single, non-missing column name.")
  }
  if (!is.character(predictors) || length(predictors) == 0L) {
    stop("predictors must be a non-empty character vector.")
  }

  predictors <- unique(predictors)
  required_cols <- c(outcome, predictors)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0L) {
    stop(
      "The following columns are missing from data: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (outcome %in% predictors) {
    stop("outcome must not also be included in predictors.")
  }
  if (!is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop("conf_level must be a single number between 0 and 1.")
  }
  if (!is.numeric(digits) ||
      length(digits) != 1L ||
      is.na(digits) ||
      digits < 0 ||
      digits != as.integer(digits)) {
    stop("digits must be a non-negative integer.")
  }
  if (!is.logical(show_labels) ||
      length(show_labels) != 1L ||
      is.na(show_labels)) {
    stop("show_labels must be TRUE or FALSE.")
  }
  if (!is.null(term_labels) &&
      (!is.character(term_labels) ||
       is.null(names(term_labels)) ||
       any(names(term_labels) == ""))) {
    stop("term_labels must be a named character vector.")
  }

  outcome_values <- data[[outcome]]
  observed_outcome <- outcome_values[!is.na(outcome_values)]
  valid_outcome <- if (is.factor(observed_outcome)) {
    nlevels(droplevels(observed_outcome)) == 2L
  } else if (is.logical(observed_outcome)) {
    length(unique(observed_outcome)) == 2L
  } else if (is.numeric(observed_outcome) || is.integer(observed_outcome)) {
    length(unique(observed_outcome)) == 2L &&
      all(unique(observed_outcome) %in% c(0, 1))
  } else {
    FALSE
  }
  if (!valid_outcome) {
    stop(
      "outcome must contain two observed levels; numeric outcomes must use ",
      "0 and 1."
    )
  }

  fit <- stats::glm(
    formula = stats::reformulate(predictors, response = outcome),
    data = data,
    family = stats::binomial()
  )

  coef_matrix <- summary(fit)$coefficients
  coef_names <- setdiff(rownames(coef_matrix), "(Intercept)")
  if (length(coef_names) == 0L) {
    stop("The fitted model contains no non-intercept coefficients.")
  }

  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  estimates <- coef_matrix[coef_names, "Estimate"]
  standard_errors <- coef_matrix[coef_names, "Std. Error"]
  p_values <- coef_matrix[coef_names, "Pr(>|z|)"]

  result <- data.frame(
    term = coef_names,
    estimate = exp(estimates),
    conf.low = exp(estimates - z_value * standard_errors),
    conf.high = exp(estimates + z_value * standard_errors),
    p.value = p_values,
    stringsAsFactors = FALSE,
    row.names = coef_names
  )

  result$display_term <- result$term
  if (!is.null(term_labels)) {
    matched_terms <- intersect(result$term, names(term_labels))
    result$display_term[match(matched_terms, result$term)] <-
      unname(term_labels[matched_terms])
  }

  result$stars <- ifelse(
    is.na(result$p.value),
    NA_character_,
    ifelse(
      result$p.value < 0.0001,
      "****",
      ifelse(
        result$p.value < 0.001,
        "***",
        ifelse(
          result$p.value < 0.01,
          "**",
          ifelse(result$p.value < 0.05, "*", "ns")
        )
      )
    )
  )
  result$label <- paste0(
    "OR=",
    formatC(result$estimate, format = "f", digits = digits),
    " [",
    formatC(result$conf.low, format = "f", digits = digits),
    "-",
    formatC(result$conf.high, format = "f", digits = digits),
    "] ",
    result$stars
  )

  if (is.null(title)) {
    title <- paste0("Logistic regression (", outcome, ")")
  }

  plot <- ggplot(
    result,
    aes(
      x = stats::reorder(display_term, estimate),
      y = estimate,
      ymin = conf.low,
      ymax = conf.high
    )
  ) +
    geom_pointrange(color = point_color) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    coord_flip() +
    labs(
      y = paste0(
        "Odds Ratio (",
        round(conf_level * 100),
        "% CI)"
      ),
      x = NULL,
      title = title
    ) +
    theme_test()

  if (show_labels) {
    plot <- plot +
      geom_text(
        aes(label = label),
        vjust = -1.2,
        color = "black",
        na.rm = TRUE
      )
  }

  list(
    model = fit,
    or = result,
    plot = plot
  )
}
