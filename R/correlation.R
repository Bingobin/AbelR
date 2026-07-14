# Functions for correlation analyses.

#' Plot the correlation between two genes
#'
#' Extracts two genes from an expression table, applies `log2(x + 1)`, computes
#' a Pearson correlation and linear regression, and draws a labelled scatter
#' plot.
#'
#' @param EXdata Expression data frame with a `Symbol` column followed by sample
#'   expression columns. The first three columns are treated as annotation.
#' @param Gene1,Gene2 Gene symbols to compare.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
ccor_paired_gene <- function(EXdata, Gene1, Gene2) {
  #  EXdata<-aml_all_tpm.matrix
  #  Gene1 <- "TET2"
  #  Gene2 <- "TP53"
  ##########################
  #           Symbol Length           Type    BA2000    BA2003      BA2004
  #TSPAN6     TSPAN6   4535 protein_coding   0.00000  1.016602   0.3353875
  #TNMD         TNMD   1610 protein_coding   0.00000  0.000000   0.0000000
  ########################
  corr.df <- data.frame(
    Gene1 = log2(
      as.numeric(EXdata[EXdata$Symbol == Gene1, 4:ncol(EXdata)]) + 1
    ),
    Gene2 = log2(as.numeric(EXdata[EXdata$Symbol == Gene2, 4:ncol(EXdata)]) + 1)
  )
  corr.result <- cor.test(~ Gene1 + Gene2, corr.df, method = "pearson")
  text <- paste(
    "k=",
    round(lm(Gene1 ~ Gene2, corr.df)$coefficients[2], 4),
    ",",
    "r=",
    round(corr.result$estimate, 4),
    sep = ""
  )
  pv <- corr.result$p.value
  if (pv < 0.05 & pv > 0.01) {
    sig <- "*"
  } else if (pv < 0.01 & pv > 0.001) {
    sig <- "**"
  } else if (pv < 0.001 & pv > 0.0001) {
    sig <- "***"
  } else if (pv < 0.0001) {
    sig <- "****"
  } else {
    sig <- "na"
  }
  blank <- theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank()
  )
  ggplot(corr.df, aes(y = Gene2, x = Gene1)) +
    geom_point(size = 1, alpha = 0.5, color = "grey") +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    ggtitle(text) +
    annotate(
      "text",
      label = sig,
      y = max(corr.df$Gene2, na.rm = TRUE),
      x = median(corr.df[corr.df$Gene1 != 0, ]$Gene1, na.rm = TRUE),
      size = 7
    ) +
    xlab(Gene1) +
    ylab(Gene2) +
    blank
}


#' Plot the correlation between two variables
#'
#' Computes a Pearson correlation and linear regression between two columns and
#' colours samples by a grouping column.
#'
#' @param EXdata Data frame containing the variables and grouping column.
#' @param obj1,obj2 Character column names of the variables to compare.
#' @param Group Character name of the column used to colour points.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
ccor_paired_obj <- function(EXdata, obj1, obj2, Group) {
  corr.df <- data.frame(
    obj1 = EXdata[, obj1],
    obj2 = EXdata[, obj2],
    Group = EXdata[, Group]
  )
  corr.result <- cor.test(~ obj1 + obj2, corr.df, method = "pearson")
  text <- paste(
    "k=",
    round(lm(obj1 ~ obj2, corr.df)$coefficients[2], 4),
    ",",
    "r=",
    round(corr.result$estimate, 4),
    sep = ""
  )
  pv <- corr.result$p.value
  if (pv < 0.05 & pv > 0.01) {
    sig <- "*"
  } else if (pv < 0.01 & pv > 0.001) {
    sig <- "**"
  } else if (pv < 0.001 & pv > 0.0001) {
    sig <- "***"
  } else if (pv < 0.0001) {
    sig <- "****"
  } else {
    sig <- "na"
  }
  blank <- theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank()
  )
  ggplot(corr.df, aes(x = obj2, y = obj1)) +
    geom_point(size = 3, aes(color = Group)) +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    ggtitle(text) +
    annotate(
      "text",
      label = sig,
      y = max(corr.df$obj1, na.rm = TRUE),
      x = min(corr.df$obj2) + (max(corr.df$obj2) - min(corr.df$obj2)) / 2,
      size = 7
    ) +
    theme_test() +
    xlab(obj2) +
    ylab(obj1) +
    scale_color_manual(values = c("#08537C", "#A81E2C"))
}

