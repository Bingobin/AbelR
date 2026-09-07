# Functions for deg-plots analyses.

#' Draw a heatmap of DESeq2 differentially expressed genes
#'
#' Extracts variance-stabilized expression for significant genes, scales each
#' gene across samples, labels selected and top-ranked genes, and draws a
#' ComplexHeatmap with sample-group annotation. A message is shown when either
#' direction has no genes; the heatmap is skipped when fewer than two genes are
#' available.
#'
#' @param deg_result AbelR DESeq2 result list containing `result`, `vsd`, and
#'   `design` components.
#' @param sample_name Character vector specifying the sample column order.
#' @param label_genes Optional gene symbols to label in addition to top genes.
#' @param top_n Number of top upregulated and downregulated genes to label.
#'
#' @return Invisibly returns the [ComplexHeatmap::Heatmap] object after drawing
#'   it on the active graphics device, or `NULL` when fewer than two genes are
#'   available.
#' @export
plot_deg_heatmap_for_DEGseq2 <- function(
  deg_result,
  sample_name,
  label_genes = NULL,
  top_n = 10
) {
  up_genes <- unique(deg_result$up$Row.names)
  dw_genes <- unique(deg_result$dw$Row.names)
  up_genes <- up_genes[up_genes %in% rownames(deg_result$vsd)]
  dw_genes <- dw_genes[dw_genes %in% rownames(deg_result$vsd)]
  dw_genes <- setdiff(dw_genes, up_genes)
  deg_genes <- c(up_genes, dw_genes)

  up_n <- nrow(deg_result$up)
  dw_n <- nrow(deg_result$dw)
  if (up_n == 0 && dw_n == 0) {
    message(
      sample_name,
      ": no upregulated or downregulated DEGs; skip heatmap."
    )
    return(invisible(NULL))
  }
  if (up_n == 0) {
    message(
      sample_name,
      ": no upregulated DEGs; plotting downregulated DEGs only."
    )
  }
  if (dw_n == 0) {
    message(
      sample_name,
      ": no downregulated DEGs; plotting upregulated DEGs only."
    )
  }
  if (length(deg_genes) < 2) {
    message(
      sample_name,
      ": fewer than 2 DEGs are available in the VST matrix; skip heatmap."
    )
    return(invisible(NULL))
  }

  design <- deg_result$design
  group_col <- if ("ContrastGroup" %in% colnames(design)) {
    "ContrastGroup"
  } else {
    "Group"
  }
  group_levels <- levels(factor(design[[group_col]]))
  ctr_group <- group_levels[1]
  tr_group <- group_levels[2]
  col_anno <- data.frame(
    row.names = design$SampleID,
    Group = factor(design[[group_col]], levels = group_levels)
  )
  col_anno_color <- list(
    Group = c("#C81E76", "#4DBBD5")[seq_along(group_levels)]
  )
  names(col_anno_color$Group) <- levels(col_anno$Group)

  EX_data <- SummarizedExperiment::assay(deg_result$vsd[deg_genes, ])
  sd_rows <- apply(EX_data, 1, sd)
  EX_data <- EX_data[sd_rows > 0, , drop = FALSE]
  if (nrow(EX_data) < 2) {
    message(sample_name, ": fewer than 2 variable DEGs, skip heatmap.")
    return(invisible(NULL))
  }

  EX_data <- t(scale(t(EX_data)))
  EX_data[EX_data > 2] <- 2
  EX_data[EX_data < -2] <- -2

  row_label_df <- deg_result$result |>
    select(any_of(c(
      "Row.names",
      "Symbol",
      "padj",
      "pvalue",
      "log2FoldChange"
    ))) |>
    distinct()
  row_label_df <- row_label_df[
    match(rownames(EX_data), row_label_df$Row.names),
  ]
  if ("Symbol" %in% colnames(row_label_df)) {
    row_symbols <- row_label_df$Symbol
    row_symbols[is.na(row_symbols)] <- rownames(EX_data)[is.na(row_symbols)]
  } else {
    row_symbols <- rownames(EX_data)
  }

  label_genes <- unique(label_genes)
  if (!is.null(top_n) && top_n > 0) {
    rank_cols <- intersect(c("padj", "pvalue"), colnames(row_label_df))
    fc_cols <- intersect(
      c("log2FoldChange", "avg_log2FC", "logFC"),
      colnames(row_label_df)
    )
    if (length(rank_cols) > 0 && length(fc_cols) > 0) {
      rank_col <- rank_cols[1]
      fc_col <- fc_cols[1]
      top_up_genes <- row_label_df |>
        filter(!is.na(.data[[rank_col]])) |>
        filter(!is.na(.data[[fc_col]])) |>
        filter(.data[[fc_col]] > 0) |>
        arrange(.data[[rank_col]]) |>
        slice_head(n = top_n) |>
        pull(Row.names)
      top_down_genes <- row_label_df |>
        filter(!is.na(.data[[rank_col]])) |>
        filter(!is.na(.data[[fc_col]])) |>
        filter(.data[[fc_col]] < 0) |>
        arrange(.data[[rank_col]]) |>
        slice_head(n = top_n) |>
        pull(Row.names)
      label_genes <- unique(c(label_genes, top_up_genes, top_down_genes))
    }
  }

  row_labels <- rep("", nrow(EX_data))
  label_index <- rownames(EX_data) %in%
    label_genes |
    row_symbols %in% label_genes
  row_labels[label_index] <- row_symbols[label_index]

  title <- paste0(
    sample_name,
    " | ",
    tr_group,
    " vs ",
    ctr_group,
    " | Up: ",
    up_n,
    " Down: ",
    dw_n
  )

  if (
    !requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)
  ) {
    stop("ComplexHeatmap and circlize are required for marked gene labels.")
  }

  label_at <- which(row_labels != "")
  right_anno <- NULL
  if (length(label_at) > 0) {
    right_anno <- ComplexHeatmap::rowAnnotation(
      Mark = ComplexHeatmap::anno_mark(
        at = label_at,
        labels = row_labels[label_at],
        labels_gp = grid::gpar(fontsize = 9)
      )
    )
  }

  heatmap_plot <- ComplexHeatmap::Heatmap(
    EX_data,
    name = "Z-score",
    col = circlize::colorRamp2(
      c(-2, 0, 2),
      c("#63010b", "white", "#023c5c")
    ),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      Group = col_anno$Group,
      col = col_anno_color
    ),
    right_annotation = right_anno,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    column_title = title
  )
  grid::grid.newpage()
  ComplexHeatmap::draw(
    heatmap_plot,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  invisible(heatmap_plot)
}


#' Select genes retained in a DESeq2 volcano plot
#'
#' Combines all significant genes, a random background sample, and explicitly
#' requested genes into the set used by [volcano_plot_Deseq2()].
#'
#' @param deseq2_result.df Annotated DESeq2 result data frame containing
#'   `Symbol`, `log2FoldChange`, and `padj`.
#' @param gene.list Character vector of gene symbols that must be retained.
#' @param n Maximum number of nonsignificant background genes to sample.
#' @param pc Logical; restrict candidates to protein-coding genes.
#' @param pv Adjusted P-value threshold.
#' @param fc Fold-change threshold on the linear scale.
#' @param p_col Significance column used to select genes.
#'
#' @return A character vector of selected gene symbols.
#' @export
target_for_volcano <- function(
  deseq2_result.df,
  gene.list,
  n = 5000,
  pc = FALSE,
  pv = 0.05,
  fc = 1.5,
  p_col = "padj"
) {
  # deseq2_result.df <- RELA_MUTvsWT_NON.DEGs$result
  if (pc == TRUE) {
    deseq2_result.df <- deseq2_result.df |>
      filter(Gene_Type == "protein_coding")
  }
  if (
    length(p_col) != 1 ||
      is.na(p_col) ||
      !p_col %in% colnames(deseq2_result.df)
  ) {
    stop("p_col must name a column in deseq2_result.df.")
  }
  significance <- deseq2_result.df[[p_col]]
  valid <- !is.na(deseq2_result.df$log2FoldChange) & !is.na(significance)
  deg_index <- valid &
    abs(deseq2_result.df$log2FoldChange) > log2(fc) &
    significance < pv
  deg.list <- deseq2_result.df$Symbol[deg_index]
  other.list <- deseq2_result.df$Symbol[
    valid & !deseq2_result.df$Symbol %in% deg.list
  ]
  if (
    length(n) != 1 ||
      !is.numeric(n) ||
      is.na(n) ||
      !is.finite(n) ||
      n < 0
  ) {
    stop("n must be one non-negative number.")
  }
  sample_n <- min(as.integer(n), length(other.list))
  set.seed(123)
  random.list <- if (sample_n > 0) {
    sample(other.list, sample_n)
  } else {
    character(0)
  }
  result.list <- c(deg.list, random.list)
  result.list <- c(result.list, gene.list)
  return(unique(result.list))
}


#' Draw a volcano plot from DESeq2 results
#'
#' Creates a volcano plot from a DESeq2 result table, classifies significantly
#' upregulated and downregulated genes, and labels selected genes together with
#' the most significant genes in each direction. The plot is still returned
#' when either or both DEG directions are empty; counts are reported in a
#' message and in the plot subtitle.
#' For plotting, zero P values are replaced by the smallest finite positive
#' value in the selected P-value column of the full input table, before applying
#' `max_y`. Significance classification and top-gene ranking use original values.
#' If no positive value exists, a warning is issued and zeros use the existing
#' `max_y` cap.
#'
#' @param deseq2_result.df A data frame containing at least `Symbol`,
#'   `log2FoldChange`, `pvalue`, and `padj` columns. A `Gene_Type` column is
#'   also required when `pc = TRUE`.
#' @param gene.list A character vector of gene symbols that should be retained
#'   and labelled in the plot.
#' @param n Number of nonsignificant background genes sampled for plotting.
#' @param pc Logical; if `TRUE`, restrict candidate genes to those with
#'   `Gene_Type == "protein_coding"`.
#' @param pv P-value or adjusted P-value threshold used to define significance.
#' @param fc Fold-change threshold on the linear scale. Vertical reference lines
#'   are drawn at `+/-log2(fc)`.
#' @param max_x Maximum absolute log2 fold change displayed. More extreme values
#'   are capped at this limit.
#' @param max_y Maximum displayed `-log10` significance value.
#' @param top Number of genes automatically labelled in each fold-change
#'   direction. Genes are ranked first by ascending P value, then by descending
#'   absolute log2 fold change to break P-value ties.
#' @param adjust Logical; use `padj` when `TRUE` and `pvalue` when `FALSE`.
#' @param up_color Color used for significantly upregulated genes.
#' @param down_color Color used for significantly downregulated genes.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @seealso [target_for_volcano()]
#' @export
volcano_plot_Deseq2 <- function(
  deseq2_result.df,
  gene.list,
  n = 5000,
  pc = FALSE,
  pv = 0.05,
  fc = 1.5,
  max_x = 100,
  max_y = 100,
  top = 5,
  adjust,
  up_color = "#810F7C",
  down_color = "#006D2C"
) {
  if (
    length(top) != 1L ||
      !is.numeric(top) ||
      is.na(top) ||
      !is.finite(top) ||
      top < 0 ||
      top != as.integer(top)
  ) {
    stop("top must be one non-negative integer.")
  }
  top <- as.integer(top)
  validate_color <- function(color, argument) {
    if (
      !is.character(color) ||
        length(color) != 1L ||
        is.na(color) ||
        !nzchar(color)
    ) {
      stop(argument, " must be one valid R color.")
    }
    valid <- tryCatch(
      {
        grDevices::col2rgb(color)
        TRUE
      },
      error = function(e) FALSE
    )
    if (!valid) {
      stop(argument, " must be one valid R color.")
    }
  }
  validate_color(up_color, "up_color")
  validate_color(down_color, "down_color")

  target_gene.list <- target_for_volcano(
    deseq2_result.df,
    gene.list,
    n = n,
    pc = pc,
    pv = pv,
    fc = fc,
    p_col = if (adjust) "padj" else "pvalue"
  )
  gg <- deseq2_result.df[, c("Symbol", "log2FoldChange")]
  if (adjust) {
    gg$padj <- deseq2_result.df$padj
  } else {
    gg$padj <- deseq2_result.df$pvalue
  }
  positive_pvalues <- gg$padj[is.finite(gg$padj) & gg$padj > 0]
  zero_pvalue_replacement <- if (length(positive_pvalues)) {
    min(positive_pvalues)
  } else {
    NA_real_
  }
  gg <- gg[match(target_gene.list, gg$Symbol), ]
  gg <- gg |> filter(!is.na(log2FoldChange), !is.na(padj))
  if (nrow(gg) == 0) {
    stop("No genes with valid fold changes and P values are available to plot.")
  }
  gg$group <- "no"
  up_index <- gg$log2FoldChange > log2(fc) & gg$padj < pv
  down_index <- gg$log2FoldChange < -log2(fc) & gg$padj < pv
  gg$group[up_index] <- "up"
  gg$group[down_index] <- "down"

  up_n <- sum(up_index)
  dw_n <- sum(down_index)
  count_message <- paste0(
    "Significant DEGs at P < ",
    pv,
    " and |FC| > ",
    fc,
    ": Up = ",
    up_n,
    ", Down = ",
    dw_n,
    "."
  )
  if (up_n == 0 && dw_n == 0) {
    message(count_message, " Plotting background genes only.")
  } else if (up_n == 0) {
    message(count_message, " No upregulated DEGs were detected.")
  } else if (dw_n == 0) {
    message(count_message, " No downregulated DEGs were detected.")
  }

  select_top_genes <- function(direction) {
    if (top == 0L) {
      return(character())
    }
    ranked <- gg[
      gg$group == direction & !is.na(gg$Symbol),
      ,
      drop = FALSE
    ]
    ranked <- ranked[
      order(
        ranked$padj,
        -abs(ranked$log2FoldChange),
        ranked$Symbol,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    ranked <- ranked[!duplicated(ranked$Symbol), , drop = FALSE]
    utils::head(ranked$Symbol, top)
  }
  gene.list <- c(
    select_top_genes("up"),
    select_top_genes("down"),
    gene.list
  )
  gene.list <- unique(gene.list)
  index <- match(gene.list, gg$Symbol)
  index <- na.omit(index)
  gg$color <- gg$group
  gg$color[index] <- "black"
  mycolour <- c("grey", up_color, down_color, "black")
  names(mycolour) <- c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
  # gg[gg$group == "no",]$label <- ""
  if (max(gg$log2FoldChange) > max_x) {
    gg[gg$log2FoldChange > max_x, ]$log2FoldChange <- max_x
  }
  if (min(gg$log2FoldChange) < -max_x) {
    gg[gg$log2FoldChange < -max_x, ]$log2FoldChange <- -max_x
  }
  if (any(gg$padj == 0)) {
    if (is.na(zero_pvalue_replacement)) {
      warning(
        "No finite positive P value was found in the full input ",
        if (adjust) "padj" else "pvalue",
        " column; zero P values will use the max_y cap.",
        call. = FALSE
      )
    } else {
      gg$padj[gg$padj == 0] <- zero_pvalue_replacement
    }
  }
  if (min(gg$padj) < 10^-max_y) {
    gg[gg$padj < 10^-max_y, ]$padj <- 10^-max_y
  }
  p <- ggplot(gg, aes(x = log2FoldChange, y = -log10(padj)))
  p <- p +
    geom_point_rast(
      aes(color = group),
      shape = 16,
      alpha = 0.6,
      show.legend = FALSE
    )
  #  p <- p + geom_point(color = ifelse(gg$Symbol %in% gene.list,"black", NA), shape = 1, show.legend = FALSE)
  p <- p +
    geom_point(
      color = ifelse(gg$label == "", NA, "black"),
      shape = 1,
      show.legend = FALSE
    )
  p <- p + scale_color_manual(values = mycolour)
  p <- p + scale_fill_manual(values = mycolour)
  p <- p + geom_hline(yintercept = -log10(pv), linetype = "dotted")
  p <- p + geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted")
  p <- p + labs(
    subtitle = paste0("Significant genes | Up: ", up_n, "  Down: ", dw_n)
  )
  p <- p +
    ggrepel::geom_text_repel(
      aes(label = label, color = group),
      show.legend = FALSE,
      fontface = "bold",
      size = 2.5,
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      max.overlaps = 10000
    )
  p <- p + theme_test()
  # p= p + xlim(-4,4) + ylim(0,100)
  if (adjust) {
    p <- p + ylab("-log10(padj)")
  } else {
    p <- p + ylab("-log10(pvalue)")
  }
  return(p)
}


#' Draw a volcano plot from single-cell marker results
#'
#' Classifies Seurat marker genes by adjusted P value and average log2 fold
#' change, caps extreme display values, and labels requested and top-ranked
#' genes.
#'
#' @param findmarkers.df Marker data frame containing `Symbol`, `avg_log2FC`,
#'   `p_val_adj`, and `pct.1` columns.
#' @param gene.list Character vector of gene symbols to label.
#' @param pv Adjusted P-value threshold.
#' @param fc Fold-change threshold on the linear scale.
#' @param top Number of top upregulated and downregulated genes to label.
#' @param max_x Maximum absolute average log2 fold change displayed.
#' @param max_y Maximum displayed `-log10(p_val_adj)` value.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
volcano_plot_scRNA <- function(
  findmarkers.df,
  gene.list,
  pv = 0.05,
  fc = 1.5,
  top = 5,
  max_x = 1,
  max_y = 50
) {
  #  gene.list <- c("BIRC3","CCR7", "NFKBIA", "TNFAIP3","REL","BCL3","BCL2","HSPA1A","SOD1","HSPB1", "PPP1R15A", "DNAJA1")
  #  fc <- 1.09
  #  pv <- 0.05
  #  max_x = 1
  #  max_y = 50
  #  gg <- CD4_Naive_T.MTvsWT.deg
  #  gg$Symbol <- as.character(rownames(CD4_Naive_T.MTvsWT.deg))
  gg <- findmarkers.df
  gg$Symbol <- as.character(rownames(findmarkers.df))
  top_up_genes <- gg |>
    filter(avg_log2FC > log2(fc), !is.na(p_val_adj)) |>
    slice_min(order_by = p_val_adj, n = top, with_ties = FALSE) |>
    pull(Symbol)
  top_down_genes <- gg |>
    filter(avg_log2FC < -log2(fc), !is.na(p_val_adj)) |>
    slice_min(order_by = p_val_adj, n = top, with_ties = FALSE) |>
    pull(Symbol)
  gene.list <- unique(c(top_up_genes, top_down_genes, gene.list))
  # gene.list <- c(rownames(gg %>% slice_head(n=8)),gene.list)
  index <- match(gene.list, gg$Symbol)
  index <- na.omit(index)
  gg$group <- "no"
  try(gg[gg$avg_log2FC > log2(fc) & gg$p_val_adj < pv, ]$group <- "up")
  try(gg[gg$avg_log2FC < -log2(fc) & gg$p_val_adj < pv, ]$group <- "down")
  # try(gg[gg$avg_log2FC > log2(fc) & gg$p_val < pv,]$group <- "up")
  # try(gg[gg$avg_log2FC < -log2(fc) & gg$p_val < pv,]$group <- "down")
  gg$color <- gg$group
  gg$color[index] <- "black"
  # mycolour = c("grey", "#B30000", "#08519C", "black")
  # mycolour = c("grey", "#810F7C", "#006D2C", "black")
  mycolour <- c("grey", "#A81E2C", "#08537C", "black")
  names(mycolour) <- c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
  #  gg[gg$group == "no",]$label <- ""
  #  if(max(gg$avg_log2FC) > max_x){gg[gg$avg_log2FC > max_x,]$avg_log2FC = max_x}
  #  if(min(gg$avg_log2FC) < -max_x){gg[gg$avg_log2FC < -max_x,]$avg_log2FC = -max_x}
  #  if(min(gg$p_val_adj) < 10^-max_y){gg[gg$p_val_adj < 10^-max_y,]$p_val_adj = 10^-max_y}
  p <- ggplot(gg, aes(x = avg_log2FC, y = -log10(p_val_adj)))
  # p <- ggplot(gg, aes(x = avg_log2FC, y = -log10(p_val)))
  p <- p +
    geom_point_rast(
      aes(fill = group, size = pct.1),
      shape = 21,
      alpha = 0.6,
      show.legend = TRUE
    )
  p <- p +
    geom_point_rast(
      aes(colour = color, size = pct.1),
      shape = 21,
      alpha = 0.6,
      show.legend = TRUE
    )
  p <- p + scale_color_manual(values = mycolour)
  p <- p + scale_fill_manual(values = mycolour)
  p <- p + geom_hline(yintercept = -log10(pv), linetype = "dotted")
  p <- p + geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted")
  p <- p +
    ggrepel::geom_text_repel(
      aes(label = label, color = group),
      show.legend = FALSE,
      fontface = "bold",
      size = 2.5,
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      max.overlaps = 10000
    )
  p <- p + blank
  p <- p +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    )
  #  p = p + xlim(-max_x,max_x) + ylim(0,max_y)
  return(p)
}


#' Compare two differential-expression analyses
#'
#' Merges two DEG result tables by gene identifier, compares caller-selected
#' numeric columns, classifies concordant and discordant patterns, optionally
#' samples background genes, and draws a correlation plot with selected labels.
#'
#' @param x_deg_result,y_deg_result DEG result data frames supplying the
#'   horizontal and vertical values, respectively. They must contain compatible
#'   merge identifiers, selected numeric columns, and P-value columns.
#' @param x_label,y_label Optional axis labels. When `NULL`, labels are generated
#'   from `x_col` and `y_col`.
#' @param fc Fold-change threshold on the linear scale, retained for backward
#'   compatibility. `log2(fc)` is used for an axis whose explicit threshold is
#'   `NULL`.
#' @param pvalue_cutoff Significance threshold applied to both selected P-value
#'   columns.
#' @param pc Logical; restrict `x_deg_result` to protein-coding genes. If
#'   `Gene_Type` is absent, it is added from the bundled gene annotation for
#'   `species`.
#' @param bg_num Maximum number of background genes sampled for display.
#' @param limit Default maximum absolute limit used for an axis whose explicit
#'   limits are `NULL`.
#' @param seed Random seed used for background sampling.
#' @param show_cor Logical; display correlation statistics.
#' @param plot_title Optional plot title.
#' @param goi Optional character vector of gene symbols to label. Matched genes
#'   are labelled even when they are not significant. Symbols absent from the
#'   matched, filtered plotting data are reported in a warning.
#' @param top Number of top genes labelled for each comparison group.
#' @param label_size Text size for gene labels.
#' @param x_col,y_col Numeric columns from `x_deg_result` and `y_deg_result`
#'   used for the horizontal and vertical axes, respectively.
#' @param x_pvalue_col,y_pvalue_col Numeric significance columns from
#'   `x_deg_result` and `y_deg_result`, respectively. These may name raw P-value,
#'   adjusted P-value, FDR, or q-value columns and may differ between tables.
#' @param x_threshold,y_threshold Non-negative absolute thresholds in the units
#'   of the selected axis columns. When `NULL`, `log2(fc)` is used.
#' @param x_limits,y_limits Optional numeric vectors of length two specifying
#'   independent lower and upper display limits. When `NULL`, `c(-limit, limit)`
#'   is used.
#' @param x_merge_col,y_merge_col Columns used to match rows from `x_deg_result`
#'   and `y_deg_result`, respectively. The columns may have different names but
#'   must contain compatible, non-duplicated identifiers.
#' @param species Species used to annotate `Gene_Type` when `pc = TRUE` and the
#'   column is absent from `x_deg_result`. Either `"human"` or `"mouse"`.
#'
#' @return A [ggplot2::ggplot] comparison plot.
#' @export
plot_deg_comparison <- function(
  x_deg_result,
  y_deg_result,
  x_label = NULL,
  y_label = NULL,
  fc = 1.5,
  pvalue_cutoff = 0.05,
  pc = TRUE,
  bg_num = 5000,
  limit = 5,
  seed = 614,
  show_cor = TRUE,
  plot_title = NULL,
  goi = NULL,
  top = 5,
  label_size = 2.5,
  x_col = "log2FoldChange",
  y_col = "log2FoldChange",
  x_pvalue_col = "padj",
  y_pvalue_col = "padj",
  x_threshold = NULL,
  y_threshold = NULL,
  x_limits = NULL,
  y_limits = NULL,
  x_merge_col = "Row.names",
  y_merge_col = "Row.names",
  species = c("human", "mouse")
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for plot_deg_comparison().")
  }

  .validate_column_name <- function(column, argument, data, data_name) {
    if (!is.character(column) || length(column) != 1L ||
      is.na(column) || !nzchar(column)) {
      stop(argument, " must be one non-empty column name.", call. = FALSE)
    }
    if (!column %in% colnames(data)) {
      stop(
        argument, " ('", column, "') was not found in ", data_name, ".",
        call. = FALSE
      )
    }
  }

  .validate_column_name(x_col, "x_col", x_deg_result, "x_deg_result")
  .validate_column_name(y_col, "y_col", y_deg_result, "y_deg_result")
  .validate_column_name(
    x_pvalue_col,
    "x_pvalue_col",
    x_deg_result,
    "x_deg_result"
  )
  .validate_column_name(
    y_pvalue_col,
    "y_pvalue_col",
    y_deg_result,
    "y_deg_result"
  )
  .validate_column_name(
    x_merge_col,
    "x_merge_col",
    x_deg_result,
    "x_deg_result"
  )
  .validate_column_name(
    y_merge_col,
    "y_merge_col",
    y_deg_result,
    "y_deg_result"
  )
  .validate_column_name("Symbol", "Symbol", x_deg_result, "x_deg_result")
  species <- .abel_normalize_species(species)
  if (!is.null(goi) && !is.character(goi)) {
    stop("goi must be NULL or a character vector of gene symbols.", call. = FALSE)
  }
  goi <- unique(goi[!is.na(goi) & nzchar(goi)])
  if (!is.numeric(x_deg_result[[x_col]])) {
    stop("x_col must select a numeric column in x_deg_result.", call. = FALSE)
  }
  if (!is.numeric(y_deg_result[[y_col]])) {
    stop("y_col must select a numeric column in y_deg_result.", call. = FALSE)
  }
  if (!is.numeric(x_deg_result[[x_pvalue_col]]) ||
    !is.numeric(y_deg_result[[y_pvalue_col]])) {
    stop("The selected P-value columns must be numeric.", call. = FALSE)
  }
  if (!is.numeric(pvalue_cutoff) || length(pvalue_cutoff) != 1L ||
    !is.finite(pvalue_cutoff) || pvalue_cutoff < 0 || pvalue_cutoff > 1) {
    stop("pvalue_cutoff must be one number between 0 and 1.", call. = FALSE)
  }

  fallback_threshold <- NULL
  if (is.null(x_threshold) || is.null(y_threshold)) {
    if (!is.numeric(fc) || length(fc) != 1L || !is.finite(fc) || fc < 1) {
      stop("fc must be one finite number greater than or equal to 1.", call. = FALSE)
    }
    fallback_threshold <- log2(fc)
  }
  x_threshold <- .validate_deg_axis_threshold(
    x_threshold,
    fallback_threshold,
    "x_threshold"
  )
  y_threshold <- .validate_deg_axis_threshold(
    y_threshold,
    fallback_threshold,
    "y_threshold"
  )
  x_limits <- .validate_deg_axis_limits(x_limits, limit, "x_limits")
  y_limits <- .validate_deg_axis_limits(y_limits, limit, "y_limits")

  if (is.null(x_label)) {
    x_label <- paste0(x_col, " in x_deg_result")
  }
  if (is.null(y_label)) {
    y_label <- paste0(y_col, " in y_deg_result")
  }

  x_merge_key <- as.character(x_deg_result[[x_merge_col]])
  y_merge_key <- as.character(y_deg_result[[y_merge_col]])
  valid_x_key <- !is.na(x_merge_key) & nzchar(x_merge_key)
  valid_y_key <- !is.na(y_merge_key) & nzchar(y_merge_key)
  if (anyDuplicated(x_merge_key[valid_x_key])) {
    stop(
      "x_merge_col ('", x_merge_col,
      "') contains duplicated non-missing identifiers.",
      call. = FALSE
    )
  }
  if (anyDuplicated(y_merge_key[valid_y_key])) {
    stop(
      "y_merge_col ('", y_merge_col,
      "') contains duplicated non-missing identifiers.",
      call. = FALSE
    )
  }

  gene_type <- NULL
  if ("Gene_Type" %in% colnames(x_deg_result)) {
    gene_type <- x_deg_result[["Gene_Type"]]
  } else if (pc) {
    gene_annotation <- .abel_gene_annotation(species)
    if (!"Gene_Type" %in% colnames(gene_annotation)) {
      stop(
        "The bundled ", species,
        " gene annotation does not contain a Gene_Type column.",
        call. = FALSE
      )
    }

    annotation_ids <- rownames(gene_annotation)
    normalized_x_ids <- sub("\\..*$", "", x_merge_key)
    gene_type <- gene_annotation[["Gene_Type"]][
      match(normalized_x_ids, annotation_ids)
    ]

    if ("Symbol" %in% colnames(gene_annotation)) {
      annotation_symbols <- as.character(gene_annotation[["Symbol"]])
      unique_symbol <- !is.na(annotation_symbols) &
        nzchar(annotation_symbols) &
        !duplicated(annotation_symbols) &
        !duplicated(annotation_symbols, fromLast = TRUE)
      symbol_match <- match(
        as.character(x_deg_result[["Symbol"]]),
        annotation_symbols[unique_symbol]
      )
      missing_gene_type <- is.na(gene_type)
      gene_type[missing_gene_type] <- gene_annotation[["Gene_Type"]][
        which(unique_symbol)[symbol_match[missing_gene_type]]
      ]
    }

    if (all(is.na(gene_type))) {
      stop(
        "Gene_Type could not be annotated for x_deg_result using the bundled ",
        species, " gene annotation. Check species, x_merge_col, and Symbol.",
        call. = FALSE
      )
    }
  }

  x_degs <- data.frame(
    Merge_ID = x_merge_key,
    Symbol = x_deg_result[["Symbol"]],
    X_value = x_deg_result[[x_col]],
    X_pvalue = x_deg_result[[x_pvalue_col]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(gene_type)) {
    x_degs$Gene_Type <- gene_type
  }
  y_degs <- data.frame(
    Merge_ID = y_merge_key,
    Y_value = y_deg_result[[y_col]],
    Y_pvalue = y_deg_result[[y_pvalue_col]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  x_degs <- x_degs[valid_x_key, , drop = FALSE]
  y_degs <- y_degs[valid_y_key, , drop = FALSE]

  vs.degs <- inner_join(x_degs, y_degs, by = "Merge_ID")
  if (pc) {
    vs.degs <- vs.degs |> filter(Gene_Type == "protein_coding")
  }
  vs.degs <- vs.degs |>
    filter(
      is.finite(X_value),
      is.finite(Y_value),
      is.finite(X_pvalue),
      is.finite(Y_pvalue),
      !is.na(Symbol),
      nzchar(Symbol)
    )
  if (nrow(vs.degs) < 3L) {
    stop(
      "At least three complete matched genes are required for comparison.",
      call. = FALSE
    )
  }
  if (stats::sd(vs.degs$X_value) == 0 || stats::sd(vs.degs$Y_value) == 0) {
    stop("The selected axis columns must contain variation.", call. = FALSE)
  }

  goi_found <- intersect(goi, unique(vs.degs$Symbol))
  goi_missing <- setdiff(goi, goi_found)
  if (length(goi_missing) > 0L) {
    warning(
      "The following goi gene symbols were not found in the matched, filtered ",
      "plotting data and cannot be labelled: ",
      paste(goi_missing, collapse = ", "),
      call. = FALSE
    )
  }

  classify_with_pval <- function(x, y, x_pvalue, y_pvalue) {
    case_when(
      x > x_threshold & y > y_threshold &
        x_pvalue < pvalue_cutoff & y_pvalue < pvalue_cutoff ~ "upup",
      x < -x_threshold & y < -y_threshold &
        x_pvalue < pvalue_cutoff & y_pvalue < pvalue_cutoff ~ "dwdw",
      x > x_threshold & y < -y_threshold &
        x_pvalue < pvalue_cutoff & y_pvalue < pvalue_cutoff ~ "updw",
      x < -x_threshold & y > y_threshold &
        x_pvalue < pvalue_cutoff & y_pvalue < pvalue_cutoff ~ "dwup",
      TRUE ~ "other"
    )
  }

  classify_values_only <- function(x, y) {
    case_when(
      x > x_threshold & y > y_threshold ~ "upup",
      x < -x_threshold & y < -y_threshold ~ "dwdw",
      x > x_threshold & y < -y_threshold ~ "updw",
      x < -x_threshold & y > y_threshold ~ "dwup",
      TRUE ~ "other"
    )
  }

  vs.degs <- vs.degs |>
    mutate(
      Group = classify_with_pval(X_value, Y_value, X_pvalue, Y_pvalue),
      Group2 = classify_values_only(X_value, Y_value)
    )

  x_deg_list <- vs.degs |>
    filter(abs(X_value) > x_threshold, X_pvalue < pvalue_cutoff) |>
    pull(Symbol)

  y_deg_list <- vs.degs |>
    filter(abs(Y_value) > y_threshold, Y_pvalue < pvalue_cutoff) |>
    pull(Symbol)

  deg.list <- unique(c(x_deg_list, y_deg_list))
  other.list <- vs.degs |>
    filter(!Symbol %in% deg.list) |>
    pull(Symbol)

  set.seed(seed)
  random.list <- sample(other.list, min(bg_num, length(other.list)))
  gene.list <- unique(c(goi_found, deg.list, random.list))
  vs.degs <- vs.degs |> filter(Symbol %in% gene.list)
  if (nrow(vs.degs) < 3L ||
    stats::sd(vs.degs$X_value) == 0 ||
    stats::sd(vs.degs$Y_value) == 0) {
    stop(
      "At least three plotted genes with variation on both axes are required; ",
      "increase bg_num or revise the filters.",
      call. = FALSE
    )
  }
  goi_final <- goi_found

  top_upup <- character()
  top_dwdw <- character()
  if (!is.null(top) && top > 0) {
    top_upup <- vs.degs |>
      filter(Group == "upup") |>
      slice_max(order_by = X_value * Y_value, n = top, with_ties = FALSE) |>
      pull(Symbol)
    top_dwdw <- vs.degs |>
      filter(Group == "dwdw") |>
      slice_max(order_by = X_value * Y_value, n = top, with_ties = FALSE) |>
      pull(Symbol)
  }

  goi_final <- unique(c(goi_final, top_upup, top_dwdw))
  vs.degs <- vs.degs |>
    mutate(GOI = if_else(Symbol %in% goi_final, Symbol, ""))

  cor_result <- cor.test(vs.degs$X_value, vs.degs$Y_value)

  if (show_cor) {
    cat("\nCorrelation test result:\n")
    print(cor_result)
  }

  lm_model <- lm(Y_value ~ X_value, data = vs.degs)
  slope_k <- lm_model$coefficients[2]
  cor_r <- cor_result$estimate
  cor_p <- cor_result$p.value
  sig_stars <- case_when(
    cor_p < 0.001 ~ "***",
    cor_p < 0.01 ~ "**",
    cor_p < 0.05 ~ "*",
    TRUE ~ "ns"
  )

  cor_label <- paste0(
    "italic(K)==",
    round(slope_k, 4),
    "*','~italic(r)==",
    round(cor_r, 4),
    "*' ",
    sig_stars,
    "'"
  )

  mycolour <- c(
    "dwdw" = "#3C5488",
    "upup" = "#A81E2C",
    "updw" = "#B5AA0F",
    "dwup" = "#7BA39D",
    "other" = "grey"
  )

  vs.degs <- vs.degs |>
    mutate(
      Plot_X = pmin(pmax(X_value, x_limits[1]), x_limits[2]),
      Plot_Y = pmin(pmax(Y_value, y_limits[1]), y_limits[2])
    )

  p <- vs.degs |>
    ggplot(aes(Plot_X, Plot_Y)) +
    geom_point_rast(
      aes(color = Group2),
      shape = 16,
      alpha = 0.6,
      show.legend = FALSE
    ) +
    geom_point_rast(
      data = vs.degs |> filter(Group != "other"),
      aes(color = Group),
      shape = 1,
      alpha = 0.6,
      show.legend = FALSE
    ) +
    geom_point(
      data = vs.degs |> filter(GOI != ""),
      color = "black",
      shape = 1,
      show.legend = FALSE
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_test() +
    xlab(x_label) +
    ylab(y_label) +
    ggrepel::geom_text_repel(
      data = vs.degs |> filter(GOI != ""),
      aes(label = GOI, color = Group2),
      show.legend = FALSE,
      fontface = "bold",
      size = label_size,
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      max.overlaps = 10000
    ) +
    geom_hline(yintercept = c(-y_threshold, y_threshold), linetype = "dotted") +
    geom_vline(xintercept = c(-x_threshold, x_threshold), linetype = "dotted") +
    annotate(
      "text",
      x = x_limits[1] + diff(x_limits) * 0.2,
      y = y_limits[2] - diff(y_limits) * 0.1,
      label = cor_label,
      size = 4,
      hjust = 0.5,
      parse = TRUE
    ) +
    scale_color_manual(values = mycolour)

  if (!is.null(plot_title)) {
    p <- p +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  return(p)
}


.validate_deg_axis_threshold <- function(value, fallback, argument) {
  if (is.null(value)) {
    return(fallback)
  }
  if (!is.numeric(value) || length(value) != 1L ||
    !is.finite(value) || value < 0) {
    stop(argument, " must be one finite non-negative number.", call. = FALSE)
  }
  unname(value)
}


.validate_deg_axis_limits <- function(value, limit, argument) {
  if (is.null(value)) {
    if (!is.numeric(limit) || length(limit) != 1L ||
      !is.finite(limit) || limit <= 0) {
      stop("limit must be one finite positive number.", call. = FALSE)
    }
    return(c(-limit, limit))
  }
  if (!is.numeric(value) || length(value) != 2L ||
    any(!is.finite(value))) {
    stop(argument, " must contain two finite numeric values.", call. = FALSE)
  }
  if (value[1] >= value[2]) {
    stop(argument, " must be ordered from lower to upper.", call. = FALSE)
  }
  unname(value)
}


#' Draw genomic-position DEG Manhattan plots
#'
#' Joins one or more DEG tables to human or mouse gene coordinates, calculates
#' a signed significance score (`log2FoldChange * -log10(padj)`), and draws
#' faceted genomic-position plots with separate top upregulated and
#' downregulated gene labels.
#'
#' @param deg_list Named list of DEG data frames, one per treatment or contrast.
#' @param color_map Named colour vector whose names match `deg_list`.
#' @param species Species defining bundled annotation and default chromosomes.
#' @param deg_cols Optional four-column character vector overriding
#'   `symbol_col`, `lfc_col`, `p_col`, and `padj_col`.
#' @param symbol_col,lfc_col,p_col,padj_col Column names containing gene symbol,
#'   log2 fold change, raw P value, and adjusted P value.
#' @param gene_anno_file Optional compatible gene-annotation table. Bundled
#'   annotation is used when `NULL`.
#' @param chromosome_lengths Optional named numeric vector of chromosome lengths.
#'   Annotation maxima are used when `NULL`.
#' @param gene_type_filter Optional gene type retained when `Gene_Type` exists.
#' @param remove_rik Logical; remove mouse symbols ending in `Rik`.
#' @param top_n Number of upregulated and downregulated genes labelled per
#'   treatment. Up to `top_n` genes are selected independently in each
#'   direction.
#' @param cap_value Maximum absolute signed significance score displayed.
#' @param chr_keep Optional chromosomes and order to retain.
#' @param facet_nrow Number of rows in the treatment facet layout.
#'
#' @return A list with the joined data (`deg_merge`), labelled genes
#'   (`top_genes`), and the [ggplot2::ggplot] object (`plot`).
#' @export
plot_deg_manhattan <- function(
  deg_list,
  color_map,
  species = c("mouse", "human"),
  deg_cols = NULL,
  symbol_col = "Symbol",
  lfc_col = "log2FoldChange",
  p_col = "pvalue",
  padj_col = "padj",
  gene_anno_file = NULL,
  chromosome_lengths = NULL,
  gene_type_filter = "protein_coding",
  remove_rik = TRUE,
  top_n = 10,
  cap_value = 20,
  chr_keep = NULL,
  facet_nrow = 1
) {
  species <- match.arg(species)

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for plot_deg_manhattan().")
  }
  if (
    length(top_n) != 1L ||
      !is.numeric(top_n) ||
      is.na(top_n) ||
      !is.finite(top_n) ||
      top_n < 0 ||
      top_n != as.integer(top_n)
  ) {
    stop("top_n must be one non-negative integer.")
  }
  top_n <- as.integer(top_n)

  if (species == "mouse") {
    if (is.null(chr_keep)) {
      chr_keep <- c(paste0("chr", 1:19), "chrX", "chrY")
    }
  } else if (species == "human") {
    if (is.null(chr_keep)) {
      chr_keep <- c(paste0("chr", 1:22), "chrX", "chrY")
    }
  }

  if (!is.null(deg_cols)) {
    if (length(deg_cols) != 4) {
      stop("deg_cols must contain exactly four column names.")
    }
    symbol_col <- deg_cols[1]
    lfc_col <- deg_cols[2]
    p_col <- deg_cols[3]
    padj_col <- deg_cols[4]
  }
  deg_cols <- c(symbol_col, lfc_col, p_col, padj_col)

  gene_anno <- .abel_gene_annotation(species, gene_anno_file)

  required_anno_cols <- c("Symbol", "Chr", "Start")
  missing_anno <- setdiff(required_anno_cols, colnames(gene_anno))
  if (length(missing_anno) > 0) {
    stop("gene_anno_file is missing required columns: ", paste(missing_anno, collapse = ", "))
  }

  if (is.null(chromosome_lengths)) {
    end_col <- if ("End" %in% colnames(gene_anno)) "End" else "Start"
    chromosome_lengths <- tapply(
      gene_anno[[end_col]],
      gene_anno$Chr,
      max,
      na.rm = TRUE
    )
  }
  chromosome_lengths <- chromosome_lengths[chr_keep]
  chromosome_lengths <- chromosome_lengths[!is.na(chromosome_lengths)]
  if (!length(chromosome_lengths)) {
    stop("No chromosome lengths are available for chr_keep.")
  }
  offset <- c(0, head(cumsum(as.numeric(chromosome_lengths)), -1))
  names(offset) <- names(chromosome_lengths)

  gene_anno <- gene_anno %>%
    mutate(
      Chr = as.character(Chr)
    ) %>%
    filter(Chr %in% names(offset)) %>%
    mutate(
      offset = offset[Chr],
      start_g = offset + Start
    ) %>%
    select(-offset)

  ## merge DEG tables
  tmp_degs_list <- lapply(names(deg_list), function(trt) {
    if (!trt %in% names(color_map)) {
      stop("Treatment '", trt, "' not found in color_map")
    }

    deg_df <- deg_list[[trt]]

    missing_deg_cols <- setdiff(deg_cols, colnames(deg_df))
    if (length(missing_deg_cols) > 0) {
      stop(
        "In treatment '", trt, "', DEG table is missing columns: ",
        paste(missing_deg_cols, collapse = ", ")
      )
    }

    deg_df <- deg_df[, deg_cols, drop = FALSE]
    colnames(deg_df) <- c("Symbol", "log2FoldChange", "pvalue", "padj")

    deg_df %>%
      left_join(gene_anno, by = "Symbol") %>%
      mutate(Treatment = trt)
  })

  deg_merge <- bind_rows(tmp_degs_list)

  ## score
  deg_merge <- deg_merge %>%
    mutate(
      padj_safe = dplyr::if_else(is.na(padj), NA_real_, pmax(padj, 1e-300)),
      value = log2FoldChange * -log10(padj_safe)
    ) %>%
    filter(!is.na(value), !is.na(start_g))

  ## optional filters
  if ("Gene_Type" %in% colnames(deg_merge) && !is.null(gene_type_filter)) {
    deg_merge <- deg_merge %>% filter(Gene_Type == gene_type_filter)
  }

  if (remove_rik) {
    deg_merge <- deg_merge %>% filter(!grepl("Rik$", Symbol))
  }

  deg_merge <- deg_merge %>%
    mutate(
      Regulation = case_when(
        log2FoldChange > 0 & padj < 0.05 ~ "Up",
        log2FoldChange < 0 & padj < 0.05 ~ "Down",
        TRUE ~ "NS"
      ),
      Regulation = factor(Regulation, levels = c("Up", "Down", "NS"))
    )

  ## select top upregulated and downregulated genes independently
  top_genes <- deg_merge %>%
    filter(Regulation %in% c("Up", "Down")) %>%
    arrange(
      Treatment,
      Regulation,
      desc(abs(value)),
      padj,
      desc(abs(log2FoldChange)),
      Symbol
    ) %>%
    distinct(Treatment, Regulation, Symbol, .keep_all = TRUE) %>%
    group_by(Treatment, Regulation) %>%
    slice_head(n = top_n) %>%
    ungroup()

  ## cap display values after ranking
  deg_merge <- deg_merge %>%
    mutate(
      value = ifelse(value > cap_value, cap_value, value),
      value = ifelse(value < -cap_value, -cap_value, value)
    )

  top_keys <- paste(top_genes$Symbol, top_genes$Treatment, sep = "\r")
  deg_merge <- deg_merge %>%
    mutate(
      label = if_else(
        paste(Symbol, Treatment, sep = "\r") %in% top_keys,
        Symbol,
        ""
      ),
      Treatment = factor(Treatment, levels = names(color_map))
    )

  reg_color_map <- c(
    Up = "#A81E2C",
    Down = "#3C5488",
    NS = "grey80"
  )

  ## plot
  p_man <- ggplot(deg_merge, aes(x = start_g, y = value)) +
    ggrastr::geom_point_rast(
      aes(color = Treatment),
      shape = 16, size = 0.7, alpha = 0.6
    ) +
    geom_point(
      data = deg_merge %>% filter(label != ""),
      color = "black", shape = 1, size = 0.8, show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = deg_merge %>% filter(label != ""),
      aes(label = label, color = Regulation),
      show.legend = FALSE,
      fontface = "bold",
      size = 2.5,
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      max.overlaps = Inf
    ) +
    facet_wrap(~Treatment, nrow = facet_nrow) +
    scale_color_manual(values = c(color_map, reg_color_map)) +
    theme_test() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    xlab("") +
    ylab("log2FoldChange * -log10(padj)")

  return(list(
    deg_merge = deg_merge,
    top_genes = top_genes,
    plot = p_man
  ))
}
