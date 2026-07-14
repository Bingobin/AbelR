# Functions for deg-plots analyses.

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
  if (length(deg_genes) < 2) {
    message(sample_name, ": fewer than 2 DEGs, skip heatmap.")
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


target_for_volcano <- function(
  deseq2_result.df,
  gene.list,
  n = 5000,
  pc = FALSE,
  pv = 0.05,
  fc = 1.5
) {
  #deseq2_result.df <- RELA_MUTvsWT_NON.DEGs$result
  if (pc == TRUE) {
    deseq2_result.df <- deseq2_result.df |>
      filter(Gene_Type == "protein_coding")
  }
  deg.list <- (deseq2_result.df |>
    filter(abs(log2FoldChange) > log2(fc), padj < pv))$Symbol
  other.list <- (deseq2_result.df |>
    filter(!Symbol %in% deg.list, !is.na(log2FoldChange), !is.na(padj)))$Symbol
  set.seed(123)
  random.list <- sample(other.list, n)
  result.list <- c(deg.list, random.list)
  result.list <- c(result.list, gene.list)
  return(result.list)
}


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
  adjust
) {
  #  deseq2_result.df <- WT1_MUT.deg$result
  #  target_gene.list <- WT1_status.gene$Symbol
  #  gene.list <- risk_gene.list
  #  pv <- 0.05
  #  fc <- 1.5
  target_gene.list <- target_for_volcano(
    deseq2_result.df,
    gene.list,
    n = n,
    pc = pc,
    pv = pv,
    fc = fc
  )
  gg <- deseq2_result.df[, c("Symbol", "log2FoldChange")]
  if (adjust) {
    gg$padj <- deseq2_result.df$padj
  } else {
    gg$padj <- deseq2_result.df$pvalue
  }
  gg <- gg[match(target_gene.list, gg$Symbol), ]
  gg <- gg |> filter(!is.na(log2FoldChange), !is.na(padj))
  gene.list <- c(
    (gg |> filter(log2FoldChange > 0) |> top_n(n = -top, wt = padj))$Symbol,
    gene.list
  )
  gene.list <- c(
    (gg |> filter(log2FoldChange < 0) |> top_n(n = -top, wt = padj))$Symbol,
    gene.list
  )
  index <- match(gene.list, gg$Symbol)
  index <- na.omit(index)
  gg$group <- "no"
  gg[gg$log2FoldChange > log2(fc) & gg$padj < pv, ]$group <- "up"
  gg[gg$log2FoldChange < -log2(fc) & gg$padj < pv, ]$group <- "down"
  gg$color <- gg$group
  gg$color[index] <- "black"
  #mycolour = c("grey", "#A81E2C","#08537C", "black")
  mycolour <- c("grey", "#810F7C", "#006D2C", "black")
  names(mycolour) <- c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
  #gg[gg$group == "no",]$label <- ""
  if (max(gg$log2FoldChange) > max_x) {
    gg[gg$log2FoldChange > max_x, ]$log2FoldChange <- max_x
  }
  if (min(gg$log2FoldChange) < -max_x) {
    gg[gg$log2FoldChange < -max_x, ]$log2FoldChange <- -max_x
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
  #p= p + xlim(-4,4) + ylim(0,100)
  if (adjust) {
    p <- p + ylab("-log10(padj)")
  } else {
    p <- p + ylab("-log10(pvalue)")
  }
  return(p)
}


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
  #gene.list <- c(rownames(gg %>% slice_head(n=8)),gene.list)
  index <- match(gene.list, gg$Symbol)
  index <- na.omit(index)
  gg$group <- "no"
  try(gg[gg$avg_log2FC > log2(fc) & gg$p_val_adj < pv, ]$group <- "up")
  try(gg[gg$avg_log2FC < -log2(fc) & gg$p_val_adj < pv, ]$group <- "down")
  #try(gg[gg$avg_log2FC > log2(fc) & gg$p_val < pv,]$group <- "up")
  #try(gg[gg$avg_log2FC < -log2(fc) & gg$p_val < pv,]$group <- "down")
  gg$color <- gg$group
  gg$color[index] <- "black"
  #mycolour = c("grey", "#B30000", "#08519C", "black")
  #mycolour = c("grey", "#810F7C", "#006D2C", "black")
  mycolour <- c("grey", "#A81E2C", "#08537C", "black")
  names(mycolour) <- c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
  #  gg[gg$group == "no",]$label <- ""
  #  if(max(gg$avg_log2FC) > max_x){gg[gg$avg_log2FC > max_x,]$avg_log2FC = max_x}
  #  if(min(gg$avg_log2FC) < -max_x){gg[gg$avg_log2FC < -max_x,]$avg_log2FC = -max_x}
  #  if(min(gg$p_val_adj) < 10^-max_y){gg[gg$p_val_adj < 10^-max_y,]$p_val_adj = 10^-max_y}
  p <- ggplot(gg, aes(x = avg_log2FC, y = -log10(p_val_adj)))
  #p <- ggplot(gg, aes(x = avg_log2FC, y = -log10(p_val)))
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


plot_deg_comparison <- function(
  deg_result1,
  deg_result2,
  label_1 = "Log2FoldChange in Sample1",
  label_2 = "Log2FoldChange in Sample2",
  fc = 1.5,
  pv = 0.05,
  adjust = TRUE,
  pc = TRUE,
  bg_num = 5000,
  limit = 5,
  seed = 614,
  show_cor = TRUE,
  plot_title = NULL,
  goi = NULL,
  top = 5,
  label_size = 2.5
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for plot_deg_comparison().")
  }

  # 选择P值列名
  pval_col <- if (adjust) "padj" else "pvalue"

  # 提取并重命名列
  degs_1 <- deg_result1 |>
    select(Row.names, Symbol, Gene_Type, log2FoldChange, all_of(pval_col)) |>
    rename(LFC_1 = log2FoldChange, PV_1 = all_of(pval_col))

  degs_2 <- deg_result2 |>
    select(Row.names, log2FoldChange, all_of(pval_col)) |>
    rename(LFC_2 = log2FoldChange, PV_2 = all_of(pval_col))

  # 合并数据
  vs.degs <- degs_1 |>
    inner_join(degs_2, by = "Row.names")

  # 计算log2(fc)阈值
  log2_fc <- log2(fc)

  # 分类函数(考虑P值)
  classify_with_pval <- function(lfc1, lfc2, pv1, pv2) {
    case_when(
      lfc1 > log2_fc & lfc2 > log2_fc & pv1 < pv & pv2 < pv ~ "upup",
      lfc1 < -log2_fc & lfc2 < -log2_fc & pv1 < pv & pv2 < pv ~ "dwdw",
      lfc1 > log2_fc & lfc2 < -log2_fc & pv1 < pv & pv2 < pv ~ "updw",
      lfc1 < -log2_fc & lfc2 > log2_fc & pv1 < pv & pv2 < pv ~ "dwup",
      TRUE ~ "other"
    )
  }

  # 分类函数(仅考虑FC)
  classify_fc_only <- function(lfc1, lfc2) {
    case_when(
      lfc1 > log2_fc & lfc2 > log2_fc ~ "upup",
      lfc1 < -log2_fc & lfc2 < -log2_fc ~ "dwdw",
      lfc1 > log2_fc & lfc2 < -log2_fc ~ "updw",
      lfc1 < -log2_fc & lfc2 > log2_fc ~ "dwup",
      TRUE ~ "other"
    )
  }

  # 添加分组
  vs.degs <- vs.degs |>
    mutate(
      Group = classify_with_pval(LFC_1, LFC_2, PV_1, PV_2),
      Group2 = classify_fc_only(LFC_1, LFC_2)
    )

  # 筛选protein_coding基因
  if (pc) {
    vs.degs <- vs.degs |> filter(Gene_Type == "protein_coding")
  }

  # 获取显著差异基因列表
  deg.list_1 <- vs.degs |>
    filter(abs(LFC_1) > log2_fc, PV_1 < pv) |>
    pull(Symbol)

  deg.list_2 <- vs.degs |>
    filter(abs(LFC_2) > log2_fc, PV_2 < pv) |>
    pull(Symbol)

  deg.list <- unique(c(deg.list_1, deg.list_2))

  # 获取背景基因并随机抽样
  other.list <- vs.degs |>
    filter(
      !Symbol %in% deg.list,
      !is.na(LFC_1),
      !is.na(PV_1),
      !is.na(LFC_2),
      !is.na(PV_2)
    ) |>
    pull(Symbol)

  set.seed(seed)
  random.list <- sample(other.list, min(bg_num, length(other.list)))
  gene.list <- c(deg.list, random.list)
  gene.list <- unique(c(gene.list, goi))

  # 筛选用于绘图的基因
  vs.degs <- vs.degs |> filter(Symbol %in% gene.list)

  # 选择要标注的基因
  goi_final <- goi

  # 自动添加top基因(upup组)
  top_upup <- vs.degs |>
    filter(Group == "upup") |>
    top_n(n = top, wt = LFC_1 * LFC_2) |>
    pull(Symbol)

  # 自动添加top基因(dwdw组)
  top_dwdw <- vs.degs |>
    filter(Group == "dwdw") |>
    top_n(n = top, wt = LFC_1 * LFC_2) |>
    pull(Symbol)

  goi_final <- unique(c(goi_final, top_upup, top_dwdw))

  # 添加GOI标记列
  vs.degs <- vs.degs |>
    mutate(GOI = if_else(Symbol %in% goi_final, Symbol, ""))

  # 相关性检验
  cor_result <- cor.test(vs.degs$LFC_1, vs.degs$LFC_2)

  if (show_cor) {
    cat("\nCorrelation test result:\n")
    print(cor_result)
  }

  # 线性回归计算斜率
  lm_model <- lm(LFC_2 ~ LFC_1, data = vs.degs)
  slope_k <- lm_model$coefficients[2]

  # 提取相关系数和p值
  cor_r <- cor_result$estimate
  cor_p <- cor_result$p.value

  # 确定显著性星号
  sig_stars <- case_when(
    cor_p < 0.001 ~ "***",
    cor_p < 0.01 ~ "**",
    cor_p < 0.05 ~ "*",
    TRUE ~ "ns"
  )

  # 创建标注文本(使用斜体K和r)
  cor_label <- paste0(
    "italic(K)==",
    round(slope_k, 4),
    "*','~italic(r)==",
    round(cor_r, 4),
    "*' ",
    sig_stars,
    "'"
  )

  # 定义颜色
  mycolour <- c(
    "dwdw" = "#3C5488",
    "upup" = "#A81E2C",
    "updw" = "#B5AA0F",
    "dwup" = "#7BA39D",
    "other" = "grey"
  )

  # 将超出limit的值裁剪到limit范围内
  vs.degs <- vs.degs |>
    mutate(
      LFC_1 = case_when(
        LFC_1 > limit ~ limit,
        LFC_1 < -limit ~ -limit,
        TRUE ~ LFC_1
      ),
      LFC_2 = case_when(
        LFC_2 > limit ~ limit,
        LFC_2 < -limit ~ -limit,
        TRUE ~ LFC_2
      )
    )

  # 绘图
  p <- vs.degs |>
    ggplot(aes(LFC_1, LFC_2)) +
    # 第一层:按Group2着色(仅FC标准)
    geom_point_rast(
      aes(color = Group2),
      shape = 16,
      alpha = 0.6,
      show.legend = FALSE
    ) +
    # 第二层:按Group着色(FC+PV标准)的空心点
    #geom_point_rast(aes(color = Group), shape = 1, alpha = 0.6, show.legend = FALSE) +
    geom_point_rast(
      data = vs.degs |> filter(Group != "other"),
      aes(color = Group),
      shape = 1,
      alpha = 0.6,
      show.legend = FALSE
    ) +
    # 第三层:标注感兴趣基因的黑色空心圆
    geom_point(
      color = ifelse(vs.degs$GOI == "", NA, "black"),
      shape = 1,
      show.legend = FALSE
    ) +
    # 坐标轴范围
    xlim(-limit, limit) +
    ylim(-limit, limit) +
    # 主题和标签
    theme_test() +
    xlab(label_1) +
    ylab(label_2) +
    # 基因名标注
    ggrepel::geom_text_repel(
      aes(label = GOI, color = Group2),
      show.legend = FALSE,
      fontface = "bold",
      size = label_size,
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      max.overlaps = 10000
    ) +
    # 阈值线
    geom_hline(yintercept = c(-log2_fc, log2_fc), linetype = "dotted") +
    geom_vline(xintercept = c(-log2_fc, log2_fc), linetype = "dotted") +
    # 相关性标注(左上角)
    annotate(
      "text",
      x = -limit * 0.6,
      y = limit * 0.9,
      label = cor_label,
      size = 4,
      hjust = 0.5,
      parse = TRUE
    ) +
    # 颜色设置
    scale_color_manual(values = mycolour)

  # 添加标题(如果提供)
  if (!is.null(plot_title)) {
    p <- p +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  return(p)
}


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
      stop("In treatment '", trt, "', DEG table is missing columns: ",
           paste(missing_deg_cols, collapse = ", "))
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
      padj_safe = dplyr::if_else(is.na(padj) , NA_real_,  pmax(padj, 1e-300)),
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
  
  ## cap values
  deg_merge <- deg_merge %>%
    mutate(
      value = ifelse(value > cap_value, cap_value, value),
      value = ifelse(value < -cap_value, -cap_value, value)
    )
  
  ## top genes
  top_genes <- deg_merge %>%
    group_by(Treatment) %>%
    slice_max(order_by = abs(value), n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  deg_merge <- deg_merge %>%
    mutate(
      label = if_else(
        paste0(Symbol, Treatment) %in% paste0(top_genes$Symbol, top_genes$Treatment),
        Symbol,
        ""
      ),
      Treatment = factor(Treatment, levels = names(color_map))
    )
  
  deg_merge <- deg_merge %>%
    mutate(
      Regulation = case_when(
        log2FoldChange > 0 & padj < 0.05 ~ "Up",
        log2FoldChange < 0 & padj < 0.05 ~ "Down",
        TRUE ~ "NS"
      ),
      Regulation = factor(Regulation, levels = c("Up", "Down", "NS"))
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
    facet_wrap(~ Treatment, nrow = facet_nrow) +
    scale_color_manual(values = c(color_map,reg_color_map)) +
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
