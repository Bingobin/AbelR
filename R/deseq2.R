# Functions for deseq2 analyses.

DESeq2_DEG_analysis_batch <- function(
  count.matrix,
  tpm.matrix = NULL,
  design.df,
  library = NULL,
  batch = NULL,
  tr = NULL,
  ctr = NULL,
  tr_filter = NULL,
  ctr_filter = NULL,
  tr_name = NULL,
  ctr_name = NULL,
  species = c("human", "mouse"),
  group_col = "Group",
  library_col = "Library",
  batch_col = "Batch",
  gene_anno_file = NULL,
  gene_id_col = "GID"
) {
  species <- tolower(species[1])
  species <- switch(
    species,
    hsa = "human",
    "homo sapiens" = "human",
    mmu = "mouse",
    "mus musculus" = "mouse",
    species
  )
  species <- match.arg(species, c("human", "mouse"))
  organism <- if (species == "human") "hsa" else "mmu"
  option_name <- paste0("AbelR.", species, "_gene_anno_file")
  if (is.null(gene_anno_file)) {
    gene_anno_file <- getOption(option_name)
  }
  if (is.null(gene_anno_file) || !file.exists(gene_anno_file)) {
    stop(
      "A valid gene annotation file is required. Supply gene_anno_file or set ",
      "options(", option_name, " = '/path/to/gene_annotation.txt')."
    )
  }

  if (xor(is.null(tr_filter), is.null(ctr_filter))) {
    stop("tr_filter and ctr_filter must be supplied together.")
  }
  if (is.null(tr_filter)) {
    if (is.null(tr) || is.null(ctr)) {
      stop("tr and ctr are required when custom filters are not supplied.")
    }
    tr_filter <- stats::setNames(list(tr), group_col)
    ctr_filter <- stats::setNames(list(ctr), group_col)

    if (!is.null(library)) {
      tr_filter[[library_col]] <- library
      ctr_filter[[library_col]] <- library
    }
    if (!is.null(batch)) {
      tr_filter[[batch_col]] <- batch
      ctr_filter[[batch_col]] <- batch
    }
  }

  filter_columns <- unique(c(names(tr_filter), names(ctr_filter)))
  if (length(filter_columns) == 0 || any(!nzchar(filter_columns))) {
    stop("tr_filter and ctr_filter must be named lists.")
  }
  missing_filter_columns <- setdiff(filter_columns, colnames(design.df))
  if (length(missing_filter_columns) > 0) {
    stop(
      "Columns used by the sample filters were not found in design.df: ",
      paste(missing_filter_columns, collapse = ", ")
    )
  }

  if (is.null(tr_name)) {
    tr_name <- if (!is.null(tr)) {
      tr
    } else {
      paste(unlist(tr_filter), collapse = ".")
    }
  }
  if (is.null(ctr_name)) {
    ctr_name <- if (!is.null(ctr)) {
      ctr
    } else {
      paste(unlist(ctr_filter), collapse = ".")
    }
  }

  match_design_filter <- function(design, filter_list) {
    keep <- rep(TRUE, nrow(design))
    for (filter_col in names(filter_list)) {
      keep <- keep & design[[filter_col]] %in% filter_list[[filter_col]]
    }
    keep
  }

  tr_idx <- match_design_filter(design.df, tr_filter)
  ctr_idx <- match_design_filter(design.df, ctr_filter)
  if (any(tr_idx & ctr_idx)) {
    stop("tr_filter and ctr_filter matched overlapping samples.")
  }
  if (!any(tr_idx)) {
    stop("tr_filter matched no samples.")
  }
  if (!any(ctr_idx)) {
    stop("ctr_filter matched no samples.")
  }
  design <- design.df[tr_idx | ctr_idx, , drop = FALSE]
  design$ContrastGroup <- ifelse(tr_idx[tr_idx | ctr_idx], tr_name, ctr_name)
  design$ContrastGroup <- factor(
    design$ContrastGroup,
    levels = c(ctr_name, tr_name)
  )
  design$SampleID <- rownames(design)

  prepare_expression_matrix <- function(expression_matrix, sample_ids, label) {
    if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
      stop(label, " must be a matrix or data.frame.")
    }
    missing_samples <- setdiff(sample_ids, colnames(expression_matrix))
    if (length(missing_samples) > 0) {
      stop(
        label,
        " is missing samples selected from design.df: ",
        paste(missing_samples, collapse = ", ")
      )
    }

    if (!is.null(gene_id_col) && gene_id_col %in% colnames(expression_matrix)) {
      gene_ids <- as.character(expression_matrix[[gene_id_col]])
    } else {
      gene_ids <- rownames(expression_matrix)
    }
    if (is.null(gene_ids) || any(is.na(gene_ids) | gene_ids == "")) {
      stop(
        label,
        " requires gene IDs in rownames or in the gene_id_col column."
      )
    }
    if (anyDuplicated(gene_ids)) {
      stop(label, " contains duplicated gene IDs.")
    }

    output <- as.matrix(expression_matrix[, sample_ids, drop = FALSE])
    suppressWarnings(storage.mode(output) <- "numeric")
    if (anyNA(output)) {
      stop(label, " contains missing or non-numeric values.")
    }
    rownames(output) <- gene_ids
    output
  }

  counts <- prepare_expression_matrix(
    count.matrix,
    rownames(design),
    "count.matrix"
  )
  tpm.mat <- NULL
  if (!is.null(tpm.matrix)) {
    tpm.mat <- prepare_expression_matrix(
      tpm.matrix,
      rownames(design),
      "tpm.matrix"
    )
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = design,
    design = ~ContrastGroup
  )
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  inter <- resultsNames(dds)
  resOrdered <- res[order(res$padj, na.last = TRUE), ]
  DESeq2::plotMA(resOrdered, alpha = 0.05)
  deseq2_result <- data.frame(resOrdered)

  gene_id_anno <- utils::read.table(
    gene_anno_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = 1
  )
  if (species == "human") {
    gene_id_anno <- gene_id_anno[
      !grepl("_PAR_", rownames(gene_id_anno)),
      ,
      drop = FALSE
    ]
  }
  rownames(gene_id_anno) <- sub("\\..*$", "", rownames(gene_id_anno))

  deseq2_result <- merge(
    gene_id_anno,
    deseq2_result,
    by.x = 0,
    by.y = 0,
    all = FALSE
  )

  if (!is.null(tpm.mat)) {
    deseq2_result <- merge(
      deseq2_result,
      tpm.mat,
      by.x = 1,
      by.y = 0,
      all = FALSE
    )
    ctr.samples <- rownames(design[design$ContrastGroup == ctr_name, ])
    tr.samples <- rownames(design[design$ContrastGroup == tr_name, ])
    deseq2_result$TPM_ctr_mean <- round(rowMeans(
      deseq2_result[, ctr.samples, drop = FALSE]
    ), digits = 3)
    deseq2_result$TPM_tr_mean <- round(rowMeans(
      deseq2_result[, tr.samples, drop = FALSE]
    ), digits = 3)
  }

  if ("Symbol" %in% colnames(deseq2_result)) {
    deseq2_result$Entrez <- MAGeCKFlute::TransGeneID(
      deseq2_result$Symbol,
      "Symbol",
      "Entrez",
      organism = organism
    )
  } else {
    warning(
      "The annotation table has no Symbol column; Entrez IDs were not added."
    )
    deseq2_result$Entrez <- NA_character_
  }
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

  result <- list()
  result[["inter"]] <- inter
  result[["result"]] <- deseq2_result
  result[["vsd"]] <- vsd
  result[["design"]] <- design
  result[["species"]] <- species
  return(result)
}


DESeq2_DEG_analysis <- function(
  count.matrix,
  tpm.matrix = NULL,
  design.df,
  tr,
  ctr
) {
  #SampleID Group

  stopifnot(all(c("Group") %in% colnames(design.df)))

  design <- design.df |> filter(Group %in% c(ctr, tr))
  design$Group <- factor(design$Group, levels = c(ctr, tr))
  design$SampleID <- rownames(design)
  counts <- count.matrix[, rownames(design)]
  rownames(counts) <- count.matrix$GID
  tpm.mat <- tpm.matrix[, rownames(design)]
  rownames(tpm.mat) <- tpm.matrix$GID

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = design,
    design = ~Group
  )
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  inter <- resultsNames(dds)
  resOrdered <- res[order(res$padj, na.last = TRUE), ]
  DESeq2::plotMA(resOrdered, alpha = 0.05)
  deseq2_result <- data.frame(resOrdered)
  gene_id_anno <- read.table(
    "~/Bin/gene_len.v43.new.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = 1
  )
  gene_id_anno <- gene_id_anno[!grepl("_PAR_", rownames(gene_id_anno)), ]

  rownames(gene_id_anno) <- sub("\\..*$", "", rownames(gene_id_anno))
  deseq2_result <- merge(
    gene_id_anno,
    deseq2_result,
    by.x = 0,
    by.y = 0,
    all = F
  )
  deseq2_result <- merge(deseq2_result, tpm.mat, by.x = 1, by.y = 0)

  if (is.data.frame(tpm.matrix)) {
    ctr.samples <- rownames(design[design$Group == ctr, ])
    tr.samples <- rownames(design[design$Group == tr, ])
    deseq2_result$TPM_ctr_mean <- round(
      apply(deseq2_result[, ctr.samples], 1, mean),
      digits = 3
    )
    deseq2_result$TPM_tr_mean <- round(
      apply(deseq2_result[, tr.samples], 1, mean),
      digits = 3
    )
  }
  deseq2_result$Entrez <- TransGeneID(
    deseq2_result$Symbol,
    "Symbol",
    "Entrez",
    organism = "hsa"
  )
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  result <- list()
  result[["inter"]] <- inter
  result[["result"]] <- deseq2_result
  result[["vsd"]] <- vsd
  result[["design"]] <- design
  return(result)
}


DESeq2_DEG_extract <- function(
  result.df,
  design.df,
  plot,
  adjust,
  fc = 1.5,
  pv = 0.05
) {
  deseq2_result <- result.df$result
  if (adjust) {
    result.df[["up"]] <- deseq2_result |>
      filter(!is.na(padj)) |>
      filter(log2FoldChange > log2(fc) & padj < pv)
    result.df[["dw"]] <- deseq2_result |>
      filter(!is.na(padj)) |>
      filter(log2FoldChange < -log2(fc) & padj < pv)
  } else {
    result.df[["up"]] <- deseq2_result |>
      filter(!is.na(pvalue)) |>
      filter(log2FoldChange > log2(fc) & pvalue < pv)
    result.df[["dw"]] <- deseq2_result |>
      filter(!is.na(pvalue)) |>
      filter(log2FoldChange < -log2(fc) & pvalue < pv)
  }
  if (plot) {
    #col_anno=data.frame(row.names = rownames(design.df),Group=design.df$Group,WT1=design.df$WT1_tpm)
    col_anno <- data.frame(
      row.names = design.df$SampleID,
      Group = design.df$Group
    )
    col_anno$Group <- as.factor(col_anno$Group)
    #col_anno_color <- list(Group = c(CD34="#F15A24",APL="#0071BC"))
    #col_anno_color <- list(Group = c(control="#F15A24",HMGA2_sh2="#0071BC"))
    col_anno_color <- list(Group = c("#A81E2C", "#08537C"))
    names(col_anno_color$Group) <- levels(col_anno$Group)

    EX_data <- assay(result.df$vsd[c(
      result.df$up$Row.names,
      result.df$dw$Row.names
    )])
    sd_rows <- apply(EX_data, 1, sd)
    EX_data <- EX_data[sd_rows > 0, ]
    EX_data <- t(scale(t(EX_data)))
    EX_data[EX_data > 2] <- 2
    EX_data[EX_data < -2] <- -2

    #EX_data <- assay(result.df$vsd[c(result.df$up$Row.names,result.df$dw$Row.names)])
    #EX_data <- log2(EX_data+1)
    #EX_data.mean <- matrix(rep(apply(EX_data,1,mean),ncol(EX_data)),ncol=ncol(EX_data))
    #EX_data.sd <- matrix(rep(apply(EX_data,1,sd),ncol(EX_data)),ncol=ncol(EX_data))
    #EX_data <- (EX_data - EX_data.mean) / EX_data.sd
    #EX_data[EX_data > 2] =  2
    #EX_data[EX_data < -2] =  -2
    result.df[["plot"]] <- pheatmap::pheatmap(
      EX_data,
      scale = "none",
      show_colnames = T,
      show_rownames = F,
      annotation_col = col_anno,
      annotation_colors = col_anno_color,
      cluster_cols = T,
      cluster_rows = T,
      clustering_method = "complete",
      color = colorRampPalette(rev(brewer.pal(n = 11, name = "PRGn")))(100)
    )
  }
  return(result.df)
}


Compare_pairwise_Deseq2 <- function(
  deg1.df,
  deg2.df,
  label1 = "DESeq 1",
  label2 = "DESeq 2",
  gene.list = c("ZMIZ1", "MEF2D"),
  fc = 1.5,
  top = 5
) {
  #deg1.df <- shZMIZ1.DEGs.ls$MOLM13_SH4_b2$result
  #deg2.df <- shZMIZ1.DEGs.ls$MOLM13_4E_b5$result
  #label1 <- "MOLM13_SH4"
  #label2 <- "MOLM13_4E"
  #fc <- 1.5
  #top <- 5
  #gene.list <- c("ZMIZ1","MEF2D")
  deg1.df <- na.omit(deg1.df[, c(1:4, 6, 9, 10, ncol(deg1.df) - 1)])
  deg2.df <- na.omit(deg2.df[, c(
    1,
    6,
    9,
    10,
    ncol(deg2.df) - 1,
    ncol(deg2.df)
  )])
  deg.merge.df <- merge(deg1.df, deg2.df, by.x = 1, by.y = 1)

  corr.result <- cor.test(
    ~ log2FoldChange.x + log2FoldChange.y,
    deg.merge.df,
    method = "pearson"
  )
  text <- paste(
    "k=",
    round(
      lm(log2FoldChange.x ~ log2FoldChange.y, deg.merge.df)$coefficients[2],
      4
    ),
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

  max <- max(max(deg.merge.df$log2FoldChange.x), deg.merge.df$log2FoldChange.y)

  deg.merge.df$group <- "no"
  deg.merge.df[
    deg.merge.df$log2FoldChange.x > log2(fc) &
      deg.merge.df$log2FoldChange.y > log2(fc),
  ]$group <- "up"
  deg.merge.df[
    deg.merge.df$log2FoldChange.x < -log2(fc) &
      deg.merge.df$log2FoldChange.y < -log2(fc),
  ]$group <- "down"

  gene.list <- c(
    (deg.merge.df |>
      filter(group == "up") |>
      top_n(n = top, wt = log2FoldChange.y * log2FoldChange.x))$Symbol,
    gene.list
  )
  gene.list <- c(
    (deg.merge.df |>
      filter(group == "down") |>
      top_n(n = top, wt = log2FoldChange.y * log2FoldChange.x))$Symbol,
    gene.list
  )
  index <- match(gene.list, deg.merge.df$Symbol)
  index <- na.omit(index)

  deg.merge.df$color <- deg.merge.df$group
  deg.merge.df$color[index] <- "black"
  mycolour <- c("grey", "#A81E2C", "#08537C", "black")
  names(mycolour) <- c("no", "up", "down", "black")

  deg.merge.df$label <- ""
  deg.merge.df$label[index] <- deg.merge.df$Symbol[index]
  deg.merge.df[deg.merge.df$group == "no", ]$label <- ""

  p <- ggplot(deg.merge.df, aes(y = log2FoldChange.y, x = log2FoldChange.x)) +
    geom_point_rast(
      shape = 16,
      alpha = 0.6,
      show.legend = FALSE,
      aes(color = group)
    ) +
    geom_point(
      color = ifelse(deg.merge.df$label == "", NA, "black"),
      shape = 1,
      show.legend = FALSE
    ) +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    ggtitle(text) +
    ylim(-max, max) +
    xlim(-max, max) +
    geom_hline(yintercept = c(-log2(fc), log2(fc)), linetype = "dotted") +
    geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted") +
    ggrepel::geom_text_repel(
      aes(label = label, color = group),
      show.legend = FALSE,
      fontface = "bold",
      size = 2.5,
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      max.overlaps = 10000
    ) +
    annotate(
      "text",
      label = sig,
      y = max(deg.merge.df$log2FoldChange.y, na.rm = TRUE),
      x = median(
        deg.merge.df[deg.merge.df$log2FoldChange.x != 0, ]$log2FoldChange.x,
        na.rm = TRUE
      ),
      size = 7
    ) +
    xlab(paste0("log2FoldChange of ", label1)) +
    ylab(paste0("log2FoldChange of ", label2)) +
    scale_color_manual(values = mycolour) +
    scale_fill_manual(values = mycolour) +
    theme_test()
  return(p)
}
