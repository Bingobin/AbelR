# Functions for deseq2 analyses.

#' Build count, TPM, and design tables across bulk RNA-seq batches
#'
#' Reads one count table and one TPM table for each batch, appends the batch ID
#' to every sample name, joins the tables by gene ID, and constructs sample
#' metadata. Input sample columns are expected to follow
#' `<CellLine>_<TargetGene>_<Group>_<replicate>`, for example
#' `M13_KAP1_sh9_1`. `Group` is taken from the segment immediately before the
#' numeric replicate, producing `SH9` in this example.
#'
#' @param batches Character or numeric vector of batch identifiers.
#' @param target Target string used to locate input files.
#' @param library_name Value assigned to the `Library` column of the design
#'   table.
#' @param count_dir Directory containing count and TPM tables.
#' @param cellline_map Named character vector mapping cell-line prefixes to
#'   standardized names. Unmapped prefixes are retained unchanged.
#' @param id_col Gene-identifier column shared by count and TPM tables.
#'
#' @return A list containing the merged `count` and `tpm` data frames and a
#'   `design` data frame with `CellLine`, `Library`, `Batch`, and `Group`.
#' @export
build_bulkRNA_batches <- function(
  batches,
  target,
  library_name,
  count_dir = "counts",
  cellline_map = c(
    "M13" = "M13",
    "MOLM13" = "M13",
    "OA3" = "OA3"
  ),
  id_col = "GID"
) {
  if (!length(batches)) {
    stop("batches must contain at least one batch identifier.")
  }
  if (anyDuplicated(as.character(batches))) {
    stop("batches must not contain duplicate identifiers.")
  }
  if (length(target) != 1 || is.na(target) || !nzchar(target)) {
    stop("target must be one non-empty string.")
  }
  if (length(library_name) != 1 || is.na(library_name)) {
    stop("library_name must be a single value.")
  }
  if (!dir.exists(count_dir)) {
    stop("count_dir does not exist: ", count_dir)
  }
  if (is.null(names(cellline_map)) || any(!nzchar(names(cellline_map)))) {
    stop("cellline_map must be a named character vector.")
  }

  find_batch_file <- function(batch, type) {
    pattern <- file.path(
      count_dir,
      paste0(target, "*mR*", batch, ".merge.", type, ".txt")
    )
    matches <- Sys.glob(pattern)
    if (!length(matches)) {
      stop(type, " file not found for batch ", batch, ": ", pattern)
    }
    if (length(matches) > 1) {
      stop(
        "Multiple ",
        type,
        " files matched batch ",
        batch,
        ": ",
        paste(matches, collapse = ", ")
      )
    }
    matches
  }

  count_list <- list()
  tpm_list <- list()
  design_list <- list()

  for (batch in batches) {
    batch <- as.character(batch)
    count_file <- find_batch_file(batch, "count")
    tpm_file <- find_batch_file(batch, "tpm")

    count_df <- utils::read.delim(count_file, check.names = FALSE)
    tpm_df <- utils::read.delim(tpm_file, check.names = FALSE)

    if (!id_col %in% colnames(count_df)) {
      stop(id_col, " not found in: ", count_file)
    }
    if (!id_col %in% colnames(tpm_df)) {
      stop(id_col, " not found in: ", tpm_file)
    }
    if (anyDuplicated(count_df[[id_col]])) {
      stop("Duplicated ", id_col, " found in: ", count_file)
    }
    if (anyDuplicated(tpm_df[[id_col]])) {
      stop("Duplicated ", id_col, " found in: ", tpm_file)
    }

    sample_cols <- setdiff(colnames(count_df), id_col)
    if (!length(sample_cols)) {
      stop("No sample columns found in: ", count_file)
    }
    if (anyDuplicated(sample_cols)) {
      stop("Duplicated sample columns found in: ", count_file)
    }
    missing_in_tpm <- setdiff(sample_cols, colnames(tpm_df))
    if (length(missing_in_tpm)) {
      stop(
        "The following count samples are missing from TPM in batch ",
        batch,
        ": ",
        paste(missing_in_tpm, collapse = ", ")
      )
    }

    count_df <- count_df[, c(id_col, sample_cols), drop = FALSE]
    tpm_df <- tpm_df[, c(id_col, sample_cols), drop = FALSE]
    new_sample_cols <- paste0(sample_cols, ".", batch)
    colnames(count_df)[-1] <- new_sample_cols
    colnames(tpm_df)[-1] <- new_sample_cols

    raw_cellline <- sub("_.*$", "", sample_cols)
    cellline <- unname(cellline_map[raw_cellline])
    cellline[is.na(cellline)] <- raw_cellline[is.na(cellline)]

    valid_sample_name <- grepl("^[^_]+_.+_[^_]+_[0-9]+$", sample_cols)
    if (any(!valid_sample_name)) {
      stop(
        paste0(
          "Sample names must follow ",
          "<CellLine>_<TargetGene>_<Group>_<replicate>: "
        ),
        paste(sample_cols[!valid_sample_name], collapse = ", ")
      )
    }
    group <- sub("^.*_([^_]+)_[0-9]+$", "\\1", sample_cols)

    design_list[[batch]] <- data.frame(
      CellLine = cellline,
      Library = rep(as.character(library_name), length(sample_cols)),
      Batch = rep(batch, length(sample_cols)),
      Group = toupper(group),
      row.names = new_sample_cols,
      stringsAsFactors = FALSE
    )
    count_list[[batch]] <- count_df
    tpm_list[[batch]] <- tpm_df
  }

  count <- Reduce(
    function(x, y) dplyr::full_join(x, y, by = id_col),
    count_list
  )
  tpm <- Reduce(
    function(x, y) dplyr::full_join(x, y, by = id_col),
    tpm_list
  )
  design <- do.call(rbind, unname(design_list))

  list(count = count, tpm = tpm, design = design)
}


#' Plot PCA for merged bulk RNA-seq batches
#'
#' Uses the count and design tables returned by [build_bulkRNA_batches()],
#' applies the DESeq2 variance-stabilizing transformation, calculates sample
#' principal components, and draws a batch-aware PCA plot.
#'
#' @param bulk_data List returned by [build_bulkRNA_batches()] containing
#'   `count` and `design`.
#' @param id_col Gene-identifier column in `bulk_data$count`.
#' @param title Plot title.
#' @param color_by,shape_by,alpha_by Design-table columns mapped to point color,
#'   shape, and transparency.
#' @param alpha_values Optional named numeric vector defining transparency for
#'   each `alpha_by` level. When `NULL`, values from 0.6 to 1 are generated
#'   automatically.
#' @param point_size Point size.
#' @param label_samples Logical; label points with sample IDs.
#' @param label_size Sample-label size.
#' @param blind Logical passed to
#'   [DESeq2::varianceStabilizingTransformation()].
#' @param top_var_n Number of genes with the highest VST variance used for both
#'   PCA and sample-distance clustering. Use `NULL` to retain all genes.
#' @param cluster_heatmap Logical; calculate a sample-distance heatmap using the
#'   same genes as PCA.
#' @param annotation_columns Design-table columns displayed on the heatmap.
#' @param heatmap_title Sample-distance heatmap title.
#' @param interactive Logical; also create a Plotly version of the PCA plot.
#'   This requires the optional `plotly` package.
#'
#' @return A list containing the ggplot (`plot`), optional Plotly widget
#'   (`interactive_plot`), sample-distance heatmap (`heatmap`), distance matrix
#'   (`sample_distance_matrix`), `dist` object (`sample_distance`), selected
#'   genes (`selected_genes`), sample coordinates (`data`), explained variance
#'   (`variance`), the `prcomp` result (`pca`), and the transformed DESeq2 object
#'   (`vsd`).
#' @export
plot_bulkRNA_PCA <- function(
  bulk_data,
  id_col = "GID",
  title = "PCA of bulk RNA-seq samples",
  color_by = "Group",
  shape_by = "CellLine",
  alpha_by = "Batch",
  alpha_values = NULL,
  point_size = 3.5,
  label_samples = TRUE,
  label_size = 2.5,
  blind = TRUE,
  top_var_n = 5000,
  cluster_heatmap = TRUE,
  annotation_columns = c("Group", "Library", "Batch"),
  heatmap_title = "Sample-to-sample distance clustering",
  interactive = FALSE
) {
  if (!is.list(bulk_data) ||
      !all(c("count", "design") %in% names(bulk_data))) {
    stop("bulk_data must contain count and design components.")
  }
  count <- bulk_data$count
  design <- bulk_data$design
  if (!is.data.frame(count) || !is.data.frame(design)) {
    stop("bulk_data$count and bulk_data$design must be data frames.")
  }
  if (!id_col %in% colnames(count)) {
    stop(id_col, " not found in bulk_data$count.")
  }
  if (anyDuplicated(count[[id_col]])) {
    stop("Duplicated ", id_col, " values found in bulk_data$count.")
  }
  if (is.null(rownames(design)) || anyDuplicated(rownames(design))) {
    stop("bulk_data$design must have unique sample IDs as row names.")
  }

  aesthetics <- c(color_by, shape_by, alpha_by)
  missing_aesthetics <- setdiff(aesthetics, colnames(design))
  if (length(missing_aesthetics)) {
    stop(
      "Design columns not found: ",
      paste(missing_aesthetics, collapse = ", ")
    )
  }
  reserved_columns <- intersect(c("SampleID", "PC1", "PC2"), colnames(design))
  if (length(reserved_columns)) {
    stop(
      "Reserved PCA columns already exist in design: ",
      paste(reserved_columns, collapse = ", ")
    )
  }
  if (!is.null(top_var_n)) {
    if (length(top_var_n) != 1 ||
        !is.numeric(top_var_n) ||
        is.na(top_var_n) ||
        !is.finite(top_var_n) ||
        top_var_n < 1 ||
        top_var_n != floor(top_var_n) ||
        top_var_n > .Machine$integer.max) {
      stop("top_var_n must be NULL or one positive integer.")
    }
    top_var_n <- as.integer(top_var_n)
  }
  if (isTRUE(cluster_heatmap)) {
    missing_annotations <- setdiff(annotation_columns, colnames(design))
    if (length(missing_annotations)) {
      stop(
        "Heatmap annotation columns not found: ",
        paste(missing_annotations, collapse = ", ")
      )
    }
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      stop("Package 'pheatmap' is required when cluster_heatmap = TRUE.")
    }
  }

  pca_samples <- rownames(design)
  missing_samples <- setdiff(pca_samples, colnames(count))
  if (length(missing_samples)) {
    stop(
      "Design samples missing from the count table: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  if (length(pca_samples) < 2) {
    stop("At least two samples are required for PCA.")
  }
  pca_count_df <- count[, pca_samples, drop = FALSE]
  numeric_columns <- vapply(
    pca_count_df,
    function(x) is.numeric(x) || is.integer(x),
    logical(1)
  )
  if (any(!numeric_columns)) {
    stop(
      "Non-numeric count columns: ",
      paste(colnames(pca_count_df)[!numeric_columns], collapse = ", ")
    )
  }

  pca_count <- round(as.matrix(pca_count_df))
  storage.mode(pca_count) <- "integer"
  rownames(pca_count) <- as.character(count[[id_col]])
  if (anyNA(pca_count) || any(pca_count < 0)) {
    stop("Counts must be non-missing, non-negative integer-like values.")
  }
  pca_count <- pca_count[rowSums(pca_count) > 0, , drop = FALSE]
  if (nrow(pca_count) < 2) {
    stop("At least two genes with non-zero counts are required for PCA.")
  }

  pca_dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pca_count,
    colData = design,
    design = ~1
  )
  pca_vsd <- DESeq2::varianceStabilizingTransformation(
    pca_dds,
    blind = blind
  )
  vsd_mat <- SummarizedExperiment::assay(pca_vsd)
  gene_var <- matrixStats::rowVars(vsd_mat, useNames = TRUE)
  variable_order <- order(gene_var, decreasing = TRUE, na.last = NA)
  if (!is.null(top_var_n)) {
    variable_order <- utils::head(variable_order, top_var_n)
  }
  if (length(variable_order) < 2) {
    stop("Fewer than two genes with finite VST variance are available.")
  }
  selected_genes <- rownames(vsd_mat)[variable_order]
  analysis_mat <- vsd_mat[variable_order, , drop = FALSE]

  pca_res <- stats::prcomp(
    t(analysis_mat),
    scale. = FALSE
  )
  if (ncol(pca_res$x) < 2) {
    stop("PCA produced fewer than two principal components.")
  }
  pca_var <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))

  sample_pca_df <- data.frame(
    SampleID = rownames(pca_res$x),
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    design[rownames(pca_res$x), , drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  sample_pca_df[[color_by]] <- factor(sample_pca_df[[color_by]])
  sample_pca_df[[shape_by]] <- factor(sample_pca_df[[shape_by]])
  sample_pca_df[[alpha_by]] <- factor(sample_pca_df[[alpha_by]])

  alpha_levels <- levels(sample_pca_df[[alpha_by]])
  if (is.null(alpha_values)) {
    alpha_values <- if (length(alpha_levels) == 1) {
      1
    } else {
      seq(0.6, 1, length.out = length(alpha_levels))
    }
    names(alpha_values) <- alpha_levels
  } else {
    if (!is.numeric(alpha_values) || any(alpha_values < 0 | alpha_values > 1)) {
      stop("alpha_values must contain numeric values between 0 and 1.")
    }
    if (is.null(names(alpha_values))) {
      if (length(alpha_values) != length(alpha_levels)) {
        stop("Unnamed alpha_values must match the number of alpha levels.")
      }
      names(alpha_values) <- alpha_levels
    }
    missing_alpha <- setdiff(alpha_levels, names(alpha_values))
    if (length(missing_alpha)) {
      stop(
        "alpha_values missing levels: ",
        paste(missing_alpha, collapse = ", ")
      )
    }
  }

  p <- ggplot2::ggplot(
    sample_pca_df,
    ggplot2::aes(
      x = .data[["PC1"]],
      y = .data[["PC2"]],
      color = .data[[color_by]],
      shape = .data[[shape_by]],
      text = .data[["SampleID"]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(alpha = .data[[alpha_by]]),
      size = point_size
    ) +
    ggplot2::scale_alpha_manual(values = alpha_values, name = alpha_by) +
    ggplot2::theme_test() +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
      color = color_by,
      shape = shape_by
    )

  if (isTRUE(label_samples)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = .data[["SampleID"]]),
        size = label_size,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[["SampleID"]]),
        size = label_size,
        vjust = -0.8,
        show.legend = FALSE
      )
    }
  }

  interactive_plot <- NULL
  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required when interactive = TRUE.")
    }
    interactive_plot <- plotly::ggplotly(p, tooltip = "text")
  }

  sample_distance <- NULL
  sample_distance_matrix <- NULL
  heatmap <- NULL
  if (isTRUE(cluster_heatmap)) {
    sample_distance <- stats::dist(t(analysis_mat), method = "euclidean")
    sample_distance_matrix <- as.matrix(sample_distance)
    annotation <- design[colnames(analysis_mat), annotation_columns, drop = FALSE]
    heatmap <- pheatmap::pheatmap(
      sample_distance_matrix,
      clustering_distance_rows = sample_distance,
      clustering_distance_cols = sample_distance,
      annotation_col = annotation,
      annotation_row = annotation,
      color = grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(9, "Blues"))
      )(255),
      main = heatmap_title,
      fontsize = 9,
      fontsize_row = 8,
      fontsize_col = 8,
      border_color = NA,
      silent = TRUE
    )
  }

  list(
    plot = p,
    interactive_plot = interactive_plot,
    heatmap = heatmap,
    sample_distance = sample_distance,
    sample_distance_matrix = sample_distance_matrix,
    selected_genes = selected_genes,
    gene_variance = gene_var,
    top_var_n = top_var_n,
    data = sample_pca_df,
    variance = pca_var,
    pca = pca_res,
    vsd = pca_vsd
  )
}


#' Plot median TPM expression for one gene
#'
#' Uses the TPM and design tables returned by [build_bulkRNA_batches()] to
#' display one gene across cell lines, batches, groups, and libraries. Bars
#' show median `log2(TPM + 1)`, error bars show the interquartile range, and
#' points show individual samples.
#'
#' @param bulk_data List returned by [build_bulkRNA_batches()] containing
#'   `tpm` and `design` data frames.
#' @param gene One gene symbol or Ensembl gene ID. Ensembl version suffixes are
#'   ignored during matching.
#' @param id_col Gene-identifier column in `bulk_data$tpm`.
#' @param group_colors Character vector of valid R colors. Named colors are
#'   matched to `Group` levels when all levels are present; otherwise colors
#'   are recycled in group-level order.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
plot_gene_tpm_median <- function(
  bulk_data,
  gene,
  id_col = "GID",
  group_colors = get("mycol", envir = asNamespace("AbelR"))
) {
  if (
    !is.list(bulk_data) ||
      !all(c("tpm", "design") %in% names(bulk_data))
  ) {
    stop("bulk_data must contain both tpm and design.")
  }
  if (!is.character(gene) ||
      length(gene) != 1L ||
      is.na(gene) ||
      !nzchar(gene)) {
    stop("gene must be one non-empty gene symbol or Ensembl ID.")
  }

  tpm <- bulk_data$tpm
  design <- bulk_data$design
  if (!is.data.frame(tpm) || !is.data.frame(design)) {
    stop("bulk_data$tpm and bulk_data$design must be data frames.")
  }
  required_design_cols <- c("CellLine", "Batch", "Group", "Library")
  missing_design_cols <- setdiff(required_design_cols, colnames(design))
  if (length(missing_design_cols)) {
    stop(
      "Missing design columns: ",
      paste(missing_design_cols, collapse = ", ")
    )
  }
  if (!id_col %in% colnames(tpm)) {
    stop("ID column not found in tpm: ", id_col)
  }
  if (anyDuplicated(tpm[[id_col]])) {
    stop("Duplicated ", id_col, " values found in bulk_data$tpm.")
  }
  if (is.null(rownames(design)) ||
      anyNA(rownames(design)) ||
      any(!nzchar(rownames(design))) ||
      anyDuplicated(rownames(design))) {
    stop("bulk_data$design must have unique sample IDs as row names.")
  }

  sample_ids <- rownames(design)
  missing_samples <- setdiff(sample_ids, colnames(tpm))
  if (length(missing_samples)) {
    stop(
      "Design samples missing from the TPM table: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  numeric_samples <- vapply(
    tpm[, sample_ids, drop = FALSE],
    function(x) is.numeric(x) || is.integer(x),
    logical(1)
  )
  if (any(!numeric_samples)) {
    stop(
      "Non-numeric TPM sample columns: ",
      paste(sample_ids[!numeric_samples], collapse = ", ")
    )
  }

  tpm_gene_ids <- sub("\\..*$", "", as.character(tpm[[id_col]]))
  query_id <- sub("\\..*$", "", gene)
  gene_rows <- which(tpm_gene_ids == query_id)

  if (!length(gene_rows) && "Symbol" %in% colnames(tpm)) {
    gene_rows <- which(
      toupper(as.character(tpm$Symbol)) == toupper(gene)
    )
  }

  if (!length(gene_rows)) {
    if (
      !requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)
    ) {
      stop(
        "Gene symbol mapping requires AnnotationDbi and org.Hs.eg.db, ",
        "or supply an Ensembl ID directly."
      )
    }
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = toupper(gene),
      keytype = "SYMBOL",
      column = "ENSEMBL",
      multiVals = "list"
    )
    mapped_ids <- unique(stats::na.omit(as.character(unlist(mapped))))
    gene_rows <- which(tpm_gene_ids %in% mapped_ids)
  }

  if (!length(gene_rows)) {
    stop("Gene not found in bulk_data$tpm: ", gene)
  }
  if (length(gene_rows) > 1L) {
    stop(
      "Gene matched more than one TPM row: ",
      paste(tpm[[id_col]][gene_rows], collapse = ", ")
    )
  }

  plot_df <- data.frame(
    SampleID = sample_ids,
    design[sample_ids, , drop = FALSE],
    TPM = as.numeric(tpm[gene_rows, sample_ids, drop = TRUE]),
    check.names = FALSE,
    row.names = NULL
  )
  plot_df <- plot_df[!is.na(plot_df$TPM), , drop = FALSE]
  if (!nrow(plot_df)) {
    stop("All TPM values are missing for gene: ", gene)
  }
  if (any(!is.finite(plot_df$TPM)) || any(plot_df$TPM < 0)) {
    stop("TPM values must be finite and non-negative.")
  }

  plot_df$log2_TPM_plus_1 <- log2(plot_df$TPM + 1)
  plot_df$Batch <- factor(as.character(plot_df$Batch))
  group_levels <- unique(as.character(plot_df$Group))
  if ("SCR" %in% group_levels) {
    group_levels <- c("SCR", setdiff(group_levels, "SCR"))
  }
  plot_df$Group <- factor(
    as.character(plot_df$Group),
    levels = group_levels
  )

  if (!is.character(group_colors) || !length(group_colors) || anyNA(group_colors)) {
    stop("group_colors must contain at least one valid R color.")
  }
  valid_colors <- vapply(group_colors, function(color) {
    tryCatch(
      {
        grDevices::col2rgb(color)
        TRUE
      },
      error = function(e) FALSE
    )
  }, logical(1))
  if (any(!valid_colors)) {
    stop(
      "Invalid colors in group_colors: ",
      paste(group_colors[!valid_colors], collapse = ", ")
    )
  }
  if (!is.null(names(group_colors)) &&
      all(group_levels %in% names(group_colors))) {
    group_palette <- group_colors[group_levels]
  } else {
    group_palette <- stats::setNames(
      rep(group_colors, length.out = length(group_levels)),
      group_levels
    )
  }

  summary_df <- plot_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(c("CellLine", "Batch", "Library", "Group"))
      )
    ) |>
    dplyr::summarise(
      n = dplyr::n(),
      Median = stats::median(.data[["log2_TPM_plus_1"]], na.rm = TRUE),
      Q1 = stats::quantile(
        .data[["log2_TPM_plus_1"]],
        0.25,
        na.rm = TRUE
      ),
      Q3 = stats::quantile(
        .data[["log2_TPM_plus_1"]],
        0.75,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  dodge <- ggplot2::position_dodge(width = 0.8)
  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = .data[["Group"]],
      y = .data[["Median"]],
      fill = .data[["Group"]],
      group = .data[["Batch"]]
    )
  ) +
    ggplot2::geom_col(
      position = dodge,
      width = 0.7,
      color = "black",
      linewidth = 0.35,
      alpha = 0.9
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data[["Q1"]],
        ymax = .data[["Q3"]]
      ),
      position = dodge,
      width = 0.3,
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        x = .data[["Group"]],
        y = .data[["log2_TPM_plus_1"]],
        group = .data[["Batch"]]
      ),
      position = ggplot2::position_jitter(
        width = 0.08,
        height = 0
      ),
      inherit.aes = FALSE,
      shape = 16,
      size = 3,
      color = "black",
      alpha = 0.6
    ) +
    ggplot2::scale_fill_manual(values = group_palette, drop = FALSE) +
    ggplot2::facet_grid(
      ~ CellLine + Batch + Library,
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::labs(
      title = paste0(gene, " expression"),
      subtitle = "Bars: median; points: samples; error bars: IQR",
      x = "Group",
      y = expression(log[2](TPM + 1)),
      fill = "Group"
    ) +
    ggplot2::theme_test(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Run a filtered two-group DESeq2 analysis
#'
#' Selects treatment and control samples from a shared design table, fits a
#' DESeq2 model, joins human or mouse gene annotation and optional TPM values,
#' and calculates a variance-stabilized expression object.
#'
#' @param count.matrix Count data frame with a `GID` column and sample columns.
#' @param tpm.matrix Optional TPM data frame using the same `GID` and sample
#'   columns as `count.matrix`.
#' @param design.df Sample metadata whose row names match matrix sample columns.
#'   It must contain `Group`, `Batch`, and `Library`.
#' @param library,batch,tr,ctr Values selecting one library, batch, treatment
#'   group, and control group when custom filters are not supplied.
#' @param tr_filter,ctr_filter Optional named lists defining treatment and
#'   control samples from arbitrary columns in `design.df`.
#' @param tr_name,ctr_name Optional display names for the resulting contrast
#'   groups.
#' @param species Species used for annotation and ID conversion; `"human"` or
#'   `"mouse"`.
#' @param gene_anno_file Optional annotation table (`.txt` or `.txt.gz`). When
#'   `NULL`, AbelR uses its bundled species-specific annotation.
#'
#' @return A list containing DESeq2 coefficient names (`inter`), the annotated
#'   result table (`result`), the variance-stabilized object (`vsd`), selected
#'   sample metadata (`design`), and normalized `species`.
#' @export
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
  gene_anno_file = NULL
) {
  species <- .abel_normalize_species(species)
  organism <- if (species == "human") "hsa" else "mmu"
  stopifnot(all(c("Group", "Batch", "Library") %in% colnames(design.df)))
  stopifnot("GID" %in% colnames(count.matrix))
  if (!is.null(tpm.matrix)) {
    stopifnot("GID" %in% colnames(tpm.matrix))
  }

  if (xor(is.null(tr_filter), is.null(ctr_filter))) {
    stop("tr_filter and ctr_filter must be supplied together.")
  }
  if (is.null(tr_filter)) {
    if (is.null(tr) || is.null(ctr)) {
      stop("tr and ctr are required when custom filters are not supplied.")
    }
    if (is.null(library) || is.null(batch)) {
      stop("library and batch are required when custom filters are not supplied.")
    }
    tr_filter <- list(Library = library, Batch = batch, Group = tr)
    ctr_filter <- list(Library = library, Batch = batch, Group = ctr)
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

  counts <- as.matrix(count.matrix[, rownames(design), drop = FALSE])
  rownames(counts) <- count.matrix$GID
  tpm.mat <- NULL
  if (!is.null(tpm.matrix)) {
    tpm.mat <- as.matrix(tpm.matrix[, rownames(design), drop = FALSE])
    rownames(tpm.mat) <- tpm.matrix$GID
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

  gene_id_anno <- .abel_gene_annotation(species, gene_anno_file)

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
    if (!requireNamespace("MAGeCKFlute", quietly = TRUE)) {
      stop("Package 'MAGeCKFlute' is required to add Entrez IDs.")
    }
    deseq2_result$Entrez <- suppressPackageStartupMessages(
      MAGeCKFlute::TransGeneID(
        deseq2_result$Symbol,
        "Symbol",
        "Entrez",
        organism = organism
      )
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


#' Run a two-group DESeq2 analysis
#'
#' Convenience interface to [DESeq2_DEG_analysis_batch()] that selects samples
#' using only the `Group` column while retaining the standard AbelR matrix and
#' metadata schema.
#'
#' @param count.matrix Count data frame with a `GID` column and sample columns.
#' @param tpm.matrix Optional TPM data frame with the same identifiers and sample
#'   columns.
#' @param design.df Sample metadata with row names matching sample columns and
#'   required `Group`, `Batch`, and `Library` columns.
#' @param tr,ctr Treatment and control values from `design.df$Group`.
#' @param species Species used for gene annotation; `"human"` or `"mouse"`.
#' @param gene_anno_file Optional custom annotation table.
#'
#' @return The result list returned by [DESeq2_DEG_analysis_batch()].
#' @export
DESeq2_DEG_analysis <- function(
  count.matrix,
  tpm.matrix = NULL,
  design.df,
  tr,
  ctr,
  species = c("human", "mouse"),
  gene_anno_file = NULL
) {
  DESeq2_DEG_analysis_batch(
    count.matrix = count.matrix,
    tpm.matrix = tpm.matrix,
    design.df = design.df,
    tr = tr,
    ctr = ctr,
    tr_filter = list(Group = tr),
    ctr_filter = list(Group = ctr),
    species = species,
    gene_anno_file = gene_anno_file
  )
}


#' Extract significant genes from an AbelR DESeq2 result
#'
#' Divides a DESeq2 result into upregulated and downregulated tables using fold
#' change and significance thresholds and optionally constructs a scaled
#' variance-stabilized heatmap.
#'
#' @param result.df Result list returned by [DESeq2_DEG_analysis()] or
#'   [DESeq2_DEG_analysis_batch()].
#' @param design.df Sample metadata containing `SampleID` and `Group`, used for
#'   heatmap annotation.
#' @param plot Logical; create and store a heatmap.
#' @param adjust Logical; use adjusted P values when `TRUE`, otherwise raw
#'   P values.
#' @param fc Fold-change threshold on the linear scale.
#' @param pv Significance threshold.
#'
#' @return The input result list with `up` and `dw` tables and, when requested,
#'   a `plot` component containing the pheatmap result.
#' @export
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

    EX_data <- SummarizedExperiment::assay(result.df$vsd[c(
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


#' Compare fold changes from two DESeq2 analyses
#'
#' Merges two annotated DESeq2 result tables, calculates their Pearson
#' correlation and regression slope, identifies concordant upregulated and
#' downregulated genes, and labels selected or top-ranked genes.
#'
#' @param deg1.df,deg2.df Annotated DESeq2 result data frames generated with a
#'   compatible AbelR result layout.
#' @param label1,label2 Axis labels identifying the two comparisons.
#' @param gene.list Character vector of gene symbols to label.
#' @param fc Fold-change threshold on the linear scale.
#' @param top Number of top concordant genes to add from each direction.
#'
#' @return A [ggplot2::ggplot] comparison plot.
#' @export
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
