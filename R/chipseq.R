# Functions for chipseq analyses.

#' Plot ranked enhancer signals
#'
#' Draws enhancer signal against rank for one sample, marks the super-enhancer
#' boundary, and labels selected genes and the highest-ranked enhancers.
#'
#' @param label Character vector of nearby gene symbols to label.
#' @param sample Sample name selected from `rank.df$Sample`.
#' @param top Number of top-ranked enhancers to label.
#' @param rank.df Data frame containing `Sample`, `Signal`, `Type`, `Rank`, and
#'   `CLOSEST_GENE` columns.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
enhancer_rank_plot <- function(label, sample, top, rank.df) {
  ###############################################
  #  Sample    Signal Type Rank CLOSEST_GENE
  #  CHH 555644.95   SE    1         ETV6
  #  CHH 429203.43   SE    2       FNDC3B
  #  CHH 235318.70   SE    3        RREB1
  #  CHH 221833.55   SE    4      COL23A1
  #  CHH 198674.49   SE    5   BZRAP1-AS1
  #  CHH 185526.61   SE    6        CELF2
  ##############################################
  #  sig_gene <- "WT1"
  #  sample <- "CHH"
  #  top <- 5
  enhancer_signal <- rank.df %>% filter(Sample == sample)
  p <- ggplot(enhancer_signal, aes(x = Rank / 1000, y = Signal)) +
    geom_point_rast(aes(color = Type), alpha = 1, size = 0.5) +
    xlab("Rank/1000") +
    ylab("Enhancer signal") +
    ggtitle(paste0(sample, " H3K27ac")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    blank +
    scale_color_manual(
      values = c("#A81E2C", brewer.pal(n = 9, name = "Set1")[9])
    ) +
    geom_vline(
      xintercept = nrow(enhancer_signal[enhancer_signal$Type == "SE", ]) / 1000,
      lty = 4,
      col = "black",
      lwd = 0.8
    )

  #p <- p + ggrepel::geom_label_repel(data = enhancer_signal[enhancer_signal$CLOSEST_GENE %in% label | enhancer_signal$Rank %in% 1:top,],
  #aes(x = Rank/1000, y = Signal, label = CLOSEST_GENE, fill = Type),
  p <- p +
    ggrepel::geom_text_repel(
      data = enhancer_signal[
        enhancer_signal$CLOSEST_GENE %in%
          label |
          enhancer_signal$Rank %in% 1:top,
      ],
      aes(x = Rank / 1000, y = Signal, label = CLOSEST_GENE, color = Type),
      fontface = "bold",
      size = 2.5,
      box.padding = unit(1, "lines"),
      segment.color = brewer.pal(n = 9, name = "Set1")[2],
      point.padding = unit(0.3, "lines"),
      segment.size = 0.3,
      show.legend = FALSE,
      label.r = 0.5,
      nudge_x = 0,
      nudge_y = 0,
      max.overlaps = 10000
    )
  #+
  #scale_fill_manual(values=c(brewer.pal(n = 9, name ="Set3")[4], brewer.pal(n = 9, name ="Set3")[9]))
  p <- p +
    theme(
      legend.title = element_blank(),
      legend.position = c(.95, .95),
      legend.justification = c("right", "top")
    )
  p <- p +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    )
  se_num <- dim(enhancer_signal |> filter(Type == "SE"))[1]
  cutoff <- enhancer_signal$Signal[se_num + 1]
  p <- p +
    annotate(
      "text",
      x = max(enhancer_signal$Rank) / 1000 / 2,
      y = max(enhancer_signal$Signal) / 2,
      label = paste0(
        "Cutoff used: ",
        cutoff,
        "\nSuper-Enhancers identified: ",
        se_num
      ),
      hjust = 0
    )
  return(p)
}


#' Compare gene expression across risk groups
#'
#' Combines disease risk-group expression with normal-sample expression and
#' draws log2 TPM boxplots with pairwise statistical comparisons.
#'
#' @param gene Gene symbol to plot.
#' @param risk.df Risk-model data frame containing the gene-expression column,
#'   risk group, and mutation columns used by the function.
#' @param nm.df Normal-sample expression table containing `Symbol` and sample
#'   expression columns.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
expression_boxplot_byMut <- function(gene, risk.df, nm.df) {
  #  gene <- "CHRNE"
  #  risk.df <- rs.df.clinc.beat
  #  nm.df <- aml_all_tpm.matrix.nm
  index1 <- match(gene, nm.df$Symbol)
  index2 <- match(gene, colnames(risk.df))
  index3 <- match(c("WT1.x", "TET2.x", "IDH1.x", "IDH2.x"), colnames(risk.df))
  group <- rowSums(risk.df[, index3])
  group[group > 0] <- "MUT"
  group[group == 0] <- "WT"
  gg1 <- data.frame(
    Group = rep("NM", length(index1)),
    TPM = as.numeric(nm.df[index1, 4:ncol(nm.df)][1, ]),
    Risk = rep("NM", length(index1))
  )
  gg2 <- data.frame(
    Group = group,
    TPM = as.numeric(risk.df[, index2]),
    Risk = risk.df[, 4]
  )

  gg <- rbind(gg1, gg2)
  gg$Group <- factor(gg$Group, levels = c("MUT", "WT", "NM"))
  gg$Risk <- factor(gg$Risk, levels = c("High_risk", "Low_risk", "NM"))
  ggplot(gg, aes(x = Risk, y = log2(TPM), color = Risk, fill = Risk)) +
    geom_boxplot(outlier.shape = 21, show.legend = FALSE, alpha = 0.65) +
    #  geom_jitter(size=0.5,show.legend = FALSE) +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1),
      panel.background = element_blank()
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(values = c("#FF7F00", "#1F78B4", "grey50")) +
    scale_fill_manual(values = c("#FF7F00", "#1F78B4", "grey50")) +
    #    stat_compare_means(show.legend = FALSE, comparisons = list(c("MUT","WT"),c("WT","NM"),c("MUT","NM"))) +
    stat_compare_means(
      show.legend = FALSE,
      comparisons = list(
        c("High_risk", "Low_risk"),
        c("Low_risk", "NM"),
        c("High_risk", "NM")
      )
    ) +
    ggtitle(paste0(gene, " Expression")) +
    xlab("")
}


#' Plot HOMER peak annotation proportions
#'
#' Extracts broad genomic feature categories from a HOMER `Annotation` column,
#' calculates category frequencies, and draws a labelled pie chart.
#'
#' @param peak_anno.df HOMER annotation data frame containing `Annotation`.
#' @param color Colour vector used for genomic feature categories.
#' @param prefix Plot title.
#'
#' @return A ggpubr pie-chart object.
#' @export
peak_pieplot_byhomer <- function(
  peak_anno.df,
  color = brewer.pal(12, "Paired")[c(1, 2, 3, 4, 5, 7, 8, 9)],
  prefix = "pieplot"
) {
  peak_anno.df$Type <- (sapply(peak_anno.df$Annotation, function(x) {
    strsplit(x, " \\(")[[1]][1]
  }))
  df.pie <- peak_anno.df %>%
    group_by(Type) %>%
    summarise(Num = n(), Pro = round(n() / nrow(peak_anno.df) * 100, 1))
  df.pie$Type <- factor(
    df.pie$Type,
    levels = c(
      "Intergenic",
      "intron",
      "exon",
      "promoter-TSS",
      "TTS",
      "non-coding",
      "3' UTR",
      "5' UTR"
    )
  )
  df.pie$labs <- paste(df.pie$Pro, "%")
  df.pie$TypeII <- paste0(df.pie$Type, "(", df.pie$labs, ")")
  df.pie$TypeII <- fct_reorder(df.pie$TypeII, -df.pie$Num)
  #ggdonutchart(df.pie, "Num", label = "labs", fill = "Type",color = "white",lab.pos = "out", lab.font = c("black","blod",5)) + theme(legend.position = "right")  + scale_fill_manual(values = color) + ggtitle(prefix)
  ggpie(
    df.pie,
    "Num",
    label = "labs",
    fill = "Type",
    color = "white",
    lab.pos = "out",
    repel = TRUE,
    lab.font = c("black", "blod", 5)
  ) +
    theme(legend.position = "right") +
    scale_fill_manual(values = color) +
    ggtitle(prefix)
}


#' Export an ordered HOMER signal heatmap
#'
#' Orders a peak-centred signal matrix by gene or region identifiers, applies a
#' `log2(x + 1)` transformation and upper cap, and writes a heatmap and separate
#' legend image.
#'
#' @param heatmap.df Data frame whose first column is `Gene` and remaining
#'   columns are ordered signal bins.
#' @param order.list Character vector defining row order.
#' @param color High-value heatmap colour.
#' @param max Maximum transformed signal value.
#' @param prefix Filename and legend title prefix.
#' @param filedir Output directory.
#' @param width,height Output heatmap dimensions in inches.
#' @param gaps Row positions at which gaps are drawn.
#'
#' @return Invisibly returns the result of the final [pheatmap::pheatmap()] call;
#'   heatmap files are written to `filedir`.
#' @export
peak_heatmap_byhomer <- function(
  heatmap.df,
  order.list,
  color = "#08537C",
  max = 10,
  prefix = "Heammap",
  filedir = ".",
  width = 2,
  height = 10,
  gaps = 0
) {
  heatmap.df.order <- heatmap.df[match(order.list, heatmap.df$Gene), ]
  heatmap.df.order <- as.matrix(heatmap.df.order[, 2:ncol(heatmap.df.order)])
  heatmap.df.order <- log2(heatmap.df.order + 1)
  heatmap.df.order[heatmap.df.order > max] <- max
  pheatmap::pheatmap(
    heatmap.df.order,
    kmeans_k = NA,
    scale = "none",
    cellwidth = NA,
    cellheight = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_names_col = FALSE,
    annotation_legend = TRUE,
    cluster_rows = FALSE,
    legend = FALSE,
    gaps_row = gaps,
    filename = paste0(filedir, "/ChIP.heatmap.", prefix, ".png"),
    width = width,
    height = height,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", color))(100)
  )
  legend_dummy <- matrix(seq(0, max, length.out = 100), nrow = 1)
  pheatmap::pheatmap(
    legend_dummy,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_legend = FALSE,
    scale = "none",
    color = colorRampPalette(c("white", color))(100),
    filename = paste0(filedir, "/ChIP.heatmap.", prefix, ".pdf"),
    width = 3,
    height = 3,
    main = prefix
  )
}


#' Visualize regulatory mutations with gene tracks
#'
#' Constructs human hg38 gene tracks for two genes, overlays mutation
#' positions, and draws the selected regulatory interval with guide lines.
#'
#' @param chr Chromosome name such as `"chr11"`.
#' @param region Numeric coordinates defining the displayed interval and guide
#'   lines; the first and fourth values set the track range.
#' @param mut_pos Numeric mutation positions.
#' @param tf_name Character vector of two human gene symbols.
#' @param num Reserved numeric mutation-count input retained by the current
#'   interface.
#'
#' @return The value returned by [trackViewer::addGuideLine()], invisibly used
#'   for its plotting side effect.
#' @export
track_view_cre_mut <- function(chr, region, mut_pos, tf_name, num) {
  #  chr <- "chr11"
  #  region <- c(32387775,32418917,32426766,32458769)
  #  mut_pos <- c(32421395,32421396,32421397,32426110)
  #  tf_name <- c("WT1","WT1-AS")
  #  num <- c(1,1,1,1)
  for (pkg in c(
    "trackViewer",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "GenomicRanges",
    "IRanges"
  )) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for track_view_cre_mut().")
    }
  }
  symbol_map <- getExportedValue("org.Hs.eg.db", "org.Hs.egSYMBOL2EG")
  txdb <- getExportedValue(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "TxDb.Hsapiens.UCSC.hg38.knownGene"
  )
  SNPs <- GenomicRanges::GRanges(
    chr,
    IRanges::IRanges(mut_pos, width = 1),
    strand = "-"
  )
  SNPs$color <- "#D7191C"
  SNPs$border <- "#D7191C"
  SNPs$feature.height <- 0.05
  SNPs$cex <- 0.8
  TF1 <- trackViewer::geneTrack(
    get(tf_name[1], symbol_map),
    txdb
  )[[1]]
  TF1$dat2 <- SNPs
  TF2 <- trackViewer::geneTrack(
    get(tf_name[2], symbol_map),
    txdb
  )[[1]]
  TF2$dat2 <- SNPs
  gr <- GenomicRanges::GRanges(
    chr,
    IRanges::IRanges(region[1], region[4])
  )
  optSty <- trackViewer::optimizeStyle(
    trackViewer::trackList(TF2, TF1),
    theme = "bw"
  )
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  trackViewer::setTrackStyleParam(
    trackList[[1]], "ylabgp", list(cex = 0.8)
  )
  trackViewer::setTrackStyleParam(
    trackList[[2]], "ylabgp", list(cex = 0.8)
  )
  names(trackList) <- rev(tf_name)
  vp <- trackViewer::viewTracks(trackList, gr = gr, viewerStyle = viewerStyle)
  trackViewer::addGuideLine(region, vp = vp)
  #  return(vp)
}


#' Plot genomic signal tracks around a region of interest
#'
#' Parses a genomic interval, adds a human hg38 gene track, imports BAM, BigWig,
#' or BED signals listed in a metadata table, applies colours and optional axis
#' limits, and draws the combined tracks.
#'
#' @param chrom Region string such as `"chr21:34,787,801-36,004,667"`.
#' @param tf_name Gene symbol used for the annotation track.
#' @param ylim Zero for automatic limits, one shared positive limit, or a vector
#'   of per-signal limits.
#' @param extend Number of bases added to both sides of the interval.
#' @param info Data frame containing `ID` and one of `BAM_t`, `BW`, or `BED`,
#'   according to `type`.
#' @param type Input signal format: `"BAM"`, `"BigWig"`, or `"BED"`.
#' @param color Colours interpolated across signal tracks.
#'
#' @return The value returned by [trackViewer::addGuideLine()], used primarily
#'   for its plotting side effect.
#' @export
trackview_peak_roi <- function(
  chrom,
  tf_name,
  ylim = 0,
  extend = 1500,
  info,
  type = "BigWig",
  color = c(
    "#08537C",
    "#A63603",
    "#A81E2C",
    "#005824",
    "#08537C",
    "#7A0177",
    "#023858"
  )
) {
  #tf_name <- "ZMIZ1"
  #info <ID BW BAM>
  #chrom <- "chr21:34,787,801-36,004,667"
  for (pkg in c(
    "trackViewer",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "GenomicRanges",
    "IRanges"
  )) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for trackview_peak_roi().")
    }
  }
  chr <- sub(":.*", "", chrom)
  region_str <- sub(".*:", "", chrom)
  region_split <- strsplit(region_str, "-")[[1]]
  region <- as.numeric(gsub(",", "", region_split))

  symbol_map <- getExportedValue("org.Hs.eg.db", "org.Hs.egSYMBOL2EG")
  txdb <- getExportedValue(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "TxDb.Hsapiens.UCSC.hg38.knownGene"
  )
  TF1 <- trackViewer::geneTrack(
    get(tf_name[1], symbol_map),
    txdb
  )[[1]]
  #chr_tf <- as.character(unique(seqnames(TF1@dat)))
  #start_tf <- min(start(TF1@dat))
  #end_tf <-  max(start(TF1@dat))

  gr <- GenomicRanges::GRanges(
    chr,
    IRanges::IRanges(min(region) - extend, max(region) + extend)
  )
  if (type == "BAM") {
    Score <- lapply(1:nrow(info), function(x) {
      trackViewer::importBam(info$BAM_t[x], ranges = gr, pairs = FALSE)
    })
  } else if (type == "BigWig") {
    Score <- lapply(1:nrow(info), function(x) {
      trackViewer::importScore(info$BW[x], format = "BigWig", ranges = gr)
    })
  } else if (type == "BED") {
    Score <- lapply(1:nrow(info), function(x) {
      trackViewer::importScore(info$BED[x], format = "BED", ranges = gr)
    })
  } else {
    stop("Unsupported type: ", type, ". Expected 'BAM' or 'BigWig'.")
  }
  optSty <- trackViewer::optimizeStyle(
    do.call(trackViewer::trackList, c(list(TF1), Score)),
    theme = "bw"
  )
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  Color <- colorRampPalette(color)(nrow(info))
  Color <- c("black", Color)

  if (length(ylim) > 1) {
    for (i in 2:length(trackList)) {
      trackViewer::setTrackStyleParam(
        trackList[[i]], "ylim", c(0, ylim[i - 1])
      )
    }
  } else if (length(ylim) == 1) {
    if (ylim > 0) {
      for (i in 2:length(trackList)) {
        trackViewer::setTrackStyleParam(trackList[[i]], "ylim", c(0, ylim))
      }
    }
  }

  for (i in 1:length(trackList)) {
    trackViewer::setTrackStyleParam(trackList[[i]], "color", Color[i])
  }
  for (i in 1:length(trackList)) {
    trackViewer::setTrackStyleParam(trackList[[i]], "ylabgp", list(cex = .8))
  }
  trackViewer::setTrackStyleParam(trackList[[1]], "height", 0.03)
  names(trackList) <- c(tf_name, info$ID)
  vp <- trackViewer::viewTracks(trackList, gr = gr, viewerStyle = viewerStyle)
  trackViewer::addGuideLine(region, vp = vp)
}


#' Plot H3K27ac signal across a regulatory element
#'
#' Reads unfolded bedGraph-like signal files for multiple samples, caps extreme
#' scores, calculates a median profile, and draws overlapping signal areas with
#' selected genomic guide lines.
#'
#' @param cre_name Regulatory-element identifier and input subdirectory name.
#' @param maxy Maximum displayed signal value.
#' @param h3k27ac.list Character vector of sample names.
#' @param color Colour vector used for sample signal areas.
#' @param dir Parent directory containing the regulatory-element subdirectory.
#' @param region Numeric genomic positions marked by vertical lines.
#'
#' @return A list containing the ggplot object (`plot`) and combined signal data
#'   frame (`df`).
#' @export
cre_h3k27ac_area_plot <- function(
  cre_name,
  maxy,
  h3k27ac.list,
  color,
  dir,
  region
) {
  #  cre_name <- "CRE_8366_WT1"
  #  maxy <- 45
  #  col <- "red"
  # region <- c()
  result <- list()
  tmp.ls <- lapply(h3k27ac.list, function(x) {
    tmp.bdg <- read.table(
      paste0(dir, cre_name, "/", x, "_signal.", cre_name, ".unfold.bed"),
      header = FALSE
    )
    tmp.bdg$Group <- x
    tmp.bdg
  })
  h3k27ac.tmp.df <- do.call(rbind, tmp.ls)
  colnames(h3k27ac.tmp.df) <- c("Chr", "Start", "End", "Score", "Group")
  h3k27ac.tmp.df[h3k27ac.tmp.df$Score > maxy, ]$Score <- maxy
  line.plot <- h3k27ac.tmp.df %>%
    group_by(Start) %>%
    summarise(Score = median(Score))
  line.plot$Group <- h3k27ac.list[1]
  p1 <- ggplot(h3k27ac.tmp.df, aes(x = Start, y = Score, fill = Group)) +
    geom_area(alpha = 2 / 16, show.legend = FALSE, position = "identity") +
    blank +
    ylim(0, maxy) +
    scale_fill_manual(values = rep(color, 16)) +
    geom_line(data = line.plot, aes(x = Start, y = Score), size = 0.3) +
    geom_vline(xintercept = region, linetype = 2) +
    ggtitle(cre_name)
  result[["plot"]] <- p1
  result[["df"]] <- h3k27ac.tmp.df
  return(result)
}


#' Sort a peak-centred heatmap by central signal density
#'
#' Calculates row sums across the central half of a binned signal matrix and
#' orders rows from highest to lowest central density.
#'
#' @param heatmap.df Data frame containing signal-bin columns and one trailing
#'   annotation column.
#'
#' @return The reordered data frame with an added `Center_density` column.
#' @export
heatmap_sort <- function(heatmap.df) {
  #heatmap.df<-heatmap.df.1
  bin_num <- (ncol(heatmap.df) - 1)
  summit <- ceiling(bin_num / 2) + 1
  start <- summit - ceiling(bin_num * 0.5) + 1
  end <- summit + ceiling(bin_num * 0.5) - 1
  colnames(heatmap.df)[start]
  colnames(heatmap.df)[summit]
  colnames(heatmap.df)[end]
  heatmap.df$Center_density <- rowSums(heatmap.df[, start:end])
  heatmap.df.order <- heatmap.df %>% dplyr::arrange(-Center_density)
}


#' Export a combined ChIP-seq signal heatmap
#'
#' Reads three peak-centred signal matrices, sorts each by central density,
#' combines them with row gaps, transforms and caps signal values, and exports a
#' single heatmap image.
#'
#' @param file.input Character vector of exactly three tab-delimited matrix
#'   paths.
#' @param brew.color Name of an RColorBrewer palette.
#' @param max Maximum `log2(x + 1)` signal displayed.
#' @param filename Output image filename.
#' @param width,height Output dimensions in inches.
#'
#' @return Invisibly returns the [pheatmap::pheatmap()] result; the heatmap is
#'   written to `filename`.
#' @export
ChIPseq_heatmap_plot <- function(
  file.input,
  brew.color = "Oranges",
  max = 10,
  filename = "ChIPseq_heatmap.png",
  width = 2,
  height = 10
) {
  heatmap.df.1 <- read.table(
    file.input[1],
    header = TRUE,
    sep = "\t",
    quote = ""
  )
  heatmap.df.2 <- read.table(
    file.input[2],
    header = TRUE,
    sep = "\t",
    quote = ""
  )
  heatmap.df.3 <- read.table(
    file.input[3],
    header = TRUE,
    sep = "\t",
    quote = ""
  )

  heatmap.df.1 <- heatmap_sort(heatmap.df.1)
  heatmap.df.2 <- heatmap_sort(heatmap.df.2)
  heatmap.df.3 <- heatmap_sort(heatmap.df.3)
  heatmap.df.order <- rbind(heatmap.df.1, heatmap.df.2, heatmap.df.3)
  bin_num <- (ncol(heatmap.df.order) - 2)
  heatmap.df.order <- as.matrix(heatmap.df.order[, 1:bin_num + 1])
  heatmap.df.order <- log2(heatmap.df.order + 1)
  heatmap.df.order[heatmap.df.order > max] <- max
  pheatmap::pheatmap(
    heatmap.df.order,
    kmeans_k = NA,
    scale = "none",
    cellwidth = NA,
    cellheight = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_names_col = FALSE,
    annotation_legend = TRUE,
    cluster_rows = FALSE,
    legend = FALSE,
    gaps_row = c(nrow(heatmap.df.1), nrow(heatmap.df.1) + nrow(heatmap.df.2)),
    filename = filename,
    width = width,
    height = height,
    cluster_cols = FALSE,
    color = colorRampPalette(brewer.pal(n = 9, name = brew.color))(100)
  )
  #dev.off()
}
