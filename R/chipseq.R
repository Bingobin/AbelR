# Functions for chipseq analyses.

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
  library("ggrastr")
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
  df.pie <- data.table::as.data.table(df.pie)
  df.pie[, `:=`(TypeII, paste0(Type, "(", labs, ")"))]
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


peak_heatmap_byhomer <- function(
  heatmap.df,
  order.list,
  color = "#08537C",
  max = 10,
  prefix = "Heammap",
  filedir = "~/Desktop",
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


track_view_cre_mut <- function(chr, region, mut_pos, tf_name, num) {
  #  chr <- "chr11"
  #  region <- c(32387775,32418917,32426766,32458769)
  #  mut_pos <- c(32421395,32421396,32421397,32426110)
  #  tf_name <- c("WT1","WT1-AS")
  #  num <- c(1,1,1,1)
  SNPs <- GRanges(chr, IRanges(mut_pos, width = 1), strand = "-")
  SNPs$color <- "#D7191C"
  SNPs$border <- "#D7191C"
  SNPs$feature.height <- 0.05
  SNPs$cex <- 0.8
  TF1 <- geneTrack(
    get(tf_name[1], org.Hs.egSYMBOL2EG),
    TxDb.Hsapiens.UCSC.hg38.knownGene
  )[[1]]
  TF1$dat2 <- SNPs
  TF2 <- geneTrack(
    get(tf_name[2], org.Hs.egSYMBOL2EG),
    TxDb.Hsapiens.UCSC.hg38.knownGene
  )[[1]]
  TF2$dat2 <- SNPs
  gr <- GRanges(chr, IRanges(region[1], region[4]))
  optSty <- optimizeStyle(trackList(TF2, TF1), theme = "bw")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex = 0.8))
  setTrackStyleParam(trackList[[2]], "ylabgp", list(cex = 0.8))
  names(trackList) <- rev(tf_name)
  vp <- viewTracks(trackList, gr = gr, viewerStyle = viewerStyle)
  addGuideLine(region, vp = vp)
  #  return(vp)
}


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
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(trackViewer)
  chr <- sub(":.*", "", chrom)
  region_str <- sub(".*:", "", chrom)
  region_split <- strsplit(region_str, "-")[[1]]
  region <- as.numeric(gsub(",", "", region_split))

  TF1 <- geneTrack(
    get(tf_name[1], org.Hs.egSYMBOL2EG),
    TxDb.Hsapiens.UCSC.hg38.knownGene
  )[[1]]
  #chr_tf <- as.character(unique(seqnames(TF1@dat)))
  #start_tf <- min(start(TF1@dat))
  #end_tf <-  max(start(TF1@dat))

  gr <- GRanges(chr, IRanges(min(region) - extend, max(region) + extend))
  if (type == "BAM") {
    Score <- lapply(1:nrow(info), function(x) {
      importBam(info$BAM_t[x], ranges = gr, pairs = FALSE)
    })
  } else if (type == "BigWig") {
    Score <- lapply(1:nrow(info), function(x) {
      importScore(info$BW[x], format = "BigWig", ranges = gr)
    })
  } else if (type == "BED") {
    Score <- lapply(1:nrow(info), function(x) {
      importScore(info$BED[x], format = "BED", ranges = gr)
    })
  } else {
    stop("Unsupported type: ", type, ". Expected 'BAM' or 'BigWig'.")
  }
  optSty <- optimizeStyle(
    do.call(trackViewer::trackList, c(list(TF1), Score)),
    theme = "bw"
  )
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  Color <- colorRampPalette(color)(nrow(info))
  Color <- c("black", Color)

  if (length(ylim) > 1) {
    for (i in 2:length(trackList)) {
      setTrackStyleParam(trackList[[i]], "ylim", c(0, ylim[i - 1]))
    }
  } else if (length(ylim) == 1) {
    if (ylim > 0) {
      for (i in 2:length(trackList)) {
        setTrackStyleParam(trackList[[i]], "ylim", c(0, ylim))
      }
    }
  }

  for (i in 1:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "color", Color[i])
  }
  for (i in 1:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex = .8))
  }
  setTrackStyleParam(trackList[[1]], "height", 0.03)
  names(trackList) <- c(tf_name, info$ID)
  vp <- viewTracks(trackList, gr = gr, viewerStyle = viewerStyle)
  addGuideLine(region, vp = vp)
}


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
  #  dir <- "/lustre/home/acct-medkkw/medlyb/project/15.APL_TF_ChIPseq/08.MYC_4C/02.bdg_unfold/"
  # region <- c()
  result <- list()
  tmp.ls <- lapply(h3k27ac.list, function(x) {
    #    tmp.bdg <- read.table(paste0("/lustre/home/acct-medkkw/medlyb/project/20.APL_AC_ChIPseq/10.CRE_h3k27ac_signal/",cre_name,"/",x,"_H3K27ac.",cre_name,".unfold.bed"), header = FALSE)
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


ChIPseq_heatmap_plot <- function(
  file.input,
  brew.color = "Oranges",
  max = 10,
  filename = "~/tmp.png",
  width = 2,
  height = 10
) {
  #  file.input <- c("/lustre/home/acct-medkkw/medlyb/wl_proj/APL_H3K27ac_WangLab/09.APL_CD34_process/heatmap/CHH_H3K27ac_peaks.APL_CD34.heatmap.txt","/lustre/home/acct-medkkw/medlyb/wl_proj/APL_H3K27ac_WangLab/09.APL_CD34_process/heatmap/CHH_H3K27ac_peaks.APL_uniq.heatmap.txt","/lustre/home/acct-medkkw/medlyb/wl_proj/APL_H3K27ac_WangLab/09.APL_CD34_process/heatmap/CHH_H3K27ac_peaks.CD34_uniq.heatmap.txt")
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
  #library(pheatmap)
  heatmap.df.order <- rbind(heatmap.df.1, heatmap.df.2, heatmap.df.3)
  bin_num <- (ncol(heatmap.df.order) - 2)
  heatmap.df.order <- as.matrix(heatmap.df.order[, 1:bin_num + 1])
  heatmap.df.order <- log2(heatmap.df.order + 1)
  heatmap.df.order[heatmap.df.order > max] <- max
  #png("~/tmp.png",width = 1000,height = 5000)
  pheatmap(
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


