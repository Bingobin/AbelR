# Functions for enrichment analyses.

enrich_combind <- function(gene, pvc = 1, qvc = 1, universe = NULL) {
  wp2gene <- read.gmt.wp(
    "/lustre/home/acct-medkkw/medlyb/wl_proj/WL234_Lib/database/gmt/wikipathways-20190510-gmt-Homo_sapiens.gmt"
  )
  #  wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
  wpid2name <- wp2gene %>% dplyr::select(wpid, name)
  ego_bp <- enrichGO(
    gene = gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe
  )
  ego_cc <- enrichGO(
    gene = gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe
  )
  ego_mf <- enrichGO(
    gene = gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe
  )
  ewp <- enricher(
    gene,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name,
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    universe = universe
  )
  ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  ekg <- enrichKEGG(
    gene,
    organism = "hsa",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    universe = universe
  )
  ekg <- setReadable(ekg, org.Hs.eg.db, keyType = "ENTREZID")
  tmp <- list()
  tmp[["ego_bp"]] <- ego_bp
  tmp[["ego_cc"]] <- ego_cc
  tmp[["ego_mf"]] <- ego_mf
  tmp[["ewp"]] <- ewp
  tmp[["ekg"]] <- ekg
  return(tmp)
}


enrich_combind_s <- function(
  gene,
  pvc = 1,
  qvc = 1,
  universe = NULL,
  species = "human"
) {
  library(clusterProfiler)
  library(dplyr)

  if (species %in% c("human", "Homo sapiens", "hsa")) {
    library(org.Hs.eg.db)
    current_OrgDb <- org.Hs.eg.db
    kegg_org <- "hsa"
  } else if (species %in% c("mouse", "Mus musculus", "mmu")) {
    library(org.Mm.eg.db)
    current_OrgDb <- org.Mm.eg.db
    kegg_org <- "mmu"
  } else {
    stop("Species must be either 'human' or 'mouse'")
  }

  ego_bp <- enrichGO(
    gene = gene,
    OrgDb = current_OrgDb,
    keyType = 'ENTREZID',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe
  )

  ego_cc <- enrichGO(
    gene = gene,
    OrgDb = current_OrgDb,
    keyType = 'ENTREZID',
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe
  )

  ego_mf <- enrichGO(
    gene = gene,
    OrgDb = current_OrgDb,
    keyType = 'ENTREZID',
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe
  )

  ekg <- enrichKEGG(
    gene,
    organism = kegg_org,
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    universe = universe
  )
  if (!is.null(ekg)) {
    ekg <- setReadable(ekg, current_OrgDb, keyType = "ENTREZID")
  }

  tmp <- list()
  tmp[["ego_bp"]] <- ego_bp
  tmp[["ego_cc"]] <- ego_cc
  tmp[["ego_mf"]] <- ego_mf
  tmp[["ekg"]] <- ekg

  return(tmp)
}


enrich_combind_s2 <- function(
  gene,
  pvc = 1,
  qvc = 1,
  universe = NULL,
  species = "human",
  reactome = TRUE,
  hallmark = TRUE,
  minGSSize = 10,
  maxGSSize = 500
) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(dplyr)
  })

  ## -------- species setup --------
  if (species %in% c("human", "Homo sapiens", "hsa")) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    current_OrgDb <- org.Hs.eg.db
    kegg_org <- "hsa"
    reactome_org <- "human"
    msig_species <- "Homo sapiens"
  } else if (species %in% c("mouse", "Mus musculus", "mmu")) {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    current_OrgDb <- org.Mm.eg.db
    kegg_org <- "mmu"
    reactome_org <- "mouse"
    msig_species <- "Mus musculus"
  } else {
    stop("Species must be either 'human' or 'mouse'")
  }

  ## -------- helper: keep ENTREZ only --------
  gene <- unique(as.character(gene))
  gene <- gene[!is.na(gene) & gene != ""]

  ## -------- universe handling (ENTREZ universe strongly recommended) --------
  if (!is.null(universe)) {
    universe <- unique(as.character(universe))
    universe <- universe[!is.na(universe) & universe != ""]
  }

  ## -------- GO --------
  ego_bp <- enrichGO(
    gene = gene,
    OrgDb = current_OrgDb,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )

  ego_cc <- enrichGO(
    gene = gene,
    OrgDb = current_OrgDb,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )

  ego_mf <- enrichGO(
    gene = gene,
    OrgDb = current_OrgDb,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    readable = TRUE,
    universe = universe,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )

  ## -------- KEGG --------
  ekg <- enrichKEGG(
    gene = gene,
    organism = kegg_org,
    pvalueCutoff = pvc,
    qvalueCutoff = qvc,
    universe = universe,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0) {
    ekg <- setReadable(ekg, current_OrgDb, keyType = "ENTREZID")
  }

  ## -------- Reactome --------
  erct <- NULL
  if (isTRUE(reactome)) {
    suppressPackageStartupMessages(library(ReactomePA))
    erct <- enrichPathway(
      gene = gene,
      organism = reactome_org,
      pvalueCutoff = pvc,
      pAdjustMethod = "BH",
      qvalueCutoff = qvc,
      universe = universe,
      readable = TRUE,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
    if (!is.null(erct) && nrow(as.data.frame(erct)) == 0) erct <- NULL
  }

  ## -------- MSigDB Hallmark (ORA via enricher) --------
  ehall <- NULL
  if (isTRUE(hallmark)) {
    suppressPackageStartupMessages({
      library(msigdbr)
      library(tidyr)
    })

    # msigdbr returns gene symbols; convert to ENTREZID to match 'gene'
    msig_h <- msigdbr(species = msig_species, category = "H") %>%
      dplyr::select(gs_name, gene_symbol) %>%
      distinct()

    # SYMBOL -> ENTREZ
    sym2ent <- bitr(
      msig_h$gene_symbol,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = current_OrgDb
    ) %>%
      distinct(SYMBOL, ENTREZID)

    msig_h2 <- msig_h %>%
      left_join(sym2ent, by = c("gene_symbol" = "SYMBOL")) %>%
      filter(!is.na(ENTREZID)) %>%
      dplyr::select(gs_name, ENTREZID) %>%
      distinct()

    # TERM2GENE format for enricher
    term2gene_h <- msig_h2 %>% dplyr::rename(term = gs_name, gene = ENTREZID)

    ehall <- enricher(
      gene = gene,
      TERM2GENE = term2gene_h,
      universe = universe,
      pAdjustMethod = "BH",
      pvalueCutoff = pvc,
      qvalueCutoff = qvc,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )

    # add "readable" gene symbols for convenience
    if (!is.null(ehall) && nrow(as.data.frame(ehall)) > 0) {
      ehall <- setReadable(ehall, current_OrgDb, keyType = "ENTREZID")
    } else {
      ehall <- NULL
    }
  }

  ## -------- output --------
  out <- list(
    ego_bp = ego_bp,
    ego_cc = ego_cc,
    ego_mf = ego_mf,
    ekg = ekg,
    reactome = erct,
    hallmark = ehall
  )
  return(out)
}


enricher_plot <- function(
  enricher,
  bp = 1:5,
  cc = 1:5,
  mf = 1:5,
  wp = 1:5,
  kg = 1:5,
  value = 6
) {
  df <- data.frame(
    Description = enricher$ego_bp@result[bp, 2],
    Pvalue = enricher$ego_bp@result[bp, value],
    Type = rep("GO_BP", length(bp))
  )
  df.tmp <- data.frame(
    Description = enricher$ego_cc@result[cc, 2],
    Pvalue = enricher$ego_cc@result[cc, value],
    Type = rep("GO_CC", length(cc))
  )
  df <- rbind(df, df.tmp)
  df.tmp <- data.frame(
    Description = enricher$ego_mf@result[mf, 2],
    Pvalue = enricher$ego_mf@result[mf, value],
    Type = rep("GO_MF", length(mf))
  )
  df <- rbind(df, df.tmp)
  df.tmp <- data.frame(
    Description = enricher$ewp@result[wp, 2],
    Pvalue = enricher$ewp@result[wp, value],
    Type = rep("WikiPath", length(wp))
  )
  df <- rbind(df, df.tmp)
  df.tmp <- data.frame(
    Description = enricher$ekg@result[kg, 2],
    Pvalue = enricher$ekg@result[kg, value],
    Type = rep("KEGG", length(kg))
  )
  df <- rbind(df, df.tmp)
  #ggplot(df,aes(x=Description,y=-log(Pvalue),fill=Type)) + geom_bar(stat = "identity",show.legend = TRUE) + coord_flip()
  p <- ggplot(
    df,
    aes(-log10(Pvalue), fct_reorder(Description, -log10(Pvalue)))
  ) +
    geom_segment(
      aes(xend = 0, yend = Description, color = Type),
      linetype = 2,
      show.legend = FALSE
    ) +
    geom_point(aes(color = Type), size = 5, show.legend = FALSE) +
    scale_color_manual(values = brewer.pal(5, "Dark2")) +
    facet_grid(Type ~ ., scales = 'free', space = 'free_y', switch = "x") +
    blank +
    ylab("")
  return(p)
}


gsea_plot_custorm <- function(gsea_ob, select_term, color, xpos = 3000) {
  library("ggrastr")
  #gsea_ob <- aml_phenolyzer.gsea.crc
  #select_term <- 1
  #color <- "#08537C"
  nes <- round(gsea_ob@result[select_term, "NES"], digits = 2)
  pv <- formatC(
    gsea_ob@result[select_term, "p.adjust"],
    format = "e",
    digits = 2
  )
  #pv <- round(gsea_ob@result[select_term,"p.adjust"], digits = 6)
  #pv <- round(gsea_ob@result[select_term,"pvalue"], digits = 6)
  gsdata <- do.call(
    rbind,
    lapply(select_term, enrichplot:::gsInfo, object = gsea_ob)
  )
  ypos <- sign(nes) * max(abs(gsdata$runningScore)) / 2
  title_text <- gsea_ob@result[select_term, "Description"]
  title_text <- gsub("^HALLMARK[_ ]*", "", title_text, ignore.case = TRUE)
  title_text <- gsub("_+", " ", title_text)
  title_text <- tools::toTitleCase(tolower(title_text))

  thin_step <- 5
  gsdata_thin <- gsdata %>%
    dplyr::group_by(Description) %>%
    dplyr::slice(seq(1, dplyr::n(), by = thin_step)) %>%
    dplyr::ungroup()

  p_gsea_1 <- ggplot(gsdata_thin, aes(x = x, y = runningScore)) +
    geom_line(aes(color = Description), size = 1, show.legend = FALSE) +
    blank +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = color) +
    #theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    xlab("") +
    ylab("Running enrichment score") +
    annotate("text", label = paste0("NES = ", nes), x = xpos, y = ypos) +
    annotate(
      "text",
      label = paste0("p.adjust = ", pv),
      x = xpos,
      y = ypos - 0.05
    ) +
    #annotate("text", label = paste0("P value = ", pv), x=xpos,y= ypos -0.05) +
    theme(plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"))

  p_gsea_2 <- ggplot(gsdata, aes_(x = ~x)) +
    rasterise(
      geom_linerange(
        aes_(ymin = ~ymin, ymax = ~ymax),
        color = color,
        show.legend = FALSE,
        alpha = 0.6,
        size = 0.4
      ),
      dpi = 300
    ) +
    #geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax),color=color, show.legend = FALSE,alpha =0.6,size = 0.4) +
    blank +
    xlab(NULL) +
    theme(axis.ticks = element_blank(), axis.text = element_blank()) +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    ) +
    #geom_hline(yintercept = 0) +
    theme(plot.margin = margin(t = -0.5, r = 0.2, b = 0, l = 1.2, unit = "cm"))

  p_gsea_3 <- ggplot(gsdata, aes_(x = ~x, y = ~geneList)) +
    #rasterise(geom_segment(aes_(xend = ~x, yend = 0), color = color, show.legend = FALSE)) +
    rasterise(
      geom_area(color = color, fill = color, show.legend = FALSE),
      dpi = 300
    ) +
    #geom_area(color = color, fill = color, show.legend = FALSE) +
    #  scale_colour_gradient(low= brewer.pal(9,"Blues")[6], high =  brewer.pal(9,"Blues")[9]) +
    ylab("Ranked list matric") +
    xlab("Rank in ordered dataset") +
    # scale_y_continuous(n.breaks = 3)+
    blank +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    ) +
    theme(
      plot.margin = margin(t = -0.05, r = 0.2, b = 0.2, l = 0.2, unit = "cm")
    )

  gsea_body <- plot_grid(
    p_gsea_1,
    p_gsea_2,
    p_gsea_3,
    nrow = 3,
    rel_heights = c(8, 0.8, 4),
    align = "v"
  )
  #plot_grid(p_gsea_1, p_gsea_2, nrow= 2, rel_heights = c(8,1), align = "v")

  final_plot <- ggdraw() +
    draw_label(
      title_text,
      x = 0.5,
      y = 0.98,
      hjust = 0.5,
      vjust = 1,
      fontface = "bold",
      size = 12
    ) +
    draw_plot(gsea_body, y = 0, height = 0.94)

  return(final_plot)
}


go_plot_custom <- function(Enricher, select_term, color, Label, pmin) {
  #  Enricher <- Enricher.crc_tg$ego_bp
  #  select_term <- 1:20
  #  color <- "Blues"
  #  Label <- "GO_BP"
  Enricher$GeneCounts <- as.numeric(as.data.frame(
    strsplit(Enricher$GeneRatio, "/"),
    stringsAsFactors = FALSE
  )[1, ])
  plot.df <- Enricher[select_term, c("Description", "p.adjust", "GeneCounts")]
  plot.df$Label <- Label
  plot.df[plot.df$p.adjust < pmin, ]$p.adjust <- pmin
  ggplot(
    plot.df,
    aes(y = fct_reorder(Description, -log10(p.adjust)), x = Label)
  ) +
    geom_point(aes(size = GeneCounts, color = -log10(p.adjust))) +
    blank +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    ) +
    #  theme(legend.position = "bottom") +
    ylab("") +
    xlab("") +
    scale_colour_gradient(
      low = brewer.pal(9, color)[6],
      high = brewer.pal(9, color)[9]
    )
}


GO_BP_treeplot_DESeq2 <- function(deg.df, label, pv = 0.05, lfc = log2(1.5)) {
  #  deg.df <- HMGA2_sh2.DESeq2.result$result
  #  label <- "shHMGA2"
  up.entrez <- (deg.df %>% filter(padj < pv, log2FoldChange > lfc))$Entrez
  dw.entrez <- (deg.df %>% filter(padj < pv, log2FoldChange < -lfc))$Entrez
  up.entrez <- unique(na.omit(up.entrez))
  dw.entrez <- unique(na.omit(dw.entrez))
  deg.entrez <- unique(c(up.entrez, dw.entrez))
  color <- "YlGnBu"
  up.enricher <- enrich_combind(up.entrez)
  dw.enricher <- enrich_combind(dw.entrez)
  deg.enricher <- enrich_combind(deg.entrez)
  ##################go_bp
  ego_bp.up2 <- enrichplot::pairwise_termsim(up.enricher$ego_bp)
  ego_bp.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ego_bp)
  ego_bp.deg2 <- enrichplot::pairwise_termsim(deg.enricher$ego_bp)

  p_up_bp <- enrichplot::treeplot(
    ego_bp.up2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Up-regulated Genes in GO_BP (',
        label,
        ' n= ',
        length(up.entrez),
        ')'
      )
    )
  p_dw_bp <- enrichplot::treeplot(
    ego_bp.dw2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Down-regulated Genes in GO_BP (',
        label,
        ' n= ',
        length(dw.entrez),
        ')'
      )
    )
  p_deg_bp <- enrichplot::treeplot(
    ego_bp.deg2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Diff. Expr. Genes in GO_BP (',
        label,
        ' n= ',
        length(deg.entrez),
        ')'
      )
    )
  ##################kegg
  ekg.up2 <- enrichplot::pairwise_termsim(up.enricher$ekg)
  ekg.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ekg)
  ekg.deg2 <- enrichplot::pairwise_termsim(deg.enricher$ekg)

  p_up_kegg <- enrichplot::treeplot(
    ekg.up2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Up-regulated Genes in KEGG pathway (',
        label,
        ' n= ',
        length(up.entrez),
        ')'
      )
    )
  p_dw_kegg <- enrichplot::treeplot(
    ekg.dw2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Down-regulated Genes in KEGG pathway (',
        label,
        ' n= ',
        length(dw.entrez),
        ')'
      )
    )
  p_deg_kegg <- enrichplot::treeplot(
    ekg.deg2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Diff. Expr. Genes in KEGG pathway (',
        label,
        ' n= ',
        length(deg.entrez),
        ')'
      )
    )

  ###################

  p_combind_up <- enricher_plot(up.enricher) +
    ggtitle(
      label = paste0(
        'Up-regulated Genes Enrichment Analysis(',
        label,
        ' n= ',
        length(up.entrez),
        ')'
      )
    )
  p_combind_dw <- enricher_plot(dw.enricher) +
    ggtitle(
      label = paste0(
        'Down-regulated Genes Enrichment Analysis(',
        label,
        ' n= ',
        length(dw.entrez),
        ')'
      )
    )
  p_combind_deg <- enricher_plot(deg.enricher) +
    ggtitle(
      label = paste0(
        'Diff. Expr. Genes Enrichment Analysis(',
        label,
        ' n= ',
        length(deg.entrez),
        ')'
      )
    )

  plot.ls <- list()
  plot.ls[["ego_bp_up"]] <- p_up_bp
  plot.ls[["ego_bp_dw"]] <- p_dw_bp
  plot.ls[["ego_bp_deg"]] <- p_deg_bp
  plot.ls[["ekg_up"]] <- p_up_kegg
  plot.ls[["ekg_dw"]] <- p_dw_kegg
  plot.ls[["ekg_deg"]] <- p_deg_kegg
  plot.ls[["combind_up"]] <- p_combind_up
  plot.ls[["combind_dw"]] <- p_combind_dw
  plot.ls[["combind_deg"]] <- p_combind_deg

  result <- list()
  result[["label"]] <- label
  result[["up_num"]] <- length(up.entrez)
  result[["dw_num"]] <- length(dw.entrez)
  result[["deg_num"]] <- length(deg.entrez)
  result[["up_enricher"]] <- up.enricher
  result[["dw_enricher"]] <- dw.enricher
  result[["deg_enricher"]] <- deg.enricher
  result[["plot"]] <- plot.ls

  return(result)
}


GO_BP_treeplot_scRNAseq <- function(deg.df, label, pv = 0.05, lfc = 0.25) {
  #  deg.df <- XGJ.NEMOBvsNEMOA.degs
  #  label <- "NEMOBvsNEMOA"
  deg.df$Entrez <- TransGeneID(
    rownames(deg.df),
    "Symbol",
    "Entrez",
    organism = "hsa"
  )
  up.entrez <- (deg.df %>% filter(p_val < pv, avg_log2FC > lfc))$Entrez
  dw.entrez <- (deg.df %>% filter(p_val < pv, avg_log2FC < -lfc))$Entrez
  up.entrez <- unique(na.omit(up.entrez))
  dw.entrez <- unique(na.omit(dw.entrez))
  color <- "YlGnBu"
  up.enricher <- enrich_combind(up.entrez)
  dw.enricher <- enrich_combind(dw.entrez)
  ##################go_bp
  ego_bp.up2 <- enrichplot::pairwise_termsim(up.enricher$ego_bp)
  ego_bp.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ego_bp)
  p_up_bp <- enrichplot::treeplot(
    ego_bp.up2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Up-regulated Genes in GO_BP (',
        label,
        ' n= ',
        length(up.entrez),
        ')'
      )
    )
  p_dw_bp <- enrichplot::treeplot(
    ego_bp.dw2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Down-regulated Genes in GO_BP (',
        label,
        ' n= ',
        length(dw.entrez),
        ')'
      )
    )
  ##################kegg
  ekg.up2 <- enrichplot::pairwise_termsim(up.enricher$ekg)
  ekg.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ekg)
  p_up_kegg <- enrichplot::treeplot(
    ekg.up2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Up-regulated Genes in KEGG pathway (',
        label,
        ' n= ',
        length(up.entrez),
        ')'
      )
    )
  p_dw_kegg <- enrichplot::treeplot(
    ekg.dw2,
    nCluster = 5,
    group_color = brewer.pal(5, "Set1"),
    hclust_methd = "average",
    offset = 10
  ) +
    scale_colour_gradient(
      high = brewer.pal(9, color)[6],
      low = brewer.pal(9, color)[9]
    ) +
    ggtitle(
      label = paste0(
        'Down-regulated Genes in KEGG pathway (',
        label,
        ' n= ',
        length(dw.entrez),
        ')'
      )
    )
  ###################

  p_combind_up <- enricher_plot(up.enricher) +
    ggtitle(
      label = paste0(
        'Up-regulated Genes Enrichment Analysis(',
        label,
        ' n= ',
        length(up.entrez),
        ')'
      )
    )
  p_combind_dw <- enricher_plot(dw.enricher) +
    ggtitle(
      label = paste0(
        'Down-regulated Genes Enrichment Analysis(',
        label,
        ' n= ',
        length(dw.entrez),
        ')'
      )
    )
  plot.ls <- list()
  plot.ls[["ego_bp_up"]] <- p_up_bp
  plot.ls[["ego_bp_dw"]] <- p_dw_bp
  plot.ls[["ekg_up"]] <- p_up_kegg
  plot.ls[["ekg_dw"]] <- p_dw_kegg
  plot.ls[["combind_up"]] <- p_combind_up
  plot.ls[["combind_dw"]] <- p_combind_dw

  result <- list()
  result[["label"]] <- label
  result[["up_num"]] <- length(up.entrez)
  result[["dw_num"]] <- length(dw.entrez)
  result[["up_enricher"]] <- up.enricher
  result[["dw_enricher"]] <- dw.enricher
  result[["plot"]] <- plot.ls

  return(result)
}


GSEA_DEGseq <- function(deg.df, species = "human") {
  # 加载必要的包
  library(msigdbr)
  library(clusterProfiler)
  library(dplyr)

  # 1. 根据物种参数设置对应的变量
  if (species %in% c("human", "Homo sapiens", "hsa")) {
    library(org.Hs.eg.db) # 确保加载人类数据库
    msig_species <- "Homo sapiens"
    kegg_org <- "hsa"
    current_OrgDb <- org.Hs.eg.db
  } else if (species %in% c("mouse", "Mus musculus", "mmu")) {
    library(org.Mm.eg.db) # 确保加载小鼠数据库
    msig_species <- "Mus musculus"
    kegg_org <- "mmu"
    current_OrgDb <- org.Mm.eg.db
  } else {
    stop("Species must be either 'human' or 'mouse'")
  }

  # 2. 准备基因列表 (这一步不分物种，假设输入都是Entrez ID)
  # 确保 log2FoldChange 是数值型，Entrez 是字符型
  gene.fc <- deg.df[
    !is.na(deg.df$Entrez) & !duplicated(deg.df$Entrez),
  ]$log2FoldChange
  names(gene.fc) <- deg.df[
    !is.na(deg.df$Entrez) & !duplicated(deg.df$Entrez),
  ]$Entrez
  gene.fc <- sort(gene.fc, decreasing = TRUE)

  # 3. 获取 Hallmark 基因集 (使用动态变量 msig_species)
  m_t2g.H <- msigdbr(species = msig_species, category = "H") %>%
    dplyr::select(gs_name, entrez_gene)

  # 4. GSEA Hallmark 分析
  gsea.H <- GSEA(gene.fc, TERM2GENE = m_t2g.H, pvalueCutoff = 1)
  # 使用动态变量 current_OrgDb 进行 ID 转换
  if (!is.null(gsea.H)) {
    gsea.H <- setReadable(gsea.H, OrgDb = current_OrgDb, keyType = "ENTREZID")
  }

  # 5. GSEA KEGG 分析 (使用动态变量 kegg_org)
  gsea.kegg <- gseKEGG(
    geneList = gene.fc,
    organism = kegg_org,
    keyType = "kegg",
    pvalueCutoff = 1
  )
  if (!is.null(gsea.kegg)) {
    gsea.kegg <- setReadable(
      gsea.kegg,
      OrgDb = current_OrgDb,
      keyType = "ENTREZID"
    )
  }

  # 6. GSEA GO_BP 分析 (使用动态变量 current_OrgDb)
  gsea.go_bp <- gseGO(
    gene.fc,
    OrgDb = current_OrgDb,
    ont = "BP",
    pvalueCutoff = 1,
    keyType = "ENTREZID"
  )
  if (!is.null(gsea.go_bp)) {
    gsea.go_bp <- setReadable(
      gsea.go_bp,
      OrgDb = current_OrgDb,
      keyType = "ENTREZID"
    )
  }

  # 7. 返回结果
  result <- list()
  result[["Hallmark"]] <- gsea.H
  result[["KEGG"]] <- gsea.kegg
  result[["GO_BP"]] <- gsea.go_bp

  return(result)
}


GSEA_scRNAseq <- function(deg.df, species = c("human", "mouse")) {
  species <- match.arg(species)

  library(dplyr)
  library(msigdbr)
  library(clusterProfiler)
  library("org.Mm.eg.db")

  if (!"Symbol" %in% colnames(deg.df)) {
    stop("Column 'Symbol' not found in deg.df")
  }
  if (!"avg_log2FC" %in% colnames(deg.df)) {
    stop("Column 'avg_log2FC' not found in deg.df")
  }

  ## species-specific settings
  if (species == "human") {
    msig_species <- "Homo sapiens"
    org_db <- org.Hs.eg.db
    kegg_org <- "hsa"
  } else if (species == "mouse") {
    msig_species <- "Mus musculus"
    org_db <- org.Mm.eg.db
    kegg_org <- "mmu"
  }

  ## Hallmark gene sets
  m_t2g.H <- msigdbr(species = msig_species, category = "H") %>%
    dplyr::select(gs_name, entrez_gene)

  ## Symbol -> Entrez
  deg.df$Entrez <- TransGeneID(
    deg.df$Symbol,
    fromType = "Symbol",
    toType = "Entrez",
    organism = ifelse(species == "human", "hsa", "mmu")
  )

  ## build ranked gene list
  deg.sub <- deg.df[!is.na(deg.df$Entrez), c("Entrez", "avg_log2FC")]
  deg.sub$Entrez <- as.character(deg.sub$Entrez)

  ## remove duplicated Entrez IDs, keep the one with largest absolute FC
  deg.sub <- deg.sub %>%
    group_by(Entrez) %>%
    slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
    ungroup()

  gene.fc <- deg.sub$avg_log2FC
  names(gene.fc) <- deg.sub$Entrez
  gene.fc <- sort(gene.fc, decreasing = TRUE)

  ## run GSEA
  gsea.H <- GSEA(gene.fc, TERM2GENE = m_t2g.H, pvalueCutoff = 1)
  gsea.H <- setReadable(gsea.H, OrgDb = org_db, keyType = "ENTREZID")

  gsea.kegg <- gseKEGG(
    geneList = gene.fc,
    organism = kegg_org,
    keyType = "kegg",
    pvalueCutoff = 1
  )
  gsea.kegg <- setReadable(gsea.kegg, OrgDb = org_db, keyType = "ENTREZID")

  gsea.go_bp <- gseGO(
    geneList = gene.fc,
    OrgDb = org_db,
    ont = "BP",
    pvalueCutoff = 1,
    keyType = "ENTREZID"
  )
  gsea.go_bp <- setReadable(gsea.go_bp, OrgDb = org_db, keyType = "ENTREZID")

  result <- list()
  result[["Hallmark"]] <- gsea.H
  result[["KEGG"]] <- gsea.kegg
  result[["GO_BP"]] <- gsea.go_bp

  return(result)
}


GSEA_analysis <- function(
  deg.df,
  species = "human",
  rank_by = "log2FC",
  dbs = c("Hallmark", "GO_BP"),
  gene_id_col = "Entrez",
  log2fc_col = NULL,
  stat_col = "stat",
  p_col = NULL,
  padj_col = NULL,
  collapse_dup = c("max_abs", "mean", "first"),
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  p_floor = 1e-300,
  seed = 123,
  kegg_source = c("clusterProfiler", "msigdbr"),
  readable = TRUE,
  verbose = FALSE
) {
  collapse_dup <- match.arg(collapse_dup)
  kegg_source <- match.arg(kegg_source)

  ## 0. check packages
  required_pkgs <- c("clusterProfiler", "msigdbr", "dplyr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package not installed: ", pkg)
    }
  }

  ## 1. species setting
  species_lower <- tolower(species)

  if (species_lower %in% c("human", "homo sapiens", "hsa")) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package not installed: org.Hs.eg.db")
    }
    msig_species <- "Homo sapiens"
    kegg_org <- "hsa"
    current_OrgDb <- org.Hs.eg.db::org.Hs.eg.db
  } else if (species_lower %in% c("mouse", "mus musculus", "mmu")) {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package not installed: org.Mm.eg.db")
    }
    msig_species <- "Mus musculus"
    kegg_org <- "mmu"
    current_OrgDb <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    stop(
      "species must be one of: human, Homo sapiens, hsa, mouse, Mus musculus, mmu"
    )
  }

  ## 2. automatically detect common column names
  available_cols <- colnames(deg.df)

  pick_col <- function(candidates) {
    hit <- candidates[candidates %in% available_cols]
    if (length(hit) == 0) {
      return(NULL)
    }
    hit[1]
  }

  if (is.null(log2fc_col)) {
    log2fc_col <- pick_col(c(
      "log2FoldChange",
      "avg_log2FC",
      "avg_logFC",
      "logFC",
      "log2FC"
    ))
  }

  if (is.null(p_col)) {
    p_col <- pick_col(c(
      "pvalue",
      "p_val",
      "P.Value",
      "p.value",
      "p"
    ))
  }

  if (is.null(padj_col)) {
    padj_col <- pick_col(c(
      "padj",
      "p_val_adj",
      "adj.P.Val",
      "FDR",
      "qvalue",
      "q_value"
    ))
  }

  get_numeric_col <- function(col, col_desc) {
    if (is.null(col) || !(col %in% colnames(deg.df))) {
      stop(
        "Cannot find ",
        col_desc,
        " column. Please specify the correct column name."
      )
    }
    suppressWarnings(as.numeric(as.character(deg.df[[col]])))
  }

  ## 3. construct ranking metric
  rank_by_lower <- tolower(rank_by)

  log2fc <- NULL
  if (!is.null(log2fc_col) && log2fc_col %in% colnames(deg.df)) {
    log2fc <- get_numeric_col(log2fc_col, "log2FC")
  }

  if (
    rank_by_lower %in% c("log2fc", "log2foldchange", "avg_log2fc", "avg_logfc")
  ) {
    metric <- get_numeric_col(log2fc_col, "log2FC")
    metric_name <- log2fc_col
  } else if (rank_by_lower %in% c("stat", "wald", "wald_stat")) {
    metric <- get_numeric_col(stat_col, "stat")
    metric_name <- stat_col
  } else if (rank_by_lower %in% c("signed_p", "signed_log10p", "sign_p")) {
    if (is.null(log2fc)) {
      stop("log2FC column is required for signed_p ranking.")
    }
    pval <- get_numeric_col(p_col, "p value")
    pval <- pmax(pval, p_floor, na.rm = TRUE)

    metric <- sign(log2fc) * (-log10(pval))
    metric_name <- paste0("sign(", log2fc_col, ") * -log10(", p_col, ")")
  } else if (
    rank_by_lower %in% c("signed_padj", "signed_fdr", "signed_log10padj")
  ) {
    if (is.null(log2fc)) {
      stop("log2FC column is required for signed_padj ranking.")
    }
    padj <- get_numeric_col(padj_col, "adjusted p value")
    padj <- pmax(padj, p_floor, na.rm = TRUE)

    metric <- sign(log2fc) * (-log10(padj))
    metric_name <- paste0("sign(", log2fc_col, ") * -log10(", padj_col, ")")
  } else if (rank_by_lower %in% c("log2fc_p", "log2fc_x_p", "fc_p")) {
    if (is.null(log2fc)) {
      stop("log2FC column is required for log2FC_p ranking.")
    }
    pval <- get_numeric_col(p_col, "p value")
    pval <- pmax(pval, p_floor, na.rm = TRUE)

    metric <- log2fc * (-log10(pval))
    metric_name <- paste0(log2fc_col, " * -log10(", p_col, ")")
  } else if (rank_by_lower %in% c("log2fc_padj", "log2fc_x_padj", "fc_fdr")) {
    if (is.null(log2fc)) {
      stop("log2FC column is required for log2FC_padj ranking.")
    }
    padj <- get_numeric_col(padj_col, "adjusted p value")
    padj <- pmax(padj, p_floor, na.rm = TRUE)

    metric <- log2fc * (-log10(padj))
    metric_name <- paste0(log2fc_col, " * -log10(", padj_col, ")")
  } else if (rank_by %in% colnames(deg.df)) {
    metric <- get_numeric_col(rank_by, rank_by)
    metric_name <- rank_by
  } else {
    stop(
      "Unknown rank_by: ",
      rank_by,
      "\n",
      "Allowed preset values: log2FC, stat, signed_p, signed_padj, log2FC_p, log2FC_padj\n",
      "Or directly provide a numeric column name in deg.df."
    )
  }

  ## 4. prepare geneList
  ## 如果没有 Entrez，但有 Symbol，则自动转换 Symbol -> Entrez
  if (!(gene_id_col %in% colnames(deg.df))) {
    if (gene_id_col == "Entrez" && "Symbol" %in% colnames(deg.df)) {
      gene.map <- clusterProfiler::bitr(
        unique(deg.df$Symbol),
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = current_OrgDb
      )

      deg.df <- deg.df %>%
        dplyr::left_join(gene.map, by = c("Symbol" = "SYMBOL"))

      deg.df$Entrez <- deg.df$ENTREZID
    } else {
      stop("Cannot find gene_id_col: ", gene_id_col)
    }
  }

  gene_id <- as.character(deg.df[[gene_id_col]])
  gene_id <- sub("\\.0$", "", gene_id)

  rank.df <- data.frame(
    gene = gene_id,
    metric = metric,
    stringsAsFactors = FALSE
  )

  rank.df <- rank.df[
    !is.na(rank.df$gene) &
      rank.df$gene != "" &
      !is.na(rank.df$metric) &
      is.finite(rank.df$metric),
  ]

  if (nrow(rank.df) < 10) {
    stop("Too few valid genes for GSEA after filtering.")
  }

  ## deal with duplicated Entrez IDs
  if (collapse_dup == "max_abs") {
    rank.df <- rank.df[order(rank.df$gene, -abs(rank.df$metric)), ]
    rank.df <- rank.df[!duplicated(rank.df$gene), ]
  } else if (collapse_dup == "mean") {
    rank.df <- stats::aggregate(
      metric ~ gene,
      data = rank.df,
      FUN = mean
    )
  } else if (collapse_dup == "first") {
    rank.df <- rank.df[!duplicated(rank.df$gene), ]
  }

  gene.fc <- rank.df$metric
  names(gene.fc) <- rank.df$gene
  gene.fc <- sort(gene.fc, decreasing = TRUE)

  set.seed(seed)

  ## 5. helper functions
  get_msig_t2g <- function(collection, subcollection = NULL) {
    args_new <- list(
      species = msig_species,
      collection = collection
    )
    if (!is.null(subcollection)) {
      args_new$subcollection <- subcollection
    }

    m <- try(
      do.call(msigdbr::msigdbr, args_new),
      silent = TRUE
    )

    if (inherits(m, "try-error")) {
      args_old <- list(
        species = msig_species,
        category = collection
      )
      if (!is.null(subcollection)) {
        args_old$subcategory <- subcollection
      }

      m <- do.call(msigdbr::msigdbr, args_old)
    }

    gene_col <- NULL
    if ("entrez_gene" %in% colnames(m)) {
      gene_col <- "entrez_gene"
    } else if ("ncbi_gene" %in% colnames(m)) {
      gene_col <- "ncbi_gene"
    } else {
      stop("Cannot find Entrez/NCBI gene column in msigdbr output.")
    }

    t2g <- m[, c("gs_name", gene_col)]
    colnames(t2g) <- c("term", "gene")
    t2g$gene <- as.character(t2g$gene)
    t2g <- unique(t2g[!is.na(t2g$gene) & t2g$gene != "", ])

    return(t2g)
  }

  safe_run <- function(label, fun) {
    tryCatch(
      fun(),
      error = function(e) {
        warning(label, " failed: ", conditionMessage(e))
        return(NULL)
      }
    )
  }

  make_readable <- function(x, label) {
    if (!readable || is.null(x)) {
      return(x)
    }

    tryCatch(
      clusterProfiler::setReadable(
        x,
        OrgDb = current_OrgDb,
        keyType = "ENTREZID"
      ),
      error = function(e) {
        warning("setReadable failed for ", label, ": ", conditionMessage(e))
        return(x)
      }
    )
  }

  run_t2g_gsea <- function(label, t2g) {
    res <- safe_run(label, function() {
      clusterProfiler::GSEA(
        geneList = gene.fc,
        TERM2GENE = t2g,
        pvalueCutoff = pvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        verbose = verbose,
        seed = TRUE
      )
    })

    res <- make_readable(res, label)
    return(res)
  }

  ## 6. run selected databases
  result <- list()

  for (db in dbs) {
    db_key <- toupper(gsub("[ .-]", "_", db))

    if (db_key %in% c("H", "HALLMARK", "MSIGDB_H")) {
      t2g <- get_msig_t2g(collection = "H")
      result[["Hallmark"]] <- run_t2g_gsea("Hallmark", t2g)
    } else if (db_key %in% c("GO_BP", "BP")) {
      res <- safe_run("GO_BP", function() {
        clusterProfiler::gseGO(
          geneList = gene.fc,
          OrgDb = current_OrgDb,
          ont = "BP",
          keyType = "ENTREZID",
          pvalueCutoff = pvalueCutoff,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          verbose = verbose,
          seed = TRUE
        )
      })

      result[["GO_BP"]] <- make_readable(res, "GO_BP")
    } else if (db_key %in% c("GO_CC", "CC")) {
      res <- safe_run("GO_CC", function() {
        clusterProfiler::gseGO(
          geneList = gene.fc,
          OrgDb = current_OrgDb,
          ont = "CC",
          keyType = "ENTREZID",
          pvalueCutoff = pvalueCutoff,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          verbose = verbose,
          seed = TRUE
        )
      })

      result[["GO_CC"]] <- make_readable(res, "GO_CC")
    } else if (db_key %in% c("GO_MF", "MF")) {
      res <- safe_run("GO_MF", function() {
        clusterProfiler::gseGO(
          geneList = gene.fc,
          OrgDb = current_OrgDb,
          ont = "MF",
          keyType = "ENTREZID",
          pvalueCutoff = pvalueCutoff,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          verbose = verbose,
          seed = TRUE
        )
      })

      result[["GO_MF"]] <- make_readable(res, "GO_MF")
    } else if (db_key %in% c("GO_ALL", "GO")) {
      for (ont in c("BP", "CC", "MF")) {
        label <- paste0("GO_", ont)

        res <- safe_run(label, function() {
          clusterProfiler::gseGO(
            geneList = gene.fc,
            OrgDb = current_OrgDb,
            ont = ont,
            keyType = "ENTREZID",
            pvalueCutoff = pvalueCutoff,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            verbose = verbose,
            seed = TRUE
          )
        })

        result[[label]] <- make_readable(res, label)
      }
    } else if (db_key %in% c("KEGG", "GSEKEGG")) {
      if (kegg_source == "clusterProfiler") {
        res <- safe_run("KEGG", function() {
          clusterProfiler::gseKEGG(
            geneList = gene.fc,
            organism = kegg_org,
            keyType = "ncbi-geneid",
            pvalueCutoff = pvalueCutoff,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            verbose = verbose
          )
        })

        result[["KEGG"]] <- make_readable(res, "KEGG")
      } else {
        t2g <- get_msig_t2g(
          collection = "C2",
          subcollection = "CP:KEGG"
        )
        result[["KEGG_msigdb"]] <- run_t2g_gsea("KEGG_msigdb", t2g)
      }
    } else if (db_key %in% c("REACTOME", "C2_CP_REACTOME", "MSIGDB_REACTOME")) {
      t2g <- get_msig_t2g(
        collection = "C2",
        subcollection = "CP:REACTOME"
      )
      result[["Reactome"]] <- run_t2g_gsea("Reactome", t2g)
    } else if (db_key %in% c("C2_CP_KEGG", "MSIGDB_KEGG")) {
      t2g <- get_msig_t2g(
        collection = "C2",
        subcollection = "CP:KEGG"
      )
      result[["KEGG_msigdb"]] <- run_t2g_gsea("KEGG_msigdb", t2g)
    } else {
      warning(
        "Unknown database: ",
        db,
        ". Allowed examples: Hallmark, GO_BP, GO_CC, GO_MF, GO_ALL, KEGG, Reactome, C2_CP_KEGG"
      )
    }
  }

  result[["geneList"]] <- gene.fc
  result[["rank_info"]] <- list(
    rank_by = rank_by,
    metric_name = metric_name,
    gene_id_col = gene_id_col,
    log2fc_col = log2fc_col,
    stat_col = stat_col,
    p_col = p_col,
    padj_col = padj_col,
    collapse_dup = collapse_dup,
    n_genes = length(gene.fc),
    species = species
  )

  return(result)
}

