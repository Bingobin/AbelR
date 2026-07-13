##Part Library: Loading R package

invisible(lapply(
  c("ggplot2", "cowplot", "tidyr", "forcats", "UpSetR", "RColorBrewer", "DESeq2", "survival",
    "survminer", "survivalROC","glmnet", "org.Hs.eg.db", "MAGeCKFlute", "clusterProfiler",  "Seurat",
    "WGCNA", "dplyr", "patchwork", "tidyverse", "ggrastr"),
  function(x) suppressPackageStartupMessages(library(x, character.only = TRUE))
))


blank <- theme(panel.border = element_blank(), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.background = element_blank()
               )
#+ theme(panel.border = element_rect(fill = NA, linetype = 1, size = 1), axis.line = element_blank())

ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

mycol <- c("#3C5488","#A81E2C","#B5AA0F","#7BA39D","#E64B35","#4DBBD5","#C81E76","#00A087","#8491B4","#B8AC90", "#E1882F","#0098C9", "#128640","#7E6148","#748098","#522973","#cab2d6","#e5486e","#fb9a99")


##Part Subplement: Build Function

DESeq2_DEG_analysis <-function(count.matrix,tpm.matrix = NULL,design.df,tr,ctr){
  #SampleID Group
  
  stopifnot(all(c("Group") %in% colnames(design.df)))
  
  design <- design.df |> filter(Group %in% c(ctr,tr))
  design$Group <- factor(design$Group,levels = c(ctr,tr))
  design$SampleID <- rownames(design)
  counts <- count.matrix[,rownames(design)]
  rownames(counts) <- count.matrix$GID
  tpm.mat <- tpm.matrix[,rownames(design)]
  rownames(tpm.mat) <- tpm.matrix$GID
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design =~ Group)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  inter<-resultsNames(dds)
  resOrdered = res[order(res$padj,na.last = TRUE),]
  DESeq2::plotMA(resOrdered, alpha = 0.05)
  deseq2_result <- data.frame(resOrdered)
  gene_id_anno <- read.table("~/database/GENCODE/gene_len.v43.new.txt", header = TRUE, sep = "\t",stringsAsFactors=FALSE,row.names = 1)
  gene_id_anno <- gene_id_anno[!grepl("_PAR_", rownames(gene_id_anno)),]
  
  rownames(gene_id_anno) <- sub("\\..*$","",rownames(gene_id_anno))
  deseq2_result <- merge(gene_id_anno, deseq2_result, by.x=0, by.y=0, all = F)
  deseq2_result <- merge(deseq2_result,tpm.mat,by.x = 1, by.y=0)  
  
  if(is.data.frame(tpm.matrix)){
    ctr.samples <- rownames(design[design$Group == ctr,])
    tr.samples <- rownames(design[design$Group == tr,])
    deseq2_result$TPM_ctr_mean <- round(apply(deseq2_result[,ctr.samples],1,mean), digits = 3)
    deseq2_result$TPM_tr_mean <- round(apply(deseq2_result[,tr.samples],1,mean), digits = 3)
  }
  deseq2_result$Entrez <- TransGeneID(deseq2_result$Symbol, "Symbol", "Entrez", organism = "hsa")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  result <- list()
  result[["inter"]] <- inter
  result[["result"]] <- deseq2_result
  result[["vsd"]] <- vsd
  result[["design"]] <- design
  return(result)
}


###02.Function of extracting up and donw genes form DESeq2 DEG result

DESeq2_DEG_extract <- function(result.df,design.df,plot, adjust, fc = 1.5, pv = 0.05){
  deseq2_result <- result.df$result
  if(adjust){
    result.df[["up"]] <- deseq2_result |> filter(!is.na(padj)) |> filter(log2FoldChange > log2(fc) & padj < pv)
    result.df[["dw"]] <- deseq2_result |> filter(!is.na(padj)) |> filter(log2FoldChange < -log2(fc) & padj < pv)
  }else{
    result.df[["up"]] <- deseq2_result |> filter(!is.na(pvalue)) |> filter(log2FoldChange > log2(fc) & pvalue < pv)
    result.df[["dw"]] <- deseq2_result |> filter(!is.na(pvalue)) |> filter(log2FoldChange < -log2(fc) & pvalue < pv)  
  }
  if(plot){
    #col_anno=data.frame(row.names = rownames(design.df),Group=design.df$Group,WT1=design.df$WT1_tpm)
    col_anno <- data.frame(row.names = design.df$SampleID,Group=design.df$Group)
    col_anno$Group <- as.factor(col_anno$Group)
    #col_anno_color <- list(Group = c(CD34="#F15A24",APL="#0071BC"))
    #col_anno_color <- list(Group = c(control="#F15A24",HMGA2_sh2="#0071BC"))
    col_anno_color <- list(Group = c("#A81E2C", "#08537C"))
    names(col_anno_color$Group) <- levels(col_anno$Group)
    
    EX_data <- assay(result.df$vsd[c(result.df$up$Row.names,result.df$dw$Row.names)])
    sd_rows <- apply(EX_data, 1, sd)
    EX_data <- EX_data[sd_rows > 0, ]
    EX_data <- t(scale(t(EX_data)))
    EX_data[EX_data > 2] =  2
    EX_data[EX_data < -2] =  -2
    
    #EX_data <- assay(result.df$vsd[c(result.df$up$Row.names,result.df$dw$Row.names)])
    #EX_data <- log2(EX_data+1)
    #EX_data.mean <- matrix(rep(apply(EX_data,1,mean),ncol(EX_data)),ncol=ncol(EX_data))
    #EX_data.sd <- matrix(rep(apply(EX_data,1,sd),ncol(EX_data)),ncol=ncol(EX_data))
    #EX_data <- (EX_data - EX_data.mean) / EX_data.sd
    #EX_data[EX_data > 2] =  2
    #EX_data[EX_data < -2] =  -2
    result.df[["plot"]] <- pheatmap::pheatmap(EX_data,scale="none",
                                              show_colnames = T,
                                              show_rownames = F,
                                              annotation_col = col_anno,
                                              annotation_colors = col_anno_color,
                                              cluster_cols  = T,
                                              cluster_rows = T,
                                              clustering_method = "complete", 
                                              color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")))(100)
    )
  }
  return(result.df)
}


###03.Function of risk groping based on multivariable cox model

RiskScore_multivar_cox <- function(sample.list, gene.list, feature.list,tpm.matrix, clinic.df, cox.df = NULL) {
  #  sample.list <- aml_all_wt1.clinc.rm[aml_all_wt1.clinc.rm$Source.x == "BEAT",]$SampleID
  #  gene.list <- gene_sig.list
  #  tpm.matrix <- aml_all_tpm.matrix
  #  clinic.df <- aml_all_wt1.clinc
  #  feature.list <- c("WT1","TET2","IDH1","IDH2","ageAtDiagnosis")
  
  clinic.df$ageAtDiagnosis <- as.numeric(clinic.df$ageAtDiagnosis)
  sample.list <- gsub("-","\\.", sample.list)  
  col_index <- match(sample.list,colnames(tpm.matrix))
  row_index <- match(gene.list, tpm.matrix$Symbol)
  risk.tpm <- tpm.matrix[row_index,col_index,]
  rownames(risk.tpm) <- tpm.matrix$Symbol[row_index]
  risk.tpm <- t(risk.tpm)
  rownames(risk.tpm) <- gsub("\\.","-",rownames(risk.tpm))
  risk.cox <- merge(risk.tpm, clinic.df, by.x=0,by.y =2, all=FALSE)
  risk.cox[risk.cox$vitalStatus  == "Death",]$vitalStatus <- 1
  risk.cox[risk.cox$vitalStatus  == "Alive",]$vitalStatus <- 0
  risk.cox$vitalStatus <- as.numeric(risk.cox$vitalStatus)
  risk.cox$overallSurvival <- as.numeric(risk.cox$overallSurvival)
  risk.cox <- risk.cox %>% filter(overallSurvival > 0)
  #  survpre_mut <- data.table::as.data.table(risk.cox[,c("WT1","TET2","IDH1","IDH2","ageAtDiagnosis")])
  #  colnames(survpre_mut) <- c("WT1_mut","TET2_mut", "IDH1_mut", "IDH2_mut","Age")
  #  survpre_mut$Age <- as.numeric(survpre_mut$Age)
  survpre_mut <- data.table::as.data.table(risk.cox[,feature.list,drop = FALSE])
  survpre_sig <- risk.cox[,gene.list,drop = FALSE]
  survpre_os <-  risk.cox[,c("overallSurvival","vitalStatus")]
  colnames(survpre_os) <- c("Time","Status")
  surv_model.df <- cbind(survpre_sig,survpre_mut)
  surv_model.df <- as.data.frame(row.names = risk.cox$Row.names, cbind(surv_model.df,survpre_os))
  surv_model.df <- na.omit(surv_model.df)
  survobj <- Surv(surv_model.df$Time,surv_model.df$Status)
  survpre <- as.matrix(surv_model.df[,c(gene.list,feature.list),drop = FALSE])
  model <- coxph( survobj ~ survpre)
  #  if(is.null(cox.df)){
  cox.df.sum <-summary(model)
  cox.df <- as.data.frame(cbind(cox.df.sum$coefficients, cox.df.sum$conf.int))
  cox.df$Feature <- gsub("survpre","",rownames(cox.df))
  cox.df <- data.table::as.data.table(cox.df[,c(10,1,5,6,8,9)])
  colnames(cox.df) <- c("Feature","Coefficient","Pvalue","Hazard.ratio","HR.Lower.95","HR.Upper.95")
  cox.df <- na.omit(cox.df)
  #  }
  ls_df <- lapply(1:nrow(surv_model.df), function(x){
    sample <- rownames(surv_model.df)[x]
    score <-sum(sapply(1:nrow(cox.df),function(i){
      coef = cox.df$Coefficient[i]
      tpm = surv_model.df[x,colnames(surv_model.df) == cox.df$Feature[i]]
      coef * tpm
    }))
    data.frame(Sample = sample, Risk_score=score)
  })
  rs.df <- do.call(rbind, ls_df)
  rs.df$RS_norm<-sapply(1:nrow(rs.df),function(i){
    (rs.df$Risk_score[i] - mean(rs.df$Risk_score))/sd(rs.df$Risk_score)
  })
  rs.df$Risk_Group <- "High_risk"
  rs.df[rs.df$RS_norm < median(rs.df$RS_norm),]$Risk_Group <- "Low_risk"
  rs.df$Sample <- gsub("\\.","-",rs.df$Sample)
  rownames(surv_model.df) <- gsub("\\.","-",rownames(surv_model.df))
  rs.df.clinc <- merge(rs.df,surv_model.df, by.x = 1, by.y =0, all = FALSE)
  rs.df.clinc <- merge(rs.df.clinc, aml_all_wt1.clinc, by.x = 1, by.y =2, all = FALSE)
  result.ls <- list()
  result.ls[["HR"]] <- cox.df
  result.ls[["Clinc"]] <- rs.df.clinc
  return(result.ls)
}


###03.Function of survial of risk groping based on univ cox and TPM

Risk_grouping_cox<-function(gene.list,design.df){
  #gene.list<-risk_gene.list
  #design.df<-Risk_model.beat.df
  
  ls_df <- lapply(match(gene.list,colnames(design.df)), function(i){
    pretty_aml <- data.frame(time   = design.df$Time, 
                             status = design.df$Status, 
                             TPM    = design.df[,i])
    Gene_label <- colnames(design.df)[i]
    model <- coxph( Surv(time, status) ~ TPM, data = pretty_aml)
    cox_result <- summary(model)
    data.frame(row.names = Gene_label, cox_result$conf.int,PV=cox_result$coefficients[5])
  })
  hr.df <- do.call(rbind, ls_df)
  ls_df <- lapply(1:nrow(design.df), function(x){
    sample <- design.df$SampleID[x]
    score <-sum(sapply(1:length(gene.list),function(i){
      coef = log(hr.df[rownames(hr.df) == gene.list[i],]$exp.coef.) 
      tpm = design.df[x,colnames(design.df)==gene.list[i]]
      coef * tpm
    }))
    data.frame(Sample = sample, Risk_score=score)
  })
  rs.df <- do.call(rbind, ls_df)
  rs.df$RS_norm<-sapply(1:nrow(rs.df),function(i){
    (rs.df$Risk_score[i] - mean(rs.df$Risk_score))/sd(rs.df$Risk_score)
  })
  rs.df$Risk_Group <- "High_risk"
  rs.df[rs.df$RS_norm < median(rs.df$RS_norm),]$Risk_Group <- "Low_risk"
  rs.df$Sample <- gsub("\\.","-",rs.df$Sample)
  design.df$SampleID <- gsub("\\.","-",design.df$SampleID)
  rs.df.clinc <- merge(rs.df,design.df, by.x = 1, by.y =1, all = FALSE)
  return(rs.df.clinc)
}



###04.Function of getting risk model dataframe
Risk_model_df<-function(gene.list,design.df,tpm.df){
  #gene.list <- risk_gene.list
  #design.df <- aml_design_TCGA
  #tpm.df <- aml_combind_tpm.df
  index_col<-match(design.df$SampleID,colnames(tpm.df))
  index_row <- match(gene.list,tpm.df$Symbol)
  Risk.df <- tpm.df[index_row,index_col]
  rownames(Risk.df) <- gene.list
  Risk.df <- t(Risk.df)
  Risk.df <- merge(Risk.df,design.df,by.x=0,by.y=1)
  Risk.df<-Risk.df %>% filter(!(is.na(Status) | Status == "Unknown" | Time==0))
  colnames(Risk.df)[1] <- "SampleID"
  Risk.df$Time <- as.numeric(as.character(Risk.df$Time))
  Risk.df$Status <- as.character(Risk.df$Status)
  Risk.df[Risk.df$Status == "Death",]$Status <- 1
  Risk.df[Risk.df$Status == "Alive",]$Status <- 0
  Risk.df$Status <- as.numeric(Risk.df$Status)
  return(Risk.df)
}



###05.Function of survival-related plot basd on risk group

plot_survial_risk<-function(clinic_rs.df){
  #clinic_rs.df<-rs.df.clinc
  blank <- theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank())
  clinic_rs.df$rank <- rank(clinic_rs.df$Risk_score)
  x1<-ggplot(clinic_rs.df,aes(x=rank,y=Risk_score)) + 
    geom_point(aes(color=Risk_Group)) + 
    blank + 
    theme(legend.justification = c("right", "top"), legend.position = c(1,1),legend.title = element_blank()) + 
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(angle=90,hjust=0.5))  +
    #    scale_x_continuous(breaks=c(dim(clinic_rs.df)[1]*0.25, dim(clinic_rs.df)[1]*0.75),labels=c("Low_risk","High_risk")) +
    geom_vline(xintercept = dim(clinic_rs.df)[1]*0.5,linetype=5) +
    scale_color_manual(values=c(c("#FF7F00","#1F78B4"))) +
    ylab("Risk score") + 
    theme(panel.border = element_rect(linetype = 1, size = 1, fill = NA), axis.line = element_blank(),axis.text.x = element_blank()) + 
    theme(legend.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))
  
  x2<-ggplot(clinic_rs.df,aes(x=rank,y=Time)) + 
    geom_point(aes(color=vitalStatus,shape=vitalStatus)) + 
    blank + 
    theme(legend.justification = c("right", "top"), legend.position = c(1,1),legend.title = element_blank()) + 
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(angle=90,hjust=0.5)) +
    scale_x_continuous(breaks=c(dim(clinic_rs.df)[1]*0.25, dim(clinic_rs.df)[1]*0.75),labels=c("Low_risk","High_risk")) +
    geom_vline(xintercept = dim(clinic_rs.df)[1]*0.5,linetype=5) +
    scale_color_manual(values=c("grey40","black")) +
    scale_shape_manual(values=c(1,2)) +
    ylab("Survial days") + 
    theme(panel.border = element_rect(linetype = 1, size = 1, fill = NA), axis.line = element_blank())+ 
    theme(legend.background = element_rect(fill=NA),legend.key = element_rect(fill=NA))
  
  return(plot_grid(x1,x2,ncol=1))
}


###06.Function of Riks model ROC plot 

plot_roc_curve<-function(risk.df){
  #  risk.df <- rs.df.clinc
  #  risk.df$lp <- a
  tmp.df.1 <- survivalROC(Stime = risk.df$Time,
                          status = risk.df$Status,
                          marker = risk.df$Risk_score,
                          predict.time = 365*1,
                          span = 0.1*nrow(risk.df)^(-0.30)
                          #                 method="KM"
  )
  tmp.df.3 <- survivalROC(Stime = risk.df$Time,
                          status = risk.df$Status,
                          marker = risk.df$Risk_score,
                          predict.time = 365*3,
                          span = 0.1*nrow(risk.df)^(-0.30)
                          #                 method="KM"
  )
  tmp.df.5 <- survivalROC(Stime = risk.df$Time,
                          status = risk.df$Status,
                          marker = risk.df$Risk_score,
                          predict.time = 365*5,
                          span = 0.1*nrow(risk.df)^(-0.30)
                          #                 method="KM"
  )
  #  str(tmp.df)
  auc1 <- paste0("1 years", "(AUC=", round(tmp.df.1$AUC,4), ")")
  auc3 <- paste0("3 years", "(AUC=", round(tmp.df.3$AUC,4), ")")
  auc5 <- paste0("5 years", "(AUC=", round(tmp.df.5$AUC,4), ")")
  tmp.df.plot.1 <- data.frame(FP=sort(tmp.df.1$FP),TP=sort(tmp.df.1$TP), PT=rep(auc1,length(tmp.df.1$TP)))
  tmp.df.plot.3 <- data.frame(FP=sort(tmp.df.3$FP),TP=sort(tmp.df.3$TP), PT=rep(auc3,length(tmp.df.3$TP)))
  tmp.df.plot.5 <- data.frame(FP=sort(tmp.df.5$FP),TP=sort(tmp.df.5$TP), PT=rep(auc5,length(tmp.df.5$TP)))
  #  tmp.df.plot.1 <- data.frame(FP=tmp.df.1$FP,TP=tmp.df.1$TP, PT=rep(auc1,length(tmp.df.1$TP)))
  #  tmp.df.plot.3 <- data.frame(FP=tmp.df.3$FP,TP=tmp.df.3$TP, PT=rep(auc3,length(tmp.df.3$TP)))
  #  tmp.df.plot.5 <- data.frame(FP=tmp.df.5$FP,TP=tmp.df.5$TP, PT=rep(auc5,length(tmp.df.5$TP)))
  tmp.df.plot <- rbind(tmp.df.plot.1, tmp.df.plot.3)
  tmp.df.plot <- rbind(tmp.df.plot, tmp.df.plot.5)
  
  
  ggplot(tmp.df.plot,aes(x=FP,y=TP)) + 
    geom_line(aes(color=PT),size = 0.5) + blank + 
    geom_abline(slope = 1,linetype=2, size = 0.5) +
    xlab("1-Specificity") + ylab("Sensitivity") + 
    scale_color_manual(values = brewer.pal(9,"Set1")[c(2,4,5)]) + 
    theme(legend.key = element_blank(),legend.title = element_blank(), 
          legend.position = c(.95, .05), legend.justification = c("right", "bottom")) + 
    theme(panel.border = element_rect(linetype = 1, size = 0.8, fill = NA), 
          axis.line = element_blank())
}


###07.Build function of Volcano Plot

target_for_volcano <- function(deseq2_result.df,gene.list,n = 5000, pc = FALSE, pv = 0.05, fc = 1.5){
  #deseq2_result.df <- RELA_MUTvsWT_NON.DEGs$result
  if(pc == TRUE){ deseq2_result.df <- deseq2_result.df |> filter(Gene_Type == "protein_coding") }
  deg.list <- (deseq2_result.df |> filter(abs(log2FoldChange) > log2(fc), padj < pv))$Symbol
  other.list <- (deseq2_result.df |> filter(! Symbol %in% deg.list, ! is.na(log2FoldChange), ! is.na(padj)))$Symbol
  set.seed(123)
  random.list <- sample(other.list,n)
  result.list <- c(deg.list, random.list)
  result.list <- c(result.list,gene.list)
  return(result.list)
}

volcano_plot_Deseq2 <- function(deseq2_result.df,gene.list, n = 5000, pc = FALSE, pv = 0.05, fc = 1.5, max_x = 100, max_y = 100, top = 5, adjust){
  #  deseq2_result.df <- WT1_MUT.deg$result
  #  target_gene.list <- WT1_status.gene$Symbol
  #  gene.list <- risk_gene.list
  #  pv <- 0.05
  #  fc <- 1.5
  target_gene.list <- target_for_volcano(deseq2_result.df, gene.list,n = n, pc = pc, pv = pv, fc = fc)
  gg <- deseq2_result.df[,c("Symbol","log2FoldChange")]
  if(adjust){ gg$padj <- deseq2_result.df$padj }else{ gg$padj <- deseq2_result.df$pvalue }
  gg <- gg[match(target_gene.list,gg$Symbol),]
  gg <- gg |> filter(! is.na(log2FoldChange), !is.na(padj))
  gene.list <- c((gg |> filter(log2FoldChange > 0) |> top_n(n=-top, wt=padj))$Symbol,gene.list)
  gene.list <- c((gg |> filter(log2FoldChange < 0) |> top_n(n=-top, wt=padj))$Symbol,gene.list)
  index <- match(gene.list,gg$Symbol)
  index <- na.omit(index)
  gg$group = "no"
  gg[gg$log2FoldChange > log2(fc) & gg$padj < pv,]$group <- "up"
  gg[gg$log2FoldChange < -log2(fc) & gg$padj < pv,]$group <- "down"
  gg$color <- gg$group
  gg$color[index] <- "black"
  #mycolour = c("grey", "#A81E2C","#08537C", "black")
  mycolour = c("grey", "#810F7C", "#006D2C", "black")
  names(mycolour) = c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
  #gg[gg$group == "no",]$label <- ""
  if(max(gg$log2FoldChange) > max_x){gg[gg$log2FoldChange > max_x,]$log2FoldChange = max_x}
  if(min(gg$log2FoldChange) < -max_x){gg[gg$log2FoldChange < -max_x,]$log2FoldChange = -max_x}
  if(min(gg$padj) < 10^-max_y){gg[gg$padj < 10^-max_y,]$padj = 10^-max_y}
  p <- ggplot(gg, aes(x = log2FoldChange, y = -log10(padj)))
  p <- p + geom_point_rast(aes(color = group), shape = 16, alpha = 0.6, show.legend = FALSE)
  #  p <- p + geom_point(color = ifelse(gg$Symbol %in% gene.list,"black", NA), shape = 1, show.legend = FALSE)
  p <- p + geom_point(color = ifelse(gg$label == "" ,NA,"black"), shape = 1, show.legend = FALSE)
  p = p + scale_color_manual(values = mycolour)
  p = p + scale_fill_manual(values = mycolour)
  p = p + geom_hline(yintercept = -log10(pv), linetype = "dotted")
  p = p + geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted")
  p = p + ggrepel::geom_text_repel(aes(label=label,color=group),show.legend = FALSE,fontface = "bold", size = 2.5,box.padding = unit(0.8, "lines"),point.padding = unit(0.3, "lines"), segment.size = 0.3, max.overlaps = 10000)
  p = p + theme_test()
  #p= p + xlim(-4,4) + ylim(0,100)
  if(adjust){ p = p + ylab("-log10(padj)") }else{ p = p + ylab("-log10(pvalue)") }
  return(p)
}


###08.Build function of enrich result combind  
enrich_combind <- function(gene,pvc=1,qvc=1, universe = NULL){
  wp2gene <- read.gmt.wp("/lustre/home/acct-medkkw/medlyb/wl_proj/WL234_Lib/database/gmt/wikipathways-20190510-gmt-Homo_sapiens.gmt")
#  wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
  wpid2name <- wp2gene %>% dplyr::select(wpid, name)
  ego_bp <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = pvc, qvalueCutoff  = qvc, readable = TRUE, universe = universe)
  ego_cc <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = pvc, qvalueCutoff  = qvc, readable = TRUE, universe = universe)
  ego_mf <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = pvc, qvalueCutoff  = qvc, readable = TRUE, universe = universe)
  ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = pvc,qvalueCutoff  = qvc, universe = universe)
  ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  ekg<- enrichKEGG(gene, organism = "hsa",pvalueCutoff = pvc,qvalueCutoff  = qvc, universe = universe)
  ekg <- setReadable(ekg, org.Hs.eg.db, keyType = "ENTREZID")
  tmp <- list()
  tmp[["ego_bp"]] <- ego_bp
  tmp[["ego_cc"]] <- ego_cc
  tmp[["ego_mf"]] <- ego_mf
  tmp[["ewp"]] <- ewp
  tmp[["ekg"]] <- ekg
  return(tmp)
}

enrich_combind_s <- function(gene, pvc = 1, qvc = 1, universe = NULL, species = "human") {

  library(clusterProfiler)
  library(dplyr)
  
  if (species %in% c("human", "Homo sapiens", "hsa")) {
    
    library(org.Hs.eg.db)
    current_OrgDb <- org.Hs.eg.db
    kegg_org      <- "hsa"
    
  } else if (species %in% c("mouse", "Mus musculus", "mmu")) {
    
    library(org.Mm.eg.db)
    current_OrgDb <- org.Mm.eg.db
    kegg_org      <- "mmu"
    
  } else {
    stop("Species must be either 'human' or 'mouse'")
  }
  
  ego_bp <- enrichGO(gene = gene, OrgDb = current_OrgDb, keyType = 'ENTREZID', ont = "BP", 
                     pAdjustMethod = "BH", pvalueCutoff = pvc, qvalueCutoff = qvc, 
                     readable = TRUE, universe = universe)
  
  ego_cc <- enrichGO(gene = gene, OrgDb = current_OrgDb, keyType = 'ENTREZID', ont = "CC", 
                     pAdjustMethod = "BH", pvalueCutoff = pvc, qvalueCutoff = qvc, 
                     readable = TRUE, universe = universe)
  
  ego_mf <- enrichGO(gene = gene, OrgDb = current_OrgDb, keyType = 'ENTREZID', ont = "MF", 
                     pAdjustMethod = "BH", pvalueCutoff = pvc, qvalueCutoff = qvc, 
                     readable = TRUE, universe = universe)

  ekg <- enrichKEGG(gene, organism = kegg_org, 
                    pvalueCutoff = pvc, qvalueCutoff = qvc, universe = universe)
  if (!is.null(ekg)) {
    ekg <- setReadable(ekg, current_OrgDb, keyType = "ENTREZID")
  }
  
  tmp <- list()
  tmp[["ego_bp"]] <- ego_bp
  tmp[["ego_cc"]] <- ego_cc
  tmp[["ego_mf"]] <- ego_mf
  tmp[["ekg"]]    <- ekg
  
  return(tmp)
}

enrich_combind_s2 <- function(gene,
                             pvc = 1, qvc = 1,
                             universe = NULL,
                             species = "human",
                             reactome = TRUE,
                             hallmark = TRUE,
                             minGSSize = 10,
                             maxGSSize = 500) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(dplyr)
  })
  
  ## -------- species setup --------
  if (species %in% c("human", "Homo sapiens", "hsa")) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    current_OrgDb <- org.Hs.eg.db
    kegg_org      <- "hsa"
    reactome_org  <- "human"
    msig_species  <- "Homo sapiens"
  } else if (species %in% c("mouse", "Mus musculus", "mmu")) {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    current_OrgDb <- org.Mm.eg.db
    kegg_org      <- "mmu"
    reactome_org  <- "mouse"
    msig_species  <- "Mus musculus"
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
  ego_bp <- enrichGO(gene = gene, OrgDb = current_OrgDb, keyType = "ENTREZID", ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = pvc, qvalueCutoff = qvc,
                     readable = TRUE, universe = universe,
                     minGSSize = minGSSize, maxGSSize = maxGSSize)
  
  ego_cc <- enrichGO(gene = gene, OrgDb = current_OrgDb, keyType = "ENTREZID", ont = "CC",
                     pAdjustMethod = "BH", pvalueCutoff = pvc, qvalueCutoff = qvc,
                     readable = TRUE, universe = universe,
                     minGSSize = minGSSize, maxGSSize = maxGSSize)
  
  ego_mf <- enrichGO(gene = gene, OrgDb = current_OrgDb, keyType = "ENTREZID", ont = "MF",
                     pAdjustMethod = "BH", pvalueCutoff = pvc, qvalueCutoff = qvc,
                     readable = TRUE, universe = universe,
                     minGSSize = minGSSize, maxGSSize = maxGSSize)
  
  ## -------- KEGG --------
  ekg <- enrichKEGG(gene = gene, organism = kegg_org,
                    pvalueCutoff = pvc, qvalueCutoff = qvc,
                    universe = universe,
                    minGSSize = minGSSize, maxGSSize = maxGSSize)
  if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0) {
    ekg <- setReadable(ekg, current_OrgDb, keyType = "ENTREZID")
  }
  
  ## -------- Reactome --------
  erct <- NULL
  if (isTRUE(reactome)) {
    suppressPackageStartupMessages(library(ReactomePA))
    erct <- enrichPathway(gene = gene,
                          organism = reactome_org,
                          pvalueCutoff = pvc,
                          pAdjustMethod = "BH",
                          qvalueCutoff = qvc,
                          universe = universe,
                          readable = TRUE,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize)
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
    sym2ent <- bitr(msig_h$gene_symbol,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = current_OrgDb) %>%
      distinct(SYMBOL, ENTREZID)
    
    msig_h2 <- msig_h %>%
      left_join(sym2ent, by = c("gene_symbol" = "SYMBOL")) %>%
      filter(!is.na(ENTREZID)) %>%
      dplyr::select(gs_name, ENTREZID) %>%
      distinct()
    
    # TERM2GENE format for enricher
    term2gene_h <- msig_h2 %>% dplyr::rename(term = gs_name, gene = ENTREZID)
    
    ehall <- enricher(gene = gene,
                      TERM2GENE = term2gene_h,
                      universe = universe,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvc,
                      qvalueCutoff = qvc,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize)
    
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
    ekg    = ekg,
    reactome = erct,
    hallmark = ehall
  )
  return(out)
}

###09. Build function of enricher_plot
enricher_plot <- function(enricher,bp=1:5,cc=1:5,mf=1:5,wp=1:5,kg=1:5,value=6){
  df <- data.frame(Description=enricher$ego_bp@result[bp,2],Pvalue=enricher$ego_bp@result[bp,value],Type=rep("GO_BP",length(bp)))
  df.tmp <- data.frame(Description=enricher$ego_cc@result[cc,2],Pvalue=enricher$ego_cc@result[cc,value],Type=rep("GO_CC",length(cc)))
  df <- rbind(df,df.tmp)
  df.tmp <- data.frame(Description=enricher$ego_mf@result[mf,2],Pvalue=enricher$ego_mf@result[mf,value],Type=rep("GO_MF",length(mf)))
  df <- rbind(df,df.tmp)
  df.tmp <- data.frame(Description=enricher$ewp@result[wp,2],Pvalue=enricher$ewp@result[wp,value],Type=rep("WikiPath",length(wp)))
  df <- rbind(df,df.tmp)
  df.tmp <- data.frame(Description=enricher$ekg@result[kg,2],Pvalue=enricher$ekg@result[kg,value],Type=rep("KEGG",length(kg)))
  df <- rbind(df,df.tmp)
  #ggplot(df,aes(x=Description,y=-log(Pvalue),fill=Type)) + geom_bar(stat = "identity",show.legend = TRUE) + coord_flip()
  p<-ggplot(df,aes(-log10(Pvalue),fct_reorder(Description, -log10(Pvalue)))) +geom_segment(aes(xend=0, yend = Description,color=Type),linetype = 2,show.legend = FALSE) + geom_point(aes(color=Type),size=5,show.legend = FALSE) + scale_color_manual(values = brewer.pal(5,"Dark2")) + facet_grid(Type~.,scales = 'free',space = 'free_y', switch = "x") + blank + ylab("")
  return(p)
}

###10. Build function of paired gene expression correlation
ccor_paired_gene <- function(EXdata,Gene1,Gene2){
  #  EXdata<-aml_all_tpm.matrix
  #  Gene1 <- "TET2"
  #  Gene2 <- "TP53"
  ##########################
  #           Symbol Length           Type    BA2000    BA2003      BA2004
  #TSPAN6     TSPAN6   4535 protein_coding   0.00000  1.016602   0.3353875
  #TNMD         TNMD   1610 protein_coding   0.00000  0.000000   0.0000000
  ########################
  corr.df <- data.frame(Gene1=log2(as.numeric(EXdata[EXdata$Symbol==Gene1,4:ncol(EXdata)])+1),Gene2=log2(as.numeric(EXdata[EXdata$Symbol==Gene2,4:ncol(EXdata)])+1))
  corr.result <- cor.test(~Gene1 + Gene2,corr.df,method="pearson")
  text<-paste("k=",round(lm( Gene1~Gene2, corr.df)$coefficients[2],4),",","r=",round(corr.result$estimate,4), sep = "")
  pv <- corr.result$p.value
  if (pv < 0.05 & pv > 0.01){sig="*"}else if(pv < 0.01 & pv > 0.001){sig="**"}else if(pv < 0.001 & pv > 0.0001){sig="***"}else if(pv < 0.0001){sig="****"}else{sig="na"}
  blank <- theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank())
  ggplot(corr.df,aes(y=Gene2,x=Gene1)) + 
    geom_point(size = 1, alpha = 0.5, color = "grey") + 
    geom_smooth(method = "lm", se = TRUE,formula = y ~ x) + ggtitle(text)  +  
    annotate("text",label=sig, 
             y=max(corr.df$Gene2,na.rm = TRUE),
             x=median(corr.df[corr.df$Gene1 != 0,]$Gene1,na.rm = TRUE), size=7)  + 
    xlab(Gene1) + ylab(Gene2) + blank
}

###11. Build function of H3K27ac enhancer rank plot
enhancer_rank_plot <- function(label, sample, top, rank.df){
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
  p<-ggplot(enhancer_signal, aes(x = Rank/1000, y= Signal)) + 
    geom_point_rast(aes(color=Type), alpha=1,size=0.5)  + 
    xlab("Rank/1000") + ylab("Enhancer signal") + ggtitle(paste0(sample," H3K27ac")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    blank +
    scale_color_manual(values=c("#A81E2C", brewer.pal(n = 9, name ="Set1")[9])) +
    geom_vline(xintercept=nrow(enhancer_signal[enhancer_signal$Type == "SE", ])/1000,lty=4,col="black",lwd=0.8)
  
  #p <- p + ggrepel::geom_label_repel(data = enhancer_signal[enhancer_signal$CLOSEST_GENE %in% label | enhancer_signal$Rank %in% 1:top,], 
                                     #aes(x = Rank/1000, y = Signal, label = CLOSEST_GENE, fill = Type), 
  p <- p + ggrepel::geom_text_repel(data = enhancer_signal[enhancer_signal$CLOSEST_GENE %in% label | enhancer_signal$Rank %in% 1:top,], 
                                     aes(x = Rank/1000, y = Signal, label = CLOSEST_GENE, color = Type), 
                                     fontface = "bold", size = 2.5, 
                                     box.padding = unit(1, "lines"), 
                                     segment.color = brewer.pal(n = 9, name ="Set1")[2], 
                                     point.padding = unit(0.3, "lines"), 
                                     segment.size = 0.3, show.legend = FALSE, 
                                     label.r = 0.5, nudge_x = 0, nudge_y = 0, max.overlaps = 10000)  
  #+ 
    #scale_fill_manual(values=c(brewer.pal(n = 9, name ="Set3")[4], brewer.pal(n = 9, name ="Set3")[9]))
  p <- p + theme(legend.title = element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"))
  p <- p + theme(panel.border = element_rect(fill=NA,linetype = 1,size = 1), axis.line = element_blank())
  se_num = dim(enhancer_signal |> filter(Type == "SE"))[1]
  cutoff = enhancer_signal$Signal[se_num + 1]
  p <- p + annotate("text", x= max(enhancer_signal$Rank)/1000/2, y = max(enhancer_signal$Signal)/2, label = paste0("Cutoff used: ",cutoff, "\nSuper-Enhancers identified: ", se_num), hjust = 0)
  return(p)
}

###12. Build function of sig gene expression by mut stat(WT1 TET2 IDH1 IDH2)
expression_boxplot_byMut <- function(gene, risk.df, nm.df) {
  #  gene <- "CHRNE"
  #  risk.df <- rs.df.clinc.beat
  #  nm.df <- aml_all_tpm.matrix.nm
  index1 <- match(gene,nm.df$Symbol)
  index2 <- match(gene,colnames(risk.df))
  index3 <- match(c("WT1.x","TET2.x", "IDH1.x", "IDH2.x"),colnames(risk.df))
  group<-rowSums(risk.df[,index3])
  group[group > 0] <-  "MUT"
  group[group == 0] <-  "WT"
  gg1 <- data.frame(Group = rep("NM",length(index1)), TPM=as.numeric(nm.df[index1,4:ncol(nm.df)][1,]), Risk = rep("NM",length(index1)))
  gg2 <- data.frame(Group = group, TPM=as.numeric(risk.df[,index2]), Risk = risk.df[,4])
  
  gg <- rbind(gg1,gg2)
  gg$Group <- factor(gg$Group, levels = c("MUT","WT","NM"))
  gg$Risk <- factor(gg$Risk, levels = c("High_risk","Low_risk","NM"))
  ggplot(gg,aes(x=Risk,y=log2(TPM),color=Risk,fill=Risk)) + 
    geom_boxplot(outlier.shape = 21,show.legend = FALSE,alpha = 0.65) + 
    #  geom_jitter(size=0.5,show.legend = FALSE) + 
    theme(panel.border = element_rect(fill=NA,linetype = 1),
          panel.background = element_blank()) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title= element_text(hjust = 0.5)) + 
    scale_color_manual(values = c("#FF7F00","#1F78B4","grey50")) + 
    scale_fill_manual(values = c("#FF7F00","#1F78B4","grey50")) + 
#    stat_compare_means(show.legend = FALSE, comparisons = list(c("MUT","WT"),c("WT","NM"),c("MUT","NM"))) +
    stat_compare_means(show.legend = FALSE, comparisons = list(c("High_risk","Low_risk"),c("Low_risk","NM"),c("High_risk","NM"))) +
    ggtitle(paste0(gene, " Expression")) + xlab("")
}

###13. Build function of pieplot for homer peak anno
peak_pieplot_byhomer <- function(peak_anno.df, color = brewer.pal(12,"Paired")[c(1,2,3,4,5,7,8,9)], prefix = "pieplot"){
  peak_anno.df$Type<-(sapply(peak_anno.df$Annotation,function(x){strsplit(x," \\(")[[1]][1]}))
  df.pie <- peak_anno.df %>% group_by(Type) %>% summarise(Num=n(),Pro=round(n()/nrow(peak_anno.df)*100,1))
  df.pie$Type <- factor(df.pie$Type,levels = c("Intergenic","intron","exon","promoter-TSS", "TTS", "non-coding", "3' UTR", "5' UTR"))
  df.pie$labs <- paste(df.pie$Pro, "%")
  df.pie <- data.table::as.data.table(df.pie)
  df.pie[,`:=`(TypeII, paste0(Type,"(",labs,")"))]
  df.pie$TypeII <-  fct_reorder(df.pie$TypeII, -df.pie$Num)
  #ggdonutchart(df.pie, "Num", label = "labs", fill = "Type",color = "white",lab.pos = "out", lab.font = c("black","blod",5)) + theme(legend.position = "right")  + scale_fill_manual(values = color) + ggtitle(prefix)
  ggpie(df.pie, "Num", label = "labs", fill = "Type",color = "white",lab.pos = "out", repel = TRUE,lab.font = c("black","blod",5)) + theme(legend.position = "right")  + scale_fill_manual(values = color) + ggtitle(prefix)
}

###14. Build function of heatmap for ChIPseq peaks density by homer
peak_heatmap_byhomer <- function(heatmap.df,order.list,color="#08537C",max=10, prefix = "Heammap",filedir="~/Desktop",width = 2, height = 10, gaps = 0){
  heatmap.df.order<-heatmap.df[match(order.list,heatmap.df$Gene),]
  heatmap.df.order <- as.matrix(heatmap.df.order[,2:ncol(heatmap.df.order)])
  heatmap.df.order <- log2(heatmap.df.order+1)
  heatmap.df.order[heatmap.df.order > max] = max
  pheatmap::pheatmap(heatmap.df.order, kmeans_k = NA, 
                     scale = "none",cellwidth = NA, cellheight = NA, 
                     show_rownames=FALSE, show_colnames = FALSE,
                     annotation_names_col=FALSE, annotation_legend=TRUE, cluster_rows = FALSE, 
                     legend = FALSE,
                     gaps_row = gaps,
                     filename = paste0(filedir,"/ChIP.heatmap.",prefix,".png"),
                     width = width,
                     height = height,
                     cluster_cols = FALSE,color = colorRampPalette(c("white",color))(100))
  legend_dummy <- matrix(seq(0, max, length.out = 100), nrow = 1)
  pheatmap::pheatmap(legend_dummy,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     legend = TRUE,
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     annotation_legend = FALSE,
                     scale = "none",
                     color = colorRampPalette(c("white", color))(100),
                     filename = paste0(filedir,"/ChIP.heatmap.",prefix,".pdf"),
                     width = 3, height = 3, main = prefix)
}

###15. mutation track view
track_view_cre_mut<-function(chr, region, mut_pos, tf_name, num){
  #  chr <- "chr11"
  #  region <- c(32387775,32418917,32426766,32458769)
  #  mut_pos <- c(32421395,32421396,32421397,32426110)
  #  tf_name <- c("WT1","WT1-AS")
  #  num <- c(1,1,1,1)
  SNPs <- GRanges(chr, IRanges(mut_pos, width = 1), strand="-")
  SNPs$color <- "#D7191C"
  SNPs$border <- "#D7191C"
  SNPs$feature.height = 0.05
  SNPs$cex <- 0.8
  TF1 <- geneTrack(get(tf_name[1], org.Hs.egSYMBOL2EG), TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]
  TF1$dat2 <- SNPs
  TF2 <- geneTrack(get(tf_name[2], org.Hs.egSYMBOL2EG), TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]
  TF2$dat2 <- SNPs
  gr <- GRanges(chr, IRanges(region[1],region[4]))
  optSty <- optimizeStyle(trackList(TF2,TF1), theme="bw")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=0.8))
  setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=0.8))
  names(trackList) <- rev(tf_name)
  vp <- viewTracks(trackList, gr=gr,viewerStyle=viewerStyle)
  addGuideLine(region, vp=vp)
  #  return(vp)
}

###16.  H3K27ac  area view
cre_h3k27ac_area_plot<-function(cre_name, maxy, h3k27ac.list, color, dir, region){
  #  cre_name <- "CRE_8366_WT1"
  #  maxy <- 45
  #  col <- "red"
  #  dir <- "/lustre/home/acct-medkkw/medlyb/project/15.APL_TF_ChIPseq/08.MYC_4C/02.bdg_unfold/"
  # region <- c()
  result = list()
  tmp.ls <- lapply(h3k27ac.list, function(x){
#    tmp.bdg <- read.table(paste0("/lustre/home/acct-medkkw/medlyb/project/20.APL_AC_ChIPseq/10.CRE_h3k27ac_signal/",cre_name,"/",x,"_H3K27ac.",cre_name,".unfold.bed"), header = FALSE)
    tmp.bdg <- read.table(paste0(dir,cre_name,"/",x,"_signal.",cre_name,".unfold.bed"), header = FALSE)
    tmp.bdg$Group = x
    tmp.bdg
  })
  h3k27ac.tmp.df <- do.call(rbind, tmp.ls)
  colnames(h3k27ac.tmp.df) <- c("Chr","Start","End", "Score", "Group")
  h3k27ac.tmp.df[h3k27ac.tmp.df$Score > maxy,]$Score = maxy
  line.plot <-h3k27ac.tmp.df %>% group_by(Start) %>% summarise(Score = median(Score))
  line.plot$Group <- h3k27ac.list[1]
  p1<-ggplot(h3k27ac.tmp.df, aes(x=Start,y=Score,fill=Group)) + 
    geom_area(alpha =2/16,show.legend = FALSE,position = "identity") + 
    blank  + 
    ylim(0,maxy) + 
    scale_fill_manual(values= rep(color,16)) + 
    geom_line(data = line.plot, aes(x=Start,y=Score), size = 0.3) + 
    geom_vline(xintercept = region, linetype = 2) +
    ggtitle(cre_name)
  result[["plot"]] <- p1
  result[["df"]] <- h3k27ac.tmp.df
  return(result)
}

###17.GSEA one term plot
gsea_plot_custorm <- function(gsea_ob,  select_term, color, xpos = 3000){
  library("ggrastr")
  #gsea_ob <- aml_phenolyzer.gsea.crc
  #select_term <- 1
  #color <- "#08537C"
  nes <- round(gsea_ob@result[select_term,"NES"], digits = 2)
  pv <- formatC(gsea_ob@result[select_term, "p.adjust"], format = "e", digits = 2)
  #pv <- round(gsea_ob@result[select_term,"p.adjust"], digits = 6)
  #pv <- round(gsea_ob@result[select_term,"pvalue"], digits = 6)
  gsdata <- do.call(rbind, lapply(select_term, enrichplot:::gsInfo, object = gsea_ob))
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
  
  
  p_gsea_1 <- ggplot(gsdata_thin, aes(x = x, y=runningScore)) + 
    geom_line(aes(color = Description), size = 1,show.legend = FALSE) + 
    blank + theme(panel.border = element_rect(fill = NA, linetype = 1, size = 1), axis.line = element_blank()) +
    geom_hline(yintercept = 0, linetype = 2) + 
    scale_color_manual(values = color) + 
    #theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) +
    theme(axis.ticks.x  = element_blank(), axis.text.x = element_blank()) + 
    xlab("") + ylab("Running enrichment score") + 
    annotate("text", label = paste0("NES = ",nes), x=xpos,y= ypos ) + 
    annotate("text", label = paste0("p.adjust = ", pv), x=xpos,y= ypos - 0.05) +
    #annotate("text", label = paste0("P value = ", pv), x=xpos,y= ypos -0.05) +
    theme(plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"))

  p_gsea_2 <- ggplot(gsdata, aes_(x = ~x)) + 
    rasterise(geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax),color=color, show.legend = FALSE,alpha =0.6,size = 0.4), dpi = 300) +
    #geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax),color=color, show.legend = FALSE,alpha =0.6,size = 0.4) + 
    blank + xlab(NULL) + 
    theme(axis.ticks  = element_blank(), axis.text = element_blank()) +
    theme(panel.border = element_rect(fill = NA, linetype = 1, size = 1), axis.line = element_blank()) +
    #geom_hline(yintercept = 0) +
    theme(plot.margin = margin(t = -0.5, r = 0.2, b = 0, l = 1.2, unit = "cm"))

  p_gsea_3 <- ggplot(gsdata,aes_(x = ~x, y = ~geneList)) + 
    #rasterise(geom_segment(aes_(xend = ~x, yend = 0), color = color, show.legend = FALSE)) +
    rasterise(geom_area(color = color, fill = color, show.legend = FALSE),dpi = 300) +
    #geom_area(color = color, fill = color, show.legend = FALSE) + 
    #  scale_colour_gradient(low= brewer.pal(9,"Blues")[6], high =  brewer.pal(9,"Blues")[9]) +
    ylab("Ranked list matric") + 
    xlab("Rank in ordered dataset") + 
    # scale_y_continuous(n.breaks = 3)+
    blank + theme(panel.border = element_rect(fill = NA, linetype = 1, size = 1), axis.line = element_blank()) +
    theme(plot.margin = margin(t = -0.05, r = 0.2, b = 0.2, l = 0.2, unit = "cm")) 

  gsea_body <- plot_grid(p_gsea_1, p_gsea_2, p_gsea_3, nrow= 3, rel_heights = c(8,0.8,4), align = "v")
  #plot_grid(p_gsea_1, p_gsea_2, nrow= 2, rel_heights = c(8,1), align = "v")
  
  final_plot <- ggdraw() + draw_label( title_text, x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, fontface = "bold", size = 12) +
    draw_plot(gsea_body, y = 0, height = 0.94)
  
  return(final_plot)
}

###18. chipseq heatmap peak sort
heatmap_sort <- function(heatmap.df) {
  #heatmap.df<-heatmap.df.1
  bin_num <- (ncol(heatmap.df)-1)
  summit <- ceiling(bin_num/2) + 1
  start <- summit - ceiling(bin_num * 0.5) +1 
  end <- summit + ceiling(bin_num * 0.5) -1
  colnames(heatmap.df)[start]
  colnames(heatmap.df)[summit]
  colnames(heatmap.df)[end]
  heatmap.df$Center_density <- rowSums(heatmap.df[,start:end])
  heatmap.df.order<-heatmap.df %>% dplyr::arrange(-Center_density)
}

###19. chipseq heatmap plot (3 diff peaks)
ChIPseq_heatmap_plot <- function(file.input,brew.color="Oranges",max=10,filename="~/tmp.png",width = 2, height = 10) {
  #  file.input <- c("/lustre/home/acct-medkkw/medlyb/wl_proj/APL_H3K27ac_WangLab/09.APL_CD34_process/heatmap/CHH_H3K27ac_peaks.APL_CD34.heatmap.txt","/lustre/home/acct-medkkw/medlyb/wl_proj/APL_H3K27ac_WangLab/09.APL_CD34_process/heatmap/CHH_H3K27ac_peaks.APL_uniq.heatmap.txt","/lustre/home/acct-medkkw/medlyb/wl_proj/APL_H3K27ac_WangLab/09.APL_CD34_process/heatmap/CHH_H3K27ac_peaks.CD34_uniq.heatmap.txt")
  heatmap.df.1 <-read.table(file.input[1], header = TRUE,sep = "\t",quote = "")
  heatmap.df.2 <-read.table(file.input[2], header = TRUE,sep = "\t",quote = "")
  heatmap.df.3 <-read.table(file.input[3], header = TRUE,sep = "\t",quote = "")
  
  heatmap.df.1 <- heatmap_sort(heatmap.df.1)
  heatmap.df.2 <- heatmap_sort(heatmap.df.2)
  heatmap.df.3 <- heatmap_sort(heatmap.df.3)
  #library(pheatmap)
  heatmap.df.order <- rbind(heatmap.df.1, heatmap.df.2, heatmap.df.3)
  bin_num <- (ncol(heatmap.df.order)-2)
  heatmap.df.order <- as.matrix(heatmap.df.order[,1:bin_num+1])
  heatmap.df.order <- log2(heatmap.df.order+1)
  heatmap.df.order[heatmap.df.order > max] = max
  #png("~/tmp.png",width = 1000,height = 5000)
  pheatmap(heatmap.df.order, kmeans_k = NA, 
           scale = "none",cellwidth = NA, cellheight = NA, 
           show_rownames=FALSE, show_colnames = FALSE,
           annotation_names_col=FALSE, annotation_legend=TRUE, cluster_rows = FALSE, 
           legend = FALSE,
           gaps_row = c(nrow(heatmap.df.1), nrow(heatmap.df.1) + nrow(heatmap.df.2)),
           filename = filename,
           width = width,
           height = height,
           cluster_cols = FALSE,color = colorRampPalette(brewer.pal(n = 9, name =brew.color))(100)
  )
  #dev.off()
}

###20.GO plot custom
go_plot_custom <- function(Enricher, select_term, color, Label, pmin){
#  Enricher <- Enricher.crc_tg$ego_bp
#  select_term <- 1:20
#  color <- "Blues"
#  Label <- "GO_BP"
  Enricher$GeneCounts <- as.numeric(as.data.frame(strsplit(Enricher$GeneRatio, "/"),stringsAsFactors=FALSE)[1,])
  plot.df <- Enricher[select_term,c("Description", "p.adjust", "GeneCounts")]
  plot.df$Label  <- Label
  plot.df[plot.df$p.adjust < pmin,]$p.adjust <- pmin
  ggplot(plot.df,aes(y=fct_reorder(Description, -log10(p.adjust)),x=Label)) + 
  geom_point(aes(size = GeneCounts, color = -log10(p.adjust))) + 
  blank + theme(panel.border = element_rect(fill = NA, linetype = 1, size = 1), axis.line = element_blank()) + 
#  theme(legend.position = "bottom") +
  ylab("") + xlab("") + scale_colour_gradient(low= brewer.pal(9,color)[6], high =  brewer.pal(9,color)[9])
}

##21. GO BP treeplot for DESeq2 result
GO_BP_treeplot_DESeq2 <- function( deg.df, label, pv = 0.05, lfc = log2(1.5)){
  #  deg.df <- HMGA2_sh2.DESeq2.result$result
  #  label <- "shHMGA2"
  up.entrez <- (deg.df %>% filter(padj < pv, log2FoldChange > lfc))$Entrez
  dw.entrez <- (deg.df %>% filter(padj < pv, log2FoldChange < -lfc))$Entrez
  up.entrez <- unique(na.omit(up.entrez))
  dw.entrez <- unique(na.omit(dw.entrez))
  deg.entrez <- unique(c(up.entrez, dw.entrez))
  color = "YlGnBu"
  up.enricher <- enrich_combind(up.entrez)
  dw.enricher <- enrich_combind(dw.entrez)
  deg.enricher <- enrich_combind(deg.entrez)
  ##################go_bp
  ego_bp.up2 <- enrichplot::pairwise_termsim(up.enricher$ego_bp)
  ego_bp.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ego_bp)
  ego_bp.deg2 <- enrichplot::pairwise_termsim(deg.enricher$ego_bp)
  
  p_up_bp <- enrichplot::treeplot(ego_bp.up2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Up-regulated Genes in GO_BP (',label,' n= ',length(up.entrez),')'))
  p_dw_bp <- enrichplot::treeplot(ego_bp.dw2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Down-regulated Genes in GO_BP (',label,' n= ',length(dw.entrez),')'))
  p_deg_bp <- enrichplot::treeplot(ego_bp.deg2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Diff. Expr. Genes in GO_BP (',label,' n= ',length(deg.entrez),')'))
  ##################kegg
  ekg.up2 <- enrichplot::pairwise_termsim(up.enricher$ekg)
  ekg.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ekg)
  ekg.deg2 <- enrichplot::pairwise_termsim(deg.enricher$ekg)
  
  p_up_kegg <- enrichplot::treeplot(ekg.up2, nCluster = 5, group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Up-regulated Genes in KEGG pathway (',label,' n= ',length(up.entrez),')'))
  p_dw_kegg <- enrichplot::treeplot(ekg.dw2, nCluster = 5, group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Down-regulated Genes in KEGG pathway (',label,' n= ',length(dw.entrez),')'))
  p_deg_kegg <- enrichplot::treeplot(ekg.deg2, nCluster = 5, group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Diff. Expr. Genes in KEGG pathway (',label,' n= ',length(deg.entrez),')'))
  
  ###################
  
  p_combind_up <- enricher_plot(up.enricher) + ggtitle(label = paste0('Up-regulated Genes Enrichment Analysis(',label,' n= ',length(up.entrez),')'))
  p_combind_dw <- enricher_plot(dw.enricher) + ggtitle(label = paste0('Down-regulated Genes Enrichment Analysis(',label,' n= ',length(dw.entrez),')'))
  p_combind_deg <- enricher_plot(deg.enricher) + ggtitle(label = paste0('Diff. Expr. Genes Enrichment Analysis(',label,' n= ',length(deg.entrez),')'))
  
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

##22. GSEA for DESeq2 result
GSEA_DEGseq <- function(deg.df, species = "human"){
  # 加载必要的包
  library(msigdbr)
  library(clusterProfiler)
  library(dplyr)
  
  # 1. 根据物种参数设置对应的变量
  if (species %in% c("human", "Homo sapiens", "hsa")) {
    
    library(org.Hs.eg.db)      # 确保加载人类数据库
    msig_species <- "Homo sapiens"
    kegg_org     <- "hsa"
    current_OrgDb <- org.Hs.eg.db
    
  } else if (species %in% c("mouse", "Mus musculus", "mmu")) {
    
    library(org.Mm.eg.db)      # 确保加载小鼠数据库
    msig_species <- "Mus musculus"
    kegg_org     <- "mmu"
    current_OrgDb <- org.Mm.eg.db
    
  } else {
    stop("Species must be either 'human' or 'mouse'")
  }
  
  # 2. 准备基因列表 (这一步不分物种，假设输入都是Entrez ID)
  # 确保 log2FoldChange 是数值型，Entrez 是字符型
  gene.fc <- deg.df[!is.na(deg.df$Entrez) & !duplicated(deg.df$Entrez), ]$log2FoldChange
  names(gene.fc) <- deg.df[!is.na(deg.df$Entrez) & !duplicated(deg.df$Entrez), ]$Entrez
  gene.fc <- sort(gene.fc, decreasing = TRUE)
  
  # 3. 获取 Hallmark 基因集 (使用动态变量 msig_species)
  m_t2g.H <- msigdbr(species = msig_species, category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  # 4. GSEA Hallmark 分析
  gsea.H <- GSEA(gene.fc, TERM2GENE = m_t2g.H, pvalueCutoff = 1)
  # 使用动态变量 current_OrgDb 进行 ID 转换
  if(!is.null(gsea.H)) {
      gsea.H <- setReadable(gsea.H, OrgDb = current_OrgDb, keyType = "ENTREZID")
  }

  # 5. GSEA KEGG 分析 (使用动态变量 kegg_org)
  gsea.kegg <- gseKEGG(gene.fc, organism = kegg_org, keyType = "kegg", pvalueCutoff = 1)
  if(!is.null(gsea.kegg)) {
      gsea.kegg <- setReadable(gsea.kegg, OrgDb = current_OrgDb, keyType = "ENTREZID")
  }

  # 6. GSEA GO_BP 分析 (使用动态变量 current_OrgDb)
  gsea.go_bp <- gseGO(gene.fc, OrgDb = current_OrgDb, ont = "BP", pvalueCutoff = 1, keyType = "ENTREZID")
  if(!is.null(gsea.go_bp)) {
      gsea.go_bp <- setReadable(gsea.go_bp, OrgDb = current_OrgDb, keyType = "ENTREZID")
  }
  
  # 7. 返回结果
  result <- list()
  result[["Hallmark"]] <- gsea.H
  result[["KEGG"]]     <- gsea.kegg
  result[["GO_BP"]]    <- gsea.go_bp
  
  return(result)
}

###24. scRNA  SCT normalization
scRNA_SCT_norm <- function(object_list){
  mut.list <- lapply(X = object_list, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = mut.list, nfeatures = 3000)
  mut.list <- PrepSCTIntegration(object.list = mut.list, anchor.features = features)
  mut.anchors <- FindIntegrationAnchors(object.list = mut.list, normalization.method = "SCT",anchor.features = features)
  combined.sct <- IntegrateData(anchorset = mut.anchors, normalization.method = "SCT")
  DefaultAssay(combined.sct) <- "integrated"
  combined.sct <- ScaleData(combined.sct, verbose = FALSE)
  combined.sct <- RunPCA(combined.sct, npcs = 50, verbose = FALSE)
  combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:40)
  combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:40)
#  combined.sct <- FindClusters(combined.sct, resolution = c(0.2,0.5,0.8,1,1.5,2))
  combined.sct <- FindClusters(combined.sct, resolution = c(0.5,1))
  return(combined.sct)
}

###25. sc-type annotate of immune system
scTYPE_annotation <- function(scRNA_object, assay, tissue = "Immune system"){
  #  scRNA_object <- pbmc.combined.sct
  #  assay <- "integrated"
  # tissue = "Immune system"
  
  ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f",
            "#e4b680","#7ba39d","#b15928","#ffff99","#6a3d9a","#cab2d6","#ff7f00",
            "#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
  lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
  source("/lustre/home/acct-medty/medty-c/script/scTYPE/gene_sets_prepare.R")
  source("/lustre/home/acct-medty/medty-c/script/scTYPE/sctype_score_.R")
  db_ = "/lustre/home/acct-medty/medty-c/script/scTYPE/ScTypeDB_full.xlsx"

  gs_list = gene_sets_prepare(db_, tissue)
  if(assay == "integrated"){
    ex.max = sctype_score(scRNAseqData = scRNA_object@assays$integrated@scale.data, 
                          scale = TRUE, 
                          gs = gs_list$gs_positive, 
                          gs2= gs_list$gs_negative)
  }else if(assay == "RNA"){
    ex.max = sctype_score(scRNAseqData = scRNA_object@assays$RNA@scale.data, 
                          scale = TRUE, 
                          gs = gs_list$gs_positive,
                          gs2= gs_list$gs_negative)
  }else if(assay == "SCT"){
    ex.max = sctype_score(scRNAseqData = scRNA_object@assays$SCT@scale.data, 
                          scale = TRUE, 
                          gs = gs_list$gs_positive,
                          gs2= gs_list$gs_negative)
  }
  CL_results = do.call("rbind", lapply(unique(scRNA_object@meta.data$seurat_clusters),function(cl){
    es.max.cl = sort(
      rowSums(ex.max[,rownames(scRNA_object@meta.data[scRNA_object@meta.data$seurat_clusters == cl,])]), 
      decreasing = !0)
    head(data.frame(cluster = cl, 
                    type = names(es.max.cl), 
                    scores = es.max.cl, 
                    ncells = sum(scRNA_object@meta.data$seurat_clusters==cl)),10)
  }))
  sctype_scores = CL_results %>% group_by(cluster) %>% top_n(n=1, wt = scores)
  sctype_scores <- sctype_scores[! duplicated(sctype_scores$cluster),]
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  scRNA_object@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]
    scRNA_object@meta.data$customclassif[scRNA_object@meta.data$seurat_clusters == j] =
      paste0(as.character(cl_type$type[1]))
  }
  
  # lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
  # 
  # #prepare edges
  # cL_resutls <- CL_results
  # cL_resutls=cL_resutls[order(cL_resutls$cluster),]
  # edges = cL_resutls
  # edges$type = paste0(edges$type,"_",edges$cluster)
  # edges$cluster = paste0("cluster ", edges$cluster)
  # edges = edges[,c("cluster", "type")]
  # colnames(edges) = c("from", "to")
  # rownames(edges) <- NULL
  # 
  # #prepare nodes
  # nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]
  # nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster)
  # nodes_lvl1$Colour = "#f1f1ef"
  # nodes_lvl1$ord = 1
  # nodes_lvl1$realname = nodes_lvl1$cluster
  # nodes_lvl1 = as.data.frame(nodes_lvl1)
  # nodes_lvl2 = c()
  # for (i in 1:length(unique(cL_resutls$cluster))){
  #   dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]
  #   nodes_lvl2 = rbind(nodes_lvl2, 
  #                      data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), 
  #                                 ncells = dt_tmp$scores, 
  #                                 Colour = ccolss[i], ord = 2, 
  #                                 realname = dt_tmp$type))
  # }
  # nodes = rbind(nodes_lvl1, nodes_lvl2)
  # nodes$ncells[nodes$ncells<1] = 1
  # files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]
  # files_db = unique(files_db)
  # nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
  # nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]
  # nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  # mygraph <- graph_from_data_frame(edges, vertices=nodes)
  # ##plot
  # gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  #   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) +
  #   geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  #   theme_void() + 
  #   geom_node_text(aes(filter=ord==2, 
  #                      label=shortName, colour=I("#ffffff"), 
  #                      fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1))) + 
  #   geom_node_label(aes(filter=ord==1,  label=shortName, 
  #                       colour=I("#000000"), size = I(3), fill="white", parse = T), 
  #                   repel = !0, segment.linetype="dotted")
  result <- list()
  result[["scRNA_object"]] <- scRNA_object 
#  result[["gggr_plot"]] <- gggr
  result[["ccolss"]] <- ccolss
  result[["sctype_scores"]] <- sctype_scores
  result[["CL_results"]] <- CL_results
  
  return(result)
}

###26. volcano plot DEGs for scRNA
volcano_plot_scRNA<- function(findmarkers.df,  gene.list, pv = 0.05, fc = 1.5, top =5, max_x = 1, max_y=50){
  #  gene.list <- c("BIRC3","CCR7", "NFKBIA", "TNFAIP3","REL","BCL3","BCL2","HSPA1A","SOD1","HSPB1", "PPP1R15A", "DNAJA1")
  #  fc <- 1.09
  #  pv <- 0.05
  #  max_x = 1
  #  max_y = 50
  #  gg <- CD4_Naive_T.MTvsWT.deg
  #  gg$Symbol <- as.character(rownames(CD4_Naive_T.MTvsWT.deg))
  gg <- findmarkers.df
  gg$Symbol <- as.character(rownames(findmarkers.df))
#  gg %>% slice_head(n=top)
  gene.list <- c((gg |> filter(avg_log2FC > fc) |> top_n(n=-top, wt=p_val))$Symbol,gene.list)
  gene.list <- c((gg |> filter(avg_log2FC < fc) |> top_n(n=-top, wt=p_val))$Symbol,gene.list)
  #gene.list <- c((gg |> filter(avg_log2FC > 0) |> slice_head(n=top))$Symbol,gene.list)
  #gene.list <- c((gg |> filter(avg_log2FC < 0) |> slice_head(n=top))$Symbol,gene.list)
  #gene.list <- c(rownames(gg %>% slice_head(n=8)),gene.list)
  index <- match(gene.list,gg$Symbol)
  index <- na.omit(index)
  gg$group = "no"
  try(gg[gg$avg_log2FC > log2(fc) & gg$p_val_adj < pv,]$group <- "up")
  try(gg[gg$avg_log2FC < -log2(fc) & gg$p_val_adj < pv,]$group <- "down")
  #try(gg[gg$avg_log2FC > log2(fc) & gg$p_val < pv,]$group <- "up")
  #try(gg[gg$avg_log2FC < -log2(fc) & gg$p_val < pv,]$group <- "down")
  gg$color <- gg$group
  gg$color[index] <- "black"
  #mycolour = c("grey", "#B30000", "#08519C", "black")
  #mycolour = c("grey", "#810F7C", "#006D2C", "black")
  mycolour = c("grey", "#A81E2C","#08537C", "black")
  names(mycolour) = c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
#  gg[gg$group == "no",]$label <- ""
#  if(max(gg$avg_log2FC) > max_x){gg[gg$avg_log2FC > max_x,]$avg_log2FC = max_x}
#  if(min(gg$avg_log2FC) < -max_x){gg[gg$avg_log2FC < -max_x,]$avg_log2FC = -max_x}
#  if(min(gg$p_val_adj) < 10^-max_y){gg[gg$p_val_adj < 10^-max_y,]$p_val_adj = 10^-max_y}
  p <- ggplot(gg, aes(x = avg_log2FC, y = -log10(p_val_adj)))
  #p <- ggplot(gg, aes(x = avg_log2FC, y = -log10(p_val)))
  p <- p + geom_point_rast(aes(fill=group, size = pct.1), shape = 21, alpha = 0.6, show.legend = TRUE)
  p <- p + geom_point_rast(aes(colour = color, size = pct.1), shape = 21, alpha = 0.6, show.legend = TRUE)
  p = p + scale_color_manual(values = mycolour)
  p = p + scale_fill_manual(values = mycolour)
  p = p + geom_hline(yintercept = -log10(pv), linetype = "dotted")
  p = p + geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted")
  p = p + ggrepel::geom_text_repel(aes(label=label,color=group),show.legend = FALSE,fontface = "bold", size = 2.5,box.padding = unit(0.8, "lines"),point.padding = unit(0.3, "lines"), segment.size = 0.3, max.overlaps=10000)
  p = p + blank
  p = p + theme(panel.border = element_rect(fill = NA, linetype = 1,size = 1 ), axis.line = element_blank())
#  p = p + xlim(-max_x,max_x) + ylim(0,max_y)
  return(p)
}


#27. GO BP treeplot for scRNAseq result
GO_BP_treeplot_scRNAseq <- function( deg.df, label, pv = 0.05, lfc = 0.25){
  #  deg.df <- XGJ.NEMOBvsNEMOA.degs
  #  label <- "NEMOBvsNEMOA"
  deg.df$Entrez <- TransGeneID(rownames(deg.df), "Symbol", "Entrez", organism = "hsa")
  up.entrez <- (deg.df %>% filter(p_val < pv, avg_log2FC > lfc))$Entrez
  dw.entrez <- (deg.df %>% filter(p_val < pv, avg_log2FC < -lfc))$Entrez
  up.entrez <- unique(na.omit(up.entrez))
  dw.entrez <- unique(na.omit(dw.entrez))
  color = "YlGnBu"
  up.enricher <- enrich_combind(up.entrez)
  dw.enricher <- enrich_combind(dw.entrez)
  ##################go_bp
  ego_bp.up2 <- enrichplot::pairwise_termsim(up.enricher$ego_bp)
  ego_bp.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ego_bp)
  p_up_bp <- enrichplot::treeplot(ego_bp.up2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Up-regulated Genes in GO_BP (',label,' n= ',length(up.entrez),')'))
  p_dw_bp <- enrichplot::treeplot(ego_bp.dw2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Down-regulated Genes in GO_BP (',label,' n= ',length(dw.entrez),')'))
  ##################kegg
  ekg.up2 <- enrichplot::pairwise_termsim(up.enricher$ekg)
  ekg.dw2 <- enrichplot::pairwise_termsim(dw.enricher$ekg)
  p_up_kegg <- enrichplot::treeplot(ekg.up2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Up-regulated Genes in KEGG pathway (',label,' n= ',length(up.entrez),')'))
  p_dw_kegg <- enrichplot::treeplot(ekg.dw2, nCluster = 5,group_color  = brewer.pal(5,"Set1"), hclust_methd = "average", offset = 10) + scale_colour_gradient(high= brewer.pal(9,color)[6], low =  brewer.pal(9,color)[9]) + ggtitle(label = paste0('Down-regulated Genes in KEGG pathway (',label,' n= ',length(dw.entrez),')'))
  ###################
  
  p_combind_up <- enricher_plot(up.enricher) + ggtitle(label = paste0('Up-regulated Genes Enrichment Analysis(',label,' n= ',length(up.entrez),')'))
  p_combind_dw <- enricher_plot(dw.enricher) + ggtitle(label = paste0('Down-regulated Genes Enrichment Analysis(',label,' n= ',length(dw.entrez),')'))
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

##28. GSEA for scRNAseq result
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
  
#  gsea.kegg <- gseKEGG(
#    geneList = gene.fc,
#    organism = kegg_org,
#    keyType = "kegg",
#    pvalueCutoff = 1
#  )
#  gsea.kegg <- setReadable(gsea.kegg, OrgDb = org_db, keyType = "ENTREZID")
  
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

##29. Correlation analysis
ccor_paired_obj <- function(EXdata,obj1,obj2, Group){
  corr.df <- data.frame(obj1=EXdata[,obj1],obj2=EXdata[,obj2], Group = EXdata[,Group])
  corr.result <- cor.test(~obj1 + obj2,corr.df,method="pearson")
  text<-paste("k=",round(lm( obj1~obj2, corr.df)$coefficients[2],4),",","r=",round(corr.result$estimate,4), sep = "")
  pv <- corr.result$p.value
  if (pv < 0.05 & pv > 0.01){sig="*"}else if(pv < 0.01 & pv > 0.001){sig="**"}else if(pv < 0.001 & pv > 0.0001){sig="***"}else if(pv < 0.0001){sig="****"}else{sig="na"}
  blank <- theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank())
  ggplot(corr.df,aes(x=obj2,y=obj1)) + 
    geom_point(size = 3, aes(color = Group)) + 
    geom_smooth(method = "lm", se = TRUE,formula = y ~ x) + ggtitle(text)  +  
    annotate("text",label=sig, 
             y=max(corr.df$obj1,na.rm = TRUE),
             x=min(corr.df$obj2) + (max(corr.df$obj2)-min(corr.df$obj2))/2, size=7)  + 
    theme_test() + xlab(obj2) + ylab(obj1) +
    scale_color_manual(values = c("#08537C","#A81E2C"))
}

##30. trackview peak region of interest
trackview_peak_roi <-function(chrom, tf_name, ylim = 0,  extend = 1500, info, type = "BigWig", color = c("#08537C", "#A63603", "#A81E2C", "#005824", "#08537C", "#7A0177", "#023858")){
  
  #tf_name <- "ZMIZ1"
  #info <ID BW BAM>
  #chrom <- "chr21:34,787,801-36,004,667"
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(trackViewer)
  chr <- sub(":.*", "", chrom)
  region_str <- sub(".*:", "", chrom)
  region_split <- strsplit(region_str, "-")[[1]]
  region <- as.numeric(gsub(",", "", region_split))
  
  TF1 <- geneTrack(get(tf_name[1], org.Hs.egSYMBOL2EG), TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]
  #chr_tf <- as.character(unique(seqnames(TF1@dat)))
  #start_tf <- min(start(TF1@dat))
  #end_tf <-  max(start(TF1@dat))
  
  gr <- GRanges(chr, IRanges( min(region) - extend,max(region) + extend ))
  if(type == "BAM"){
    Score <- lapply(1:nrow(info), function(x){
      importBam(info$BAM_t[x], ranges = gr, pairs = FALSE)
    }) 
  } else if (type == "BigWig") {
    Score <- lapply(1:nrow(info), function(x){
      importScore(info$BW[x], format = "BigWig", ranges = gr)
    })
  } else if (type == "BED") {
    Score <- lapply(1:nrow(info), function(x){
      importScore(info$BED[x], format = "BED", ranges = gr)
    })
  } else {
    stop("Unsupported type: ", type, ". Expected 'BAM' or 'BigWig'.")
  }
  optSty <- optimizeStyle(do.call(trackViewer::trackList, c(list(TF1), Score)), theme="bw")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  
  Color <- colorRampPalette(color)(nrow(info))
  Color <- c("black",Color)
  
  if(length(ylim) >1){
    for( i in 2: length(trackList)){
      setTrackStyleParam(trackList[[i]], "ylim", c(0,ylim[i-1]))
    }
  }else if(length(ylim) == 1){
    if(ylim > 0)
      for( i in 2: length(trackList)){
        setTrackStyleParam(trackList[[i]], "ylim", c(0,ylim))
      }
  }
  
  for( i in 1:length(trackList)){
    setTrackStyleParam(trackList[[i]], "color", Color[i])
  }
  for(i in 1:length(trackList)){
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=.8))
  }
  setTrackStyleParam(trackList[[1]], "height", 0.03)
  names(trackList) <- c(tf_name, info$ID)
  vp <- viewTracks(trackList, gr=gr,viewerStyle=viewerStyle)
  addGuideLine(region, vp=vp)
}


##31. Compare pairwise DESeq2 results
Compare_pairwise_Deseq2 <- function(deg1.df, deg2.df, label1 = "DESeq 1",label2 = "DESeq 2", gene.list = c("ZMIZ1","MEF2D"), fc = 1.5, top = 5){
  #deg1.df <- shZMIZ1.DEGs.ls$MOLM13_SH4_b2$result
  #deg2.df <- shZMIZ1.DEGs.ls$MOLM13_4E_b5$result
  #label1 <- "MOLM13_SH4"
  #label2 <- "MOLM13_4E"
  #fc <- 1.5
  #top <- 5
  #gene.list <- c("ZMIZ1","MEF2D")
  deg1.df <- na.omit(deg1.df[,c(1:4,6,9,10,ncol(deg1.df)-1)])
  deg2.df <- na.omit(deg2.df[,c(1,6,9,10,ncol(deg2.df)-1,ncol(deg2.df))])
  deg.merge.df <- merge(deg1.df, deg2.df, by.x = 1, by.y = 1)
  
  corr.result <- cor.test(~log2FoldChange.x + log2FoldChange.y,deg.merge.df,method="pearson")
  text<-paste("k=",round(lm( log2FoldChange.x~log2FoldChange.y, deg.merge.df)$coefficients[2],4),",","r=",round(corr.result$estimate,4), sep = "")
  
  pv <- corr.result$p.value
  if (pv < 0.05 & pv > 0.01){sig="*"}else if(pv < 0.01 & pv > 0.001){sig="**"}else if(pv < 0.001 & pv > 0.0001){sig="***"}else if(pv < 0.0001){sig="****"}else{sig="na"}
  
  max <- max(max(deg.merge.df$log2FoldChange.x),deg.merge.df$log2FoldChange.y)
  
  deg.merge.df$group = "no"
  deg.merge.df[deg.merge.df$log2FoldChange.x > log2(fc) & deg.merge.df$log2FoldChange.y > log2(fc),]$group <- "up"
  deg.merge.df[deg.merge.df$log2FoldChange.x < -log2(fc) & deg.merge.df$log2FoldChange.y < -log2(fc),]$group <- "down"
  
  
  gene.list <- c((deg.merge.df |> filter(group == "up") |> top_n(n=top, wt=log2FoldChange.y * log2FoldChange.x))$Symbol,gene.list)
  gene.list <- c((deg.merge.df |> filter(group == "down") |> top_n(n=top, wt=log2FoldChange.y * log2FoldChange.x))$Symbol,gene.list)
  index <- match(gene.list,deg.merge.df$Symbol)
  index <- na.omit(index)
  
  deg.merge.df$color <- deg.merge.df$group
  deg.merge.df$color[index] <- "black"
  mycolour = c("grey", "#A81E2C","#08537C", "black")
  names(mycolour) = c("no", "up", "down", "black")
  
  deg.merge.df$label <- ""
  deg.merge.df$label[index] <- deg.merge.df$Symbol[index]
  deg.merge.df[deg.merge.df$group == "no",]$label <- ""
  
  p <- ggplot(deg.merge.df,aes(y=log2FoldChange.y,x=log2FoldChange.x)) + 
    geom_point_rast(shape = 16, alpha = 0.6, show.legend = FALSE, aes(color = group)) + 
    geom_point(color = ifelse(deg.merge.df$label == "" ,NA,"black"), shape = 1, show.legend = FALSE) + 
    geom_smooth(method = "lm", se = TRUE,formula = y ~ x) + 
    ggtitle(text)  +  ylim(-max,max) + xlim(-max,max) + 
    geom_hline(yintercept = c(-log2(fc), log2(fc)), linetype = "dotted") +
    geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted") +
    ggrepel::geom_text_repel(aes(label=label,color=group),show.legend = FALSE,
                             fontface = "bold", size = 2.5,box.padding = unit(0.8, "lines"),
                             point.padding = unit(0.3, "lines"), segment.size = 0.3, max.overlaps = 10000) +
    annotate("text",label=sig, 
             y=max(deg.merge.df$log2FoldChange.y,na.rm = TRUE),
             x=median(deg.merge.df[deg.merge.df$log2FoldChange.x != 0,]$log2FoldChange.x,na.rm = TRUE), size=7)  + 
    xlab(paste0("log2FoldChange of ",label1)) + ylab(paste0("log2FoldChange of ",label2)) + 
    scale_color_manual(values = mycolour) + scale_fill_manual(values = mycolour) + 
    theme_test()
  return(p)
}

#' 比较两个DEG分析结果并生成散点图
#'
#' @param deg_result1 第一个DEG分析结果数据框
#' @param deg_result2 第二个DEG分析结果数据框
#' @param label_1 X轴标签(默认: "Log2FoldChange in Sample1")
#' @param label_2 Y轴标签(默认: "Log2FoldChange in Sample2")
#' @param fc 倍数变化阈值(默认: 1.5)
#' @param pv P值阈值(默认: 0.05)
#' @param adjust 是否使用校正后的P值(默认: FALSE)
#' @param pc 是否只保留protein_coding基因(默认: TRUE)
#' @param bg_num 背景基因随机抽样数量(默认: 5000)
#' @param limit 坐标轴范围(默认: 9)
#' @param seed 随机种子(默认: 614)
#' @param show_cor 是否显示相关性检验结果(默认: TRUE)
#' @param plot_title 图表标题(默认: NULL,不显示标题)
#' @param goi 感兴趣的基因列表,会被标注在图上(默认: NULL)
#' @param top 自动标注每组(upup/dwdw)中top N个基因(默认: 5)
#' @param label_size 基因标签字体大小(默认: 2.5)
#'
#' @return ggplot对象
#' @export
#'
#' @examples
#' plot_deg_comparison(
#'   deg_result1 = H2A_OV.DEGS.ls$NALM6_b2$DEG$result,
#'   deg_result2 = H2A_OV.DEGS.ls$RCV_b2$DEG$result,
#'   label_1 = "Log2FoldChange in NALM6",
#'   label_2 = "Log2FoldChange in RCV",
#'   goi = c("CD96"),
#'   plot_title = "DEG Comparison: NALM6 vs RCV"
#' )

plot_deg_comparison <- function(deg_result1, 
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
                                label_size = 2.5) {
  
  # 加载必需的包
  require(dplyr)
  require(ggplot2)
  require(ggrastr)
  require(ggrepel)
  
  # 选择P值列名
  pval_col <- if(adjust) "padj" else "pvalue"
  
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
      lfc1 > log2_fc  & lfc2 > log2_fc  & pv1 < pv & pv2 < pv ~ "upup",
      lfc1 < -log2_fc & lfc2 < -log2_fc & pv1 < pv & pv2 < pv ~ "dwdw",
      lfc1 > log2_fc  & lfc2 < -log2_fc & pv1 < pv & pv2 < pv ~ "updw",
      lfc1 < -log2_fc & lfc2 > log2_fc  & pv1 < pv & pv2 < pv ~ "dwup",
      TRUE ~ "other"
    )
  }
  
  # 分类函数(仅考虑FC)
  classify_fc_only <- function(lfc1, lfc2) {
    case_when(
      lfc1 > log2_fc  & lfc2 > log2_fc   ~ "upup",
      lfc1 < -log2_fc & lfc2 < -log2_fc  ~ "dwdw",
      lfc1 > log2_fc  & lfc2 < -log2_fc  ~ "updw",
      lfc1 < -log2_fc & lfc2 > log2_fc   ~ "dwup",
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
  if(pc) {
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
      !is.na(LFC_1), !is.na(PV_1),
      !is.na(LFC_2), !is.na(PV_2)
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
  
  if(show_cor) {
    cat("\n相关性检验结果:\n")
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
  cor_label <- paste0("italic(K)==", round(slope_k, 4), "*','~italic(r)==", round(cor_r, 4), "*' ", sig_stars, "'")
  
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
    geom_point_rast(aes(color = Group2), shape = 16, alpha = 0.6, show.legend = FALSE) + 
    # 第二层:按Group着色(FC+PV标准)的空心点
    #geom_point_rast(aes(color = Group), shape = 1, alpha = 0.6, show.legend = FALSE) + 
    geom_point_rast(data = vs.degs |> filter(Group != "other"),
                    aes(color = Group), shape = 1, alpha = 0.6, show.legend = FALSE) +
    # 第三层:标注感兴趣基因的黑色空心圆
    geom_point(color = ifelse(vs.degs$GOI == "", NA, "black"), shape = 1, show.legend = FALSE) + 
    # 坐标轴范围
    xlim(-limit, limit) + 
    ylim(-limit, limit) + 
    # 主题和标签
    theme_test() + 
    xlab(label_1) + 
    ylab(label_2) + 
    # 基因名标注
    geom_text_repel(
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
    annotate("text", 
             x = -limit * 0.6, 
             y = limit * 0.9, 
             label = cor_label, 
             size = 4, 
             hjust = 0.5,
             parse = TRUE) + 
    # 颜色设置
    scale_color_manual(values = mycolour)
  
  # 添加标题(如果提供)
  if(!is.null(plot_title)) {
    p <- p + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
}

## DEGs manhattan
plot_deg_manhattan <- function(
    deg_list,
    color_map,
    species = c("mouse", "human"),
    deg_cols = c("Symbol", "log2FoldChange", "pvalue", "padj"),
    symbol_col = "Symbol",
    lfc_col = "log2FoldChange",
    p_col = "pvalue",
    padj_col = "padj",
    gene_type_filter = "protein_coding",
    remove_rik = TRUE,
    top_n = 10,
    cap_value = 20,
    chr_keep = NULL,
    facet_nrow = 1
) {
  species <- match.arg(species)
  
  ## required packages
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  require(ggrastr)
  
  ## species-specific genome
  if (species == "mouse") {
    gene_anno_file <- "~/database/GRCm39/gene_len.vM38.txt"
    require(BSgenome.Mmusculus.UCSC.mm39)
    bsgenome_obj <- BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39
    if (is.null(chr_keep)) {
      chr_keep <- c(paste0("chr", 1:19), "chrX", "chrY")
    }
  } else if (species == "human") {
    gene_anno_file <- "~/database/GENCODE/gene_len.v43.new.txt"
    require(BSgenome.Hsapiens.UCSC.hg38)
    bsgenome_obj <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    if (is.null(chr_keep)) {
      chr_keep <- c(paste0("chr", 1:22), "chrX", "chrY")
    }
  }
  
  ## check deg_cols
  if (length(deg_cols) != 4) {
    stop("deg_cols must contain exactly 4 column names: Symbol, log2FoldChange, pvalue, padj")
  }
  
  symbol_col <- deg_cols[1]
  lfc_col    <- deg_cols[2]
  p_col      <- deg_cols[3]
  padj_col   <- deg_cols[4]
  
  ## read gene annotation
  gene_anno <- read.table(gene_anno_file, header = TRUE, sep = "\t", check.names = FALSE)
  
  required_anno_cols <- c("Symbol", "Chr", "Start")
  missing_anno <- setdiff(required_anno_cols, colnames(gene_anno))
  if (length(missing_anno) > 0) {
    stop("gene_anno_file is missing required columns: ", paste(missing_anno, collapse = ", "))
  }
  
  seqlen <- GenomeInfoDb::seqlengths(bsgenome_obj)[chr_keep]
  seqlen <- seqlen[!is.na(seqlen)]
  offset <- c(0, cumsum(as.numeric(seqlen))[-length(seqlen)])
  names(offset) <- names(seqlen)
  
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


## monocle3
run_monocle3_from_seurat_umap <- function(
    seu,
    ref_seu = NULL,
    assay = "SCT",
    layer = "counts",
    reduction = "umap",
    alignment_group = "SampleID",
    seurat_cluster_col = "seurat_clusters",
    root_clusters = NULL,
    num_dim = 50,
    use_partition = FALSE,
    cluster_cells_first = TRUE
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(monocle3)
    library(SingleCellExperiment)
  })
  # 1. prepare counts
  counts <- GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )
  gene_anno <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts)
  )
  cell_meta <- seu@meta.data
  # make sure rownames of metadata match cells
  cell_meta <- cell_meta[colnames(counts), , drop = FALSE]
  cds <- new_cell_data_set(
    counts,
    cell_metadata = cell_meta,
    gene_metadata = gene_anno
  )
  # 2. preprocess and align
  cds <- preprocess_cds(
    cds,
    num_dim = num_dim
  )
  if (!is.null(alignment_group)) {
    if (!alignment_group %in% colnames(colData(cds))) {
      stop("alignment_group not found in cell metadata: ", alignment_group)
    }
    cds <- align_cds(
      cds,
      alignment_group = alignment_group
    )
  }
  # 3. use external Seurat UMAP
  if (is.null(ref_seu)) {
    ref_seu <- seu
  }
  if (!reduction %in% Reductions(ref_seu)) {
    stop("Reduction not found in ref_seu: ", reduction)
  }
  umap_coord <- Embeddings(ref_seu, reduction = reduction)
  missing_cells <- setdiff(colnames(cds), rownames(umap_coord))
  if (length(missing_cells) > 0) {
    stop(
      "Some cells in cds are missing from reference UMAP coordinates. n = ",
      length(missing_cells)
    )
  }
  umap_coord <- umap_coord[colnames(cds), , drop = FALSE]
  reducedDims(cds)$UMAP <- umap_coord
  # 4. cluster and learn graph
  if (cluster_cells_first) {
    cds <- cluster_cells(
      cds,
      reduction_method = "UMAP"
    )
  }
  cds <- learn_graph(
    cds,
    use_partition = use_partition
  )
  # 5. order cells
  if (!is.null(root_clusters)) {
    if (!seurat_cluster_col %in% colnames(colData(cds))) {
      stop("seurat_cluster_col not found in metadata: ", seurat_cluster_col)
    }
    root_cells <- colnames(cds)[
      as.character(colData(cds)[[seurat_cluster_col]]) %in% as.character(root_clusters)
    ]
    if (length(root_cells) == 0) {
      stop("No root cells found. Check root_clusters and seurat_cluster_col.")
    }
    cds <- order_cells(
      cds,
      root_cells = root_cells
    )
  }
  return(cds)
}

## scRNAseq SCT
SCT_METHOD_V3 <- function(
    seu,
    # SCTransform / PCA
    var_features_n = 3000,
    sct_method = "glmGamPoi",
    npcs = 30,
    
    # Integration / graph
    dims = 1:30,
    k.anchor = 20,
    k.weight = 20,
    resolution = 0.5,
    
    # UMAP (影响“分得太开”的视觉表现很大)
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    umap_spread = 1.0,
    
    # housekeeping
    join_layers = TRUE,
    maxSize_GB = 35,
    verbose = FALSE
){
  options(future.globals.maxSize = maxSize_GB * 1024^3)
  
  # ensure SampleID
  seu$SampleID <- factor(seu$SampleID)
  
  # RNA layers: join then split by SampleID
  if (join_layers) {
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$SampleID)
  
  # SCT
  seu <- SCTransform(
    seu,
    variable.features.n = var_features_n,
    method = sct_method,
    assay = "RNA",
    verbose = verbose
  )
  
  DefaultAssay(seu) <- "SCT"
  
  # PCA
  seu <- RunPCA(seu, npcs = max(dims), verbose = verbose, assay = "SCT")
  # 说明：用 max(dims) 保证后面 dims 有足够PC；你也可用 npcs 参数
  # 若你想严格用 npcs，可改成 RunPCA(..., npcs = npcs) 并确保 max(dims) <= npcs
  
  # Integration
  seu <- IntegrateLayers(
    object = seu,
    method = RPCAIntegration,
    normalization.method = "SCT",
    orig.reduction = "pca",
    new.reduction = "integrated.dr",
    verbose = verbose,
    k.anchor = k.anchor,
    k.weight = k.weight,
    dims = dims
  )
  
  # Neighbors / clusters
  seu <- FindNeighbors(seu, reduction = "integrated.dr", dims = dims, graph.name = "integrated_snn")
  seu <- FindClusters(seu, graph.name = "integrated_snn", resolution = resolution)
  
  # UMAP
  seu <- RunUMAP(
    seu,
    reduction = "integrated.dr",
    dims = dims,
    n.neighbors = umap_n_neighbors,
    min.dist = umap_min_dist,
    spread = umap_spread,
    verbose = verbose
  )
  
  return(seu)
}

## CellChat
Build_CellChat_object <- function(
    seu,
    species = c("human", "mouse"),
    group.by = "CellType",
    sample.by = "SampleID",
    cluster.by = NULL,
    assay = "RNA",
    layer = "data",
    min.cells = 10,
    workers = 4,
    maxSize = 40 * 1024^3,
    type = "triMean",
    use_parallel = TRUE,
    verbose = TRUE
) {
  species <- match.arg(species)
  suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(future)
  })
  if (!group.by %in% colnames(seu@meta.data)) {
    stop("group.by column not found in seu@meta.data: ", group.by)
  }
  if (!sample.by %in% colnames(seu@meta.data)) {
    stop("sample.by column not found in seu@meta.data: ", sample.by)
  }
  if (!is.null(cluster.by) && !cluster.by %in% colnames(seu@meta.data)) {
    stop("cluster.by column not found in seu@meta.data: ", cluster.by)
  }
  data.input <- GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )
  meta <- data.frame(
    labels = seu@meta.data[[group.by]],
    samples = seu@meta.data[[sample.by]],
    row.names = colnames(seu)
  )
  if (!is.null(cluster.by)) {
    meta$clusters <- seu@meta.data[[cluster.by]]
  }
  meta$labels <- as.factor(meta$labels)
  meta$samples <- as.factor(meta$samples)
  if (verbose) {
    message("CellChat grouping column: ", group.by)
    message("Sample column: ", sample.by)
    message("Number of cells: ", ncol(data.input))
    message("Number of groups: ", length(unique(meta$labels)))
    print(table(meta$labels))
  }
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels"
  )
  if (species == "human") {
    cellchat@DB <- CellChatDB.human
  } else {
    cellchat@DB <- CellChatDB.mouse
  }
  cellchat <- subsetData(cellchat)
  if (use_parallel) {
    options(future.globals.maxSize = maxSize)
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = type
  )
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min.cells
  )
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )
  future::plan(future::sequential)
  return(cellchat)
}

## 
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

