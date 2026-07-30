<pre align="center">
 █████╗ ██████╗ ███████╗██╗     ██████╗
██╔══██╗██╔══██╗██╔════╝██║     ██╔══██╗
███████║██████╔╝█████╗  ██║     ██████╔╝
██╔══██║██╔══██╗██╔══╝  ██║     ██╔══██╗
██║  ██║██████╔╝███████╗███████╗██║  ██║
╚═╝  ╚═╝╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═╝

Bioinformatics analysis, modelling, and visualizatio 
(•̀ᴗ•́)و  Ready to explore your data!
</pre>

<p align="center">
  Personal R framework for bioinformatics analysis, statistical modeling and
  scientific visualization.
</p>

## Development

Load the package directly from the repository during development:

```r
# When the current working directory is the AbelR repository root
devtools::load_all(export_all = FALSE)
```

`load_all()` uses the current directory (`"."`) by default. When the R console
is already running from the AbelR repository root, no path is needed. From a
different working directory, supply the package directory explicitly:

```r
devtools::load_all(path = "/path/to/AbelR", export_all = FALSE)
```

This loads the package source into the current R session without installing it.
`export_all = FALSE` respects `NAMESPACE` and exposes only AbelR's public
functions. Run it again after changing the source code to test the latest
version.

Install it locally and use it like a regular package:

```r
devtools::install()
library(AbelR)
```

The current modules cover DESeq2, DEG visualization, survival analysis,
enrichment analysis, correlation analysis, ChIP-seq, scRNA-seq, CytoTRACE2,
Monocle3, and CellChat workflows.

## Quick start

```r
library(AbelR)

# List the objects currently exported by AbelR
ls("package:AbelR")
```

The function index reflects the current public analysis interface. Superseded
compatibility wrappers are removed after their logic is consolidated into the
corresponding unified function.

## Function index

AbelR organizes its exported functions into domain-based analysis modules.
The tables provide a short overview; detailed parameter, input, and return-value
documentation is maintained in the function-level roxygen2 comments.

### Differential expression analysis

Source: `R/deseq2.R`

| Function | Description |
|---|---|
| `build_bulkRNA_batches()` | Merge batch-specific count and TPM tables and construct `CellLine/Library/Batch/Group` sample metadata. |
| `plot_bulkRNA_PCA()` | Calculate VST-based sample PCA from `build_bulkRNA_batches()` output and draw a batch-aware static or interactive plot. |
| `DESeq2_DEG_analysis()` | Run a two-group human or mouse DESeq2 analysis using the shared species-aware implementation. |
| `DESeq2_DEG_analysis_batch()` | Run human or mouse DESeq2 analysis using the same `Group/Batch/Library` design and `GID` matrix format. |
| `DESeq2_DEG_extract()` | Extract significantly upregulated and downregulated genes and optionally draw a heatmap. |
| `Compare_pairwise_Deseq2()` | Compare fold changes from two DESeq2 result tables and highlight shared directional changes. |

```r
bulk <- build_bulkRNA_batches(
  batches = c("26526", "26527"),
  target = "KAP1",
  library_name = "mRNA"
)

pca <- plot_bulkRNA_PCA(
  bulk,
  title = "PCA of KAP1 RNA-seq samples",
  top_var_n = 5000,
  interactive = TRUE
)
pca$plot
pca$heatmap
pca$interactive_plot
```

### Differential expression visualization

Source: `R/deg-plots.R`

| Function | Description |
|---|---|
| `target_for_volcano()` | Select protein-coding and user-specified genes for labelling in a DESeq2 volcano plot. |
| `volcano_plot_Deseq2()` | Draw a DESeq2 volcano plot with configurable significance thresholds and gene labels. |
| `volcano_plot_scRNA()` | Draw a volcano plot for Seurat marker results and label the most significant upregulated and downregulated genes. |
| `plot_deg_heatmap_for_DEGseq2()` | Draw a clustered DEG heatmap with optional labels for selected and top-ranked genes. |
| `plot_deg_comparison()` | Compare two DEG analyses using fold-change groups, correlation statistics, and selected gene labels. |
| `plot_deg_manhattan()` | Draw faceted Manhattan-style DEG plots based on genomic positions and signed significance scores. |

### Survival analysis

Source: `R/survival.R`

| Function | Description |
|---|---|
| `RiskScore_multivar_cox()` | Fit a multivariable Cox model and calculate normalized sample risk scores from expression and clinical features. |
| `Risk_grouping_cox()` | Fit univariable Cox models for selected genes and divide samples into high- and low-risk groups. |
| `Risk_model_df()` | Combine selected-gene expression with survival metadata to prepare a risk-model data frame. |
| `plot_survial_risk()` | Plot ranked risk scores together with patient survival status. |
| `plot_roc_curve()` | Plot time-dependent survival ROC curves from a risk-score data set. |
| `plot_logistic_forest()` | Fit a binary logistic regression model and plot coefficient-level odds ratios with confidence intervals. |

### UK Biobank analysis

Source: `R/ukb.R`

| Function | Description |
|---|---|
| `prepare_ukb_icd10_events()` | Align pipe-separated ICD-10 diagnoses with UK Biobank diagnosis-date array fields. |
| `classify_ukb_thrombosis()` | Classify ICD-10 codes into configurable thromboembolism groups and subgroups. |
| `merge_ukb_chip_callsets()` | Merge in-house and UKB-released CHIP variant calls while retaining source-specific VAF values. |
| `classify_chip_carriers()` | Create participant-level overall, gene-specific, small-clone, and large-clone CHIP indicators. |
| `make_ukb_competing_risk_data()` | Construct incident-event follow-up data with death represented as a competing event. |
| `fit_ukb_competing_risk()` | Fit crude and adjusted Fine-Gray models and calculate cumulative incidence by exposure group. |
| `run_ukb_omics_limma()` | Run covariate-adjusted limma analysis across Olink or other continuous UKB omics traits. |
| `build_ukb_multiomic_score()` | Calculate beta-weighted module scores and a combined standardized multi-omics score. |

### Enrichment and gene-set analysis

Source: `R/enrichment.R`

| Function | Description |
|---|---|
| `enrich_combind()` | Unified GO analysis with optional KEGG, Reactome, Hallmark, and WikiPathways results. |
| `enricher_plot()` | Combine selected enrichment categories into a summary visualization. |
| `gsea_plot_custorm()` | Build a customized GSEA running-score plot for selected terms. |
| `plot_gsea_dotplot()` | Plot selected or top-ranked GSEA pathways using adjusted P-value significance and normalized enrichment scores. |
| `go_plot_custom()` | Draw a customized GO enrichment plot for selected terms. |
| `GO_BP_treeplot_DESeq2()` | Run and visualize enrichment analyses separately for upregulated and downregulated DESeq2 genes. |
| `GO_BP_treeplot_scRNAseq()` | Run and visualize enrichment analyses separately for upregulated and downregulated single-cell markers. |
| `GSEA_analysis()` | Run configurable species-aware GSEA using multiple ranking metrics and gene-set databases. |

### Correlation analysis

Source: `R/correlation.R`

| Function | Description |
|---|---|
| `ccor_paired_gene()` | Test and plot the correlation between two genes in an expression matrix. |
| `ccor_paired_obj()` | Test and plot the correlation between two variables, with points coloured by sample group. |

### ChIP-seq and regulatory-element visualization

Source: `R/chipseq.R`

| Function | Description |
|---|---|
| `enhancer_rank_plot()` | Plot ranked enhancer signals and label super-enhancers or selected nearby genes. |
| `expression_boxplot_byMut()` | Compare selected-gene expression across molecular risk groups and normal samples. |
| `peak_pieplot_byhomer()` | Summarize HOMER peak annotations as a genomic-feature pie chart. |
| `peak_heatmap_byhomer()` | Draw an ordered signal heatmap from HOMER-style peak-centred matrices. |
| `track_view_cre_mut()` | Visualize a regulatory region, mutation position, and transcription-factor gene tracks. |
| `trackview_peak_roi()` | Visualize BAM or BigWig signal tracks around a selected genomic region. |
| `cre_h3k27ac_area_plot()` | Plot H3K27ac signal profiles and calculate signal area for a regulatory element. |
| `heatmap_sort()` | Sort a signal matrix using its centre-region density. |
| `ChIPseq_heatmap_plot()` | Read a ChIP-seq matrix and export a peak-centred signal heatmap. |

### Single-cell RNA-seq

Source: `R/scrna.R`

| Function | Description |
|---|---|
| `filter_seurat_qc()` | Filter a Seurat object using configurable feature, count, mitochondrial-content, and doublet thresholds. |
| `scRNA_SCT_norm()` | Run a Seurat v5 SCT and RPCA integration workflow after splitting RNA layers by a configurable sample column. |
| `run_cytotrace2_by_sample()` | Run CytoTRACE2 independently for each sample and combine the resulting cell-level metadata. |
| `run_monocle3_from_seurat_umap()` | Build a Monocle3 cell-data set using Seurat metadata and UMAP coordinates, learn a trajectory, and optionally order cells from selected root clusters. |
| `Build_CellChat_object()` | Build and run a human or mouse CellChat workflow from a Seurat expression layer and metadata. |
| `scTYPE_annotation()` | Annotate Seurat clusters using an explicit scTYPE directory or `AbelR.sctype_dir` option. |

## Optional dependencies and external resources

Some modules require software or data that are not needed by the rest of the
package:

| Feature | Additional requirement |
|---|---|
| SCT integration with `sct_method = "glmGamPoi"` | `glmGamPoi` |
| CytoTRACE2 analysis | `CytoTRACE2` |
| Monocle3 trajectory analysis | `monocle3` and `SingleCellExperiment` |
| CellChat analysis | `CellChat` and `future` |
| UKB omics association analysis | `limma` |
| UKB Fine-Gray competing-risk analysis | `cmprsk` |
| Mouse enrichment analysis | `org.Mm.eg.db` and the corresponding mouse genome resources |
| Gene-symbol-to-Entrez conversion | `MAGeCKFlute` |
| scTYPE annotation | scTYPE scripts and the scTYPE marker database |
| Genome-track visualization | `trackViewer`, genome annotations, and compatible BAM or BigWig files |
| DEG annotation and Manhattan plots | Bundled compressed gene-annotation tables, or a compatible user-supplied table |

Environment-specific resources are supplied through function arguments or
AbelR options rather than fixed paths in active package code.

## Long-running jobs on a server

The two command-line entry points under `inst/scripts/` call the same exported
AbelR functions used in an interactive R session. After installing AbelR on the
server, locate the installed scripts with:

```r
system.file("scripts", "run_sct_integrate.R", package = "AbelR")
system.file("scripts", "run_cytotrace2_by_sample.R", package = "AbelR")
```

For example, run SCT/RPCA integration from the shell with:

```bash
SCT_SCRIPT=$(Rscript -e 'cat(system.file("scripts", "run_sct_integrate.R", package = "AbelR"))')
Rscript "$SCT_SCRIPT" \
  --input_seurat input.seurat.rds \
  --out_seurat output.sct.rds \
  --sample_col SampleID \
  --maxSize_GB 120 \
  --workers 1
```

Run CytoTRACE2 separately by sample with:

```bash
CYTOTRACE_SCRIPT=$(Rscript -e 'cat(system.file("scripts", "run_cytotrace2_by_sample.R", package = "AbelR"))')
Rscript "$CYTOTRACE_SCRIPT" \
  --input_seurat output.sct.rds \
  --out_rds output.CytoTRACE2.metadata.rds \
  --sample_col SampleID \
  --species human \
  --ncores 4
```

Use `Rscript "$SCT_SCRIPT" --help` or
`Rscript "$CYTOTRACE_SCRIPT" --help` to list all thresholds and optional
arguments. By default, the CytoTRACE2 runner saves only combined metadata to
avoid retaining every large per-sample Seurat result in memory.

## Human and mouse DESeq2 configuration

`DESeq2_DEG_analysis_batch()` and `DESeq2_DEG_analysis()` support both human and
mouse RNA-seq. AbelR includes compressed annotation tables for both species and
uses them automatically:

```r
system.file(
  "extdata",
  "gene_len.v43.new.txt.gz",
  package = "AbelR"
)

system.file(
  "extdata",
  "gene_len.vM38.txt.gz",
  package = "AbelR"
)
```

Both plain `.txt` files and gzip-compressed `.txt.gz` files are supported. To
override the bundled annotations, supply `gene_anno_file` directly or configure
a path once per R session:

```r
options(
  AbelR.human_gene_anno_file = "/custom/path/gene_len.v43.new.txt.gz",
  AbelR.mouse_gene_anno_file = "/custom/path/gene_len.vM38.txt.gz"
)
```

Human and mouse analyses use exactly the same input format:

- `design.df` must contain `Group`, `Batch`, and `Library` columns.
- `count.matrix` must contain a `GID` column.
- `tpm.matrix`, when supplied, must also contain a `GID` column.

Human example:

```r
human_deg <- DESeq2_DEG_analysis_batch(
  count.matrix = human_count,
  tpm.matrix = human_tpm,
  design.df = human_design,
  library = "RNA-seq",
  batch = "Batch1",
  tr = "Treatment",
  ctr = "Control",
  species = "human"
)
```

Mouse uses the same columns and arguments; only `species` and the biological
group names change:

```r
mouse_deg <- DESeq2_DEG_analysis_batch(
  count.matrix = mouse_count,
  tpm.matrix = mouse_tpm,
  design.df = mouse_design,
  library = "RNA-seq",
  batch = "Batch1",
  tr = "DSS",
  ctr = "BLANK",
  species = "mouse"
)
```

Alternatively, `tr_filter` and `ctr_filter` can be named lists containing any
columns in `design.df`, while the required standard columns remain present. The
annotation table must use gene IDs in its first column; a `Symbol` column is
recommended so that AbelR can add Entrez IDs using the selected species.
