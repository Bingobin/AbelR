# AbelR

Personal R framework for bioinformatics analysis, statistical modeling and
scientific visualization.

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
enrichment analysis, correlation analysis, ChIP-seq, scRNA-seq, Monocle3, and
CellChat workflows.

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

AbelR currently contains 38 functions organized into nine analysis modules.
The tables provide a short overview; detailed parameter, input, and return-value
documentation should be maintained in the function-level roxygen2 comments.

### Differential expression analysis

Source: `R/deseq2.R`

| Function | Description |
|---|---|
| `DESeq2_DEG_analysis()` | Run a two-group human or mouse DESeq2 analysis using the shared species-aware implementation. |
| `DESeq2_DEG_analysis_batch()` | Run human or mouse DESeq2 analysis using the same `Group/Batch/Library` design and `GID` matrix format. |
| `DESeq2_DEG_extract()` | Extract significantly upregulated and downregulated genes and optionally draw a heatmap. |
| `Compare_pairwise_Deseq2()` | Compare fold changes from two DESeq2 result tables and highlight shared directional changes. |

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

### Enrichment and gene-set analysis

Source: `R/enrichment.R`

| Function | Description |
|---|---|
| `enrich_combind_s2()` | Unified GO and KEGG analysis with optional Reactome, Hallmark, and WikiPathways results. |
| `enricher_plot()` | Combine selected enrichment categories into a summary visualization. |
| `gsea_plot_custorm()` | Build a customized GSEA running-score plot for selected terms. |
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
| `scRNA_SCT_norm()` | Integrate a list of Seurat objects using the SCT workflow, followed by PCA, UMAP, and clustering. |
| `SCT_METHOD_V3()` | Run a Seurat v5 SCT and RPCA integration workflow with configurable graph and UMAP parameters. |
| `scTYPE_annotation()` | Annotate Seurat clusters using an explicit scTYPE directory or `AbelR.sctype_dir` option. |

### Trajectory analysis

Source: `R/monocle3.R`

| Function | Description |
|---|---|
| `run_monocle3_from_seurat_umap()` | Build a Monocle3 cell-data set using Seurat metadata and UMAP coordinates, learn a trajectory, and optionally order cells from selected root clusters. |

### Cell-cell communication

Source: `R/cellchat.R`

| Function | Description |
|---|---|
| `Build_CellChat_object()` | Build and run a human or mouse CellChat workflow from a Seurat expression layer and metadata. |

## Optional dependencies and external resources

Some modules require software or data that are not needed by the rest of the
package:

| Feature | Additional requirement |
|---|---|
| Monocle3 trajectory analysis | `monocle3` and `SingleCellExperiment` |
| CellChat analysis | `CellChat` and `future` |
| Mouse enrichment analysis | `org.Mm.eg.db` and the corresponding mouse genome resources |
| Gene-symbol-to-Entrez conversion | `MAGeCKFlute` |
| scTYPE annotation | scTYPE scripts and the scTYPE marker database |
| Genome-track visualization | `trackViewer`, genome annotations, and compatible BAM or BigWig files |
| DEG annotation and Manhattan plots | Bundled compressed gene-annotation tables, or a compatible user-supplied table |

Environment-specific resources are supplied through function arguments or
AbelR options rather than fixed paths in active package code.

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
