# Functions for enrichment analyses.

#' Run species-aware over-representation analysis
#'
#' Performs GO Biological Process, Cellular Component, and Molecular Function
#' enrichment for Entrez genes, with optional KEGG, Reactome, Hallmark, and
#' WikiPathways analyses.
#'
#' @param gene Character vector of Entrez gene identifiers.
#' @param pvc P-value cutoff passed to enrichment functions.
#' @param qvc Q-value cutoff passed to enrichment functions.
#' @param universe Optional character vector of background Entrez identifiers.
#' @param species Species name or code accepted by AbelR, such as `"human"`,
#'   `"hsa"`, `"mouse"`, or `"mmu"`.
#' @param kegg Logical; include KEGG pathway enrichment.
#' @param reactome Logical; include Reactome pathway enrichment.
#' @param hallmark Logical; include MSigDB Hallmark enrichment.
#' @param wikipathways Logical; include MSigDB WikiPathways enrichment.
#' @param minGSSize,maxGSSize Minimum and maximum gene-set sizes.
#'
#' @return A named list containing `ego_bp`, `ego_cc`, `ego_mf`, `ekg`,
#'   `reactome`, `hallmark`, and `wikipathways`; disabled or empty optional
#'   analyses are returned as `NULL`.
#' @export
enrich_combind <- function(
  gene,
  pvc = 1,
  qvc = 1,
  universe = NULL,
  species = "human",
  kegg = TRUE,
  reactome = TRUE,
  hallmark = TRUE,
  wikipathways = FALSE,
  minGSSize = 10,
  maxGSSize = 500
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required for enrichment analysis.")
  }

  ## -------- species setup --------
  species <- .abel_normalize_species(species)
  if (species == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is required for human enrichment.")
    }
    current_OrgDb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    kegg_org <- "hsa"
    reactome_org <- "human"
    msig_species <- "Homo sapiens"
  } else {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package 'org.Mm.eg.db' is required for mouse enrichment.")
    }
    current_OrgDb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    kegg_org <- "mmu"
    reactome_org <- "mouse"
    msig_species <- "Mus musculus"
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
  ego_bp <- clusterProfiler::enrichGO(
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

  ego_cc <- clusterProfiler::enrichGO(
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

  ego_mf <- clusterProfiler::enrichGO(
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
  ekg <- NULL
  if (isTRUE(kegg)) {
    ekg <- clusterProfiler::enrichKEGG(
      gene = gene,
      organism = kegg_org,
      pvalueCutoff = pvc,
      qvalueCutoff = qvc,
      universe = universe,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
    if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0) {
      ekg <- clusterProfiler::setReadable(
        ekg,
        current_OrgDb,
        keyType = "ENTREZID"
      )
    } else {
      ekg <- NULL
    }
  }

  ## -------- Reactome --------
  erct <- NULL
  if (isTRUE(reactome)) {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      stop("Package 'ReactomePA' is required when reactome = TRUE.")
    }
    erct <- ReactomePA::enrichPathway(
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
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Package 'msigdbr' is required when hallmark = TRUE.")
    }

    # msigdbr returns gene symbols; convert to ENTREZID to match 'gene'
    msig_h <- msigdbr::msigdbr(species = msig_species, category = "H") %>%
      dplyr::select(gs_name, gene_symbol) %>%
      distinct()

    # SYMBOL -> ENTREZ
    sym2ent <- clusterProfiler::bitr(
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

    ehall <- clusterProfiler::enricher(
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
      ehall <- clusterProfiler::setReadable(
        ehall,
        current_OrgDb,
        keyType = "ENTREZID"
      )
    } else {
      ehall <- NULL
    }
  }

  ## -------- MSigDB WikiPathways --------
  ewp <- NULL
  if (isTRUE(wikipathways)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Package 'msigdbr' is required when wikipathways = TRUE.")
    }
    wp <- try(
      msigdbr::msigdbr(
        species = msig_species,
        collection = "C2",
        subcollection = "CP:WIKIPATHWAYS"
      ),
      silent = TRUE
    )
    if (inherits(wp, "try-error")) {
      wp <- msigdbr::msigdbr(
        species = msig_species,
        category = "C2",
        subcategory = "CP:WIKIPATHWAYS"
      )
    }
    entrez_col <- if ("ncbi_gene" %in% colnames(wp)) {
      "ncbi_gene"
    } else {
      "entrez_gene"
    }
    term2gene_wp <- unique(data.frame(
      term = wp$gs_name,
      gene = as.character(wp[[entrez_col]])
    ))
    term2gene_wp <- term2gene_wp[!is.na(term2gene_wp$gene), ]
    ewp <- clusterProfiler::enricher(
      gene = gene,
      TERM2GENE = term2gene_wp,
      universe = universe,
      pAdjustMethod = "BH",
      pvalueCutoff = pvc,
      qvalueCutoff = qvc,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
    if (!is.null(ewp) && nrow(as.data.frame(ewp)) > 0) {
      ewp <- clusterProfiler::setReadable(
        ewp,
        current_OrgDb,
        keyType = "ENTREZID"
      )
    } else {
      ewp <- NULL
    }
  }

  ## -------- output --------
  out <- list(
    ego_bp = ego_bp,
    ego_cc = ego_cc,
    ego_mf = ego_mf,
    ekg = ekg,
    reactome = erct,
    hallmark = ehall,
    wikipathways = ewp
  )
  return(out)
}


#' Run MSigDB over-representation analysis
#'
#' Performs over-representation analysis for one or more selected native Human
#' or Mouse MSigDB collections using gene symbols. Human and mouse use distinct
#' official collection names.
#'
#' @param gene Character vector of human or mouse gene symbols.
#' @param database Character vector naming one or more species-specific MSigDB
#'   collections or subcollections. If `NULL`, `"H"` is used for human and
#'   `"MH"` for mouse. See Details for all supported names.
#' @param species Species name or code accepted by AbelR: `"human"`, `"hsa"`,
#'   `"Homo sapiens"`, `"mouse"`, `"mmu"`, or `"Mus musculus"`.
#' @param universe Optional character vector of background gene symbols from
#'   the same species as `gene`.
#' @param pvalueCutoff,qvalueCutoff P-value and q-value cutoffs passed to
#'   [clusterProfiler::enricher()].
#' @param pAdjustMethod Multiple-testing correction method.
#' @param minGSSize,maxGSSize Minimum and maximum gene-set sizes.
#' @param return_object Logical; return `enrichResult` objects instead of data
#'   frames.
#'
#' @details
#' **Human MSigDB (`species = "human"`)**
#'
#' - `H`: Hallmark gene sets representing 50 coherent biological states and
#'   processes; a useful first-pass overview.
#' - `C1`: genes grouped by human chromosome cytogenetic bands; useful for
#'   detecting regional or copy-number-associated signals.
#' - `C2`, `C2_CGP`, `C2_CP`: curated gene sets, chemical/genetic perturbation
#'   signatures, or canonical pathways, respectively.
#' - `C2_CP_BIOCARTA`, `C2_CP_KEGG_MEDICUS`, `C2_CP_KEGG_LEGACY`,
#'   `C2_CP_PID`, `C2_CP_REACTOME`, `C2_CP_WIKIPATHWAYS`: pathway collections
#'   from the named resources. `KEGG_LEGACY`, BioCarta, and PID are mainly
#'   useful for comparison with older analyses.
#' - `C3`, `C3_MIR_MIRDB`, `C3_MIR_MIR_LEGACY`, `C3_TFT_GTRD`,
#'   `C3_TFT_TFT_LEGACY`: predicted microRNA or transcription-factor targets;
#'   useful for generating regulatory hypotheses, not proving direct binding.
#' - `C4`, `C4_3CA`, `C4_CGN`, `C4_CM`: computational cancer gene sets,
#'   including cancer-cell metaprograms, cancer-gene neighbourhoods, and cancer
#'   modules.
#' - `C5`, `C5_GO_BP`, `C5_GO_CC`, `C5_GO_MF`, `C5_HPO`: ontology gene sets for
#'   biological processes, cellular components, molecular functions, or human
#'   disease phenotypes.
#' - `C6`: oncogenic perturbation signatures; useful for relating genes to
#'   cancer-associated pathway activity.
#' - `C7`, `C7_IMMUNESIGDB`, `C7_VAX`: immune states, immune perturbations, or
#'   vaccine-response signatures.
#' - `C8`: curated human cell-type marker signatures from single-cell studies.
#' - `C9`: computational perturbation signatures inferred from DepMap CRISPR
#'   dependency and CCLE expression profiles.
#'
#' **Mouse MSigDB (`species = "mouse"`)**
#'
#' - `MH`: mouse-ortholog Hallmark gene sets; a useful first-pass overview.
#' - `M1`: genes grouped by mouse chromosome cytogenetic bands.
#' - `M2`, `M2_CGP`: mouse-native curated sets or perturbation signatures.
#' - `M2_CP_BIOCARTA`, `M2_CP_REACTOME`, `M2_CP_WIKIPATHWAYS`: mouse canonical
#'   pathway subcollections. Mouse MSigDB does not provide the human KEGG or PID
#'   subcollections listed above.
#' - `M3`, `M3_MIRDB`, `M3_GTRD`: predicted mouse microRNA or
#'   transcription-factor targets.
#' - `M5`, `M5_GO_BP`, `M5_GO_CC`, `M5_GO_MF`: mouse Gene Ontology sets.
#' - `M5_MPT`: cancer-related terms from the Mammalian Phenotype Ontology.
#' - `M7`: mouse immune-cell states and immune perturbation signatures.
#' - `M8`: mouse cell-type marker signatures from single-cell studies.
#'
#' An enriched set indicates overlap with the supplied genes. It does not by
#' itself establish pathway activation, direct regulation, or altered cell-type
#' abundance. Interpret results with expression direction and experimental
#' context. Missing dependencies, database-loading failures, and empty
#' enrichment results generate warnings and return `NULL` rather than stopping
#' a batch workflow. Invalid species, database names, and parameter values still
#' generate errors.
#'
#' @return For one database, a data frame or `enrichResult` object. If the
#'   analysis cannot be run or has no enriched terms, a warning is issued and
#'   `NULL` is returned. For multiple databases, a named list is returned and
#'   failed or empty analyses are represented by `NULL` entries.
#' @export
.abel_msigdb_database_map <- function(species) {
  species <- .abel_normalize_species(species)

  if (species == "human") {
    return(list(
      H = c("H", NA_character_),
      C1 = c("C1", NA_character_),
      C2 = c("C2", NA_character_),
      C2_CGP = c("C2", "CGP"),
      C2_CP = c("C2", "CP"),
      C2_CP_BIOCARTA = c("C2", "CP:BIOCARTA"),
      C2_CP_KEGG_LEGACY = c("C2", "CP:KEGG_LEGACY"),
      C2_CP_KEGG_MEDICUS = c("C2", "CP:KEGG_MEDICUS"),
      C2_CP_PID = c("C2", "CP:PID"),
      C2_CP_REACTOME = c("C2", "CP:REACTOME"),
      C2_CP_WIKIPATHWAYS = c("C2", "CP:WIKIPATHWAYS"),
      C3 = c("C3", NA_character_),
      C3_MIR_MIRDB = c("C3", "MIR:MIRDB"),
      C3_MIR_MIR_LEGACY = c("C3", "MIR:MIR_LEGACY"),
      C3_TFT_GTRD = c("C3", "TFT:GTRD"),
      C3_TFT_TFT_LEGACY = c("C3", "TFT:TFT_LEGACY"),
      C4 = c("C4", NA_character_),
      C4_3CA = c("C4", "3CA"),
      C4_CGN = c("C4", "CGN"),
      C4_CM = c("C4", "CM"),
      C5 = c("C5", NA_character_),
      C5_GO_BP = c("C5", "GO:BP"),
      C5_GO_CC = c("C5", "GO:CC"),
      C5_GO_MF = c("C5", "GO:MF"),
      C5_HPO = c("C5", "HPO"),
      C6 = c("C6", NA_character_),
      C7 = c("C7", NA_character_),
      C7_IMMUNESIGDB = c("C7", "IMMUNESIGDB"),
      C7_VAX = c("C7", "VAX"),
      C8 = c("C8", NA_character_),
      C9 = c("C9", NA_character_)
    ))
  }

  list(
    MH = c("MH", NA_character_),
    M1 = c("M1", NA_character_),
    M2 = c("M2", NA_character_),
    M2_CGP = c("M2", "CGP"),
    M2_CP_BIOCARTA = c("M2", "CP:BIOCARTA"),
    M2_CP_REACTOME = c("M2", "CP:REACTOME"),
    M2_CP_WIKIPATHWAYS = c("M2", "CP:WIKIPATHWAYS"),
    M3 = c("M3", NA_character_),
    M3_MIRDB = c("M3", "MIRDB"),
    M3_GTRD = c("M3", "GTRD"),
    M5 = c("M5", NA_character_),
    M5_GO_BP = c("M5", "GO:BP"),
    M5_GO_CC = c("M5", "GO:CC"),
    M5_GO_MF = c("M5", "GO:MF"),
    M5_MPT = c("M5", "MPT"),
    M7 = c("M7", NA_character_),
    M8 = c("M8", NA_character_)
  )
}

.abel_match_msigdb_database <- function(database, species) {
  species <- .abel_normalize_species(species)
  db_map <- .abel_msigdb_database_map(species)
  database <- gsub(
    ":",
    "_",
    toupper(trimws(as.character(database))),
    fixed = TRUE
  )
  if (length(database) != 1L || is.na(database) || !nzchar(database)) {
    return(NA_character_)
  }
  if (database %in% names(db_map)) {
    return(database)
  }

  if (species == "mouse") {
    legacy_aliases <- c(
      M2_CP = "M2",
      M3_MIR_MIRDB = "M3_MIRDB",
      M3_TFT_GTRD = "M3_GTRD"
    )
    if (database %in% names(legacy_aliases)) {
      return(unname(legacy_aliases[[database]]))
    }
  }

  NA_character_
}

EnrichMSigDB <- function(
  gene,
  database = NULL,
  species = "human",
  universe = NULL,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  return_object = FALSE
) {
  species <- .abel_normalize_species(species)
  msig_species <- switch(species,
    human = "Homo sapiens",
    mouse = "Mus musculus"
  )
  db_species <- if (species == "human") "HS" else "MM"
  db_map <- .abel_msigdb_database_map(species)

  if (is.null(database)) {
    database <- if (species == "human") "H" else "MH"
  }
  if (!is.character(database) || !length(database) || anyNA(database)) {
    stop("database must contain at least one supported database name.")
  }
  database <- toupper(trimws(database))
  resolved_database <- vapply(
    database,
    .abel_match_msigdb_database,
    character(1),
    species = species,
    USE.NAMES = FALSE
  )
  unsupported <- database[is.na(resolved_database) | !nzchar(resolved_database)]
  if (length(unsupported)) {
    stop(
      "Unsupported ", species, " database: ",
      paste(unsupported, collapse = ", "),
      "\nSupported ", species, " databases: ",
      paste(names(db_map), collapse = ", ")
    )
  }
  database <- resolved_database

  gene <- unique(as.character(gene))
  gene <- gene[!is.na(gene) & nzchar(gene)]
  if (!length(gene)) {
    warning("No valid genes were provided; returning NULL.", call. = FALSE)
    return(NULL)
  }
  if (!is.null(universe)) {
    universe <- unique(as.character(universe))
    universe <- universe[!is.na(universe) & nzchar(universe)]
    if (!length(universe)) {
      warning(
        "universe contains no valid genes; using the gene-set database ",
        "background instead.",
        call. = FALSE
      )
      universe <- NULL
    }
  }
  if (!is.numeric(pvalueCutoff) || length(pvalueCutoff) != 1L ||
    is.na(pvalueCutoff) || pvalueCutoff < 0 || pvalueCutoff > 1) {
    stop("pvalueCutoff must be one number between 0 and 1.")
  }
  if (!is.numeric(qvalueCutoff) || length(qvalueCutoff) != 1L ||
    is.na(qvalueCutoff) || qvalueCutoff < 0 || qvalueCutoff > 1) {
    stop("qvalueCutoff must be one number between 0 and 1.")
  }
  if (!is.numeric(minGSSize) || length(minGSSize) != 1L ||
    is.na(minGSSize) || minGSSize < 1 || minGSSize %% 1 != 0) {
    stop("minGSSize must be one positive integer.")
  }
  if (!is.numeric(maxGSSize) || length(maxGSSize) != 1L ||
    is.na(maxGSSize) || maxGSSize < minGSSize || maxGSSize %% 1 != 0) {
    stop("maxGSSize must be an integer greater than or equal to minGSSize.")
  }
  if (!is.logical(return_object) || length(return_object) != 1L ||
    is.na(return_object)) {
    stop("return_object must be TRUE or FALSE.")
  }

  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    warning("Package 'msigdbr' is required; returning NULL.", call. = FALSE)
    return(NULL)
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    warning(
      "Package 'clusterProfiler' is required; returning NULL.",
      call. = FALSE
    )
    return(NULL)
  }
  msigdbr_args <- names(formals(msigdbr::msigdbr))
  if (species == "mouse" && !"db_species" %in% msigdbr_args) {
    warning(
      "Native Mouse MSigDB requires msigdbr >= 10.0.0. ",
      "Please update the 'msigdbr' package; returning NULL.",
      call. = FALSE
    )
    return(NULL)
  }

  load_msigdb <- function(collection, subcollection) {
    if ("collection" %in% msigdbr_args) {
      args <- list(
        db_species = db_species,
        species = msig_species,
        collection = collection
      )
      if (!is.na(subcollection)) {
        args$subcollection <- subcollection
      }
    } else {
      args <- list(
        species = msig_species,
        category = collection
      )
      if (!is.na(subcollection)) {
        args$subcategory <- subcollection
      }
    }
    do.call(msigdbr::msigdbr, args)
  }

  run_one <- function(db) {
    tryCatch(
      {
        collection <- unname(db_map[[db]][1])
        subcollection <- unname(db_map[[db]][2])
        message("Running ", species, " MSigDB enrichment: ", db)

        msig <- load_msigdb(collection, subcollection)
        required_cols <- c("gs_name", "gene_symbol")
        missing_cols <- setdiff(required_cols, colnames(msig))
        if (length(missing_cols)) {
          stop(
            "msigdbr result is missing required columns: ",
            paste(missing_cols, collapse = ", ")
          )
        }
        if (!nrow(msig)) {
          warning(
            "No gene sets were found for ", db, "; returning NULL.",
            call. = FALSE
          )
          return(NULL)
        }

        term2gene <- unique(data.frame(
          term = as.character(msig$gs_name),
          gene = as.character(msig$gene_symbol),
          stringsAsFactors = FALSE
        ))
        if ("gs_description" %in% colnames(msig)) {
          term_name <- as.character(msig$gs_description)
        } else {
          term_name <- as.character(msig$gs_name)
        }
        missing_name <- is.na(term_name) | !nzchar(term_name)
        term_name[missing_name] <- as.character(msig$gs_name[missing_name])
        term2name <- unique(data.frame(
          term = as.character(msig$gs_name),
          name = term_name,
          stringsAsFactors = FALSE
        ))

        enrichment <- clusterProfiler::enricher(
          gene = gene,
          universe = universe,
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize
        )
        if (is.null(enrichment)) {
          warning(
            "No enrichment result was returned for ", db, "; returning NULL.",
            call. = FALSE
          )
          return(NULL)
        }
        if (return_object) {
          return(enrichment)
        }

        result <- as.data.frame(enrichment)
        if (!nrow(result)) {
          warning(
            "No enriched terms were found for ", db, "; returning NULL.",
            call. = FALSE
          )
          return(NULL)
        }
        result$Database <- db
        result$Species <- species
        result <- result[, c(
          "Database",
          "Species",
          setdiff(colnames(result), c("Database", "Species"))
        ), drop = FALSE]
        result
      },
      error = function(e) {
        warning(
          "MSigDB enrichment failed for ", db, ": ",
          conditionMessage(e),
          "; returning NULL.",
          call. = FALSE
        )
        NULL
      }
    )
  }

  result <- lapply(database, run_one)
  names(result) <- database
  if (length(database) == 1L) {
    return(result[[1]])
  }
  result
}


#' Summarize multiple enrichment categories
#'
#' Selects rows from GO, WikiPathways, and KEGG enrichment results and combines
#' them in a faceted significance plot.
#'
#' @param enricher Enrichment list returned by [enrich_combind()].
#' @param bp,cc,mf,wp,kg Integer row indices selected from GO BP, GO CC, GO MF,
#'   WikiPathways, and KEGG results.
#' @param value Column index in each enrichment result used as the plotted
#'   P-value-like quantity.
#'
#' @return A [ggplot2::ggplot] enrichment summary.
#' @export
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
    Description = enricher$wikipathways@result[wp, 2],
    Pvalue = enricher$wikipathways@result[wp, value],
    Type = rep("WikiPath", length(wp))
  )
  df <- rbind(df, df.tmp)
  df.tmp <- data.frame(
    Description = enricher$ekg@result[kg, 2],
    Pvalue = enricher$ekg@result[kg, value],
    Type = rep("KEGG", length(kg))
  )
  df <- rbind(df, df.tmp)
  # ggplot(df,aes(x=Description,y=-log(Pvalue),fill=Type)) + geom_bar(stat = "identity",show.legend = TRUE) + coord_flip()
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
    facet_grid(Type ~ ., scales = "free", space = "free_y", switch = "x") +
    blank +
    ylab("")
  return(p)
}


#' Draw a customized GSEA running-score plot
#'
#' Builds a three-panel GSEA plot showing the enrichment score, ranked gene
#' positions, and ranked-list metric for one selected gene set.
#'
#' @param gsea_ob A GSEA result object compatible with `enrichplot`.
#' @param select_term Row index or term selection passed to the internal GSEA
#'   plotting data extractor.
#' @param color Colour used for the enrichment-score line and annotation.
#' @param xpos X coordinate used to position the NES and adjusted P-value label.
#'
#' @return A combined cowplot drawing object.
#' @export
gsea_plot_custorm <- function(gsea_ob, select_term, color, xpos = 3000) {
  # gsea_ob <- aml_phenolyzer.gsea.crc
  # select_term <- 1
  # color <- "#08537C"
  nes <- round(gsea_ob@result[select_term, "NES"], digits = 2)
  pv <- formatC(
    gsea_ob@result[select_term, "p.adjust"],
    format = "e",
    digits = 2
  )
  # pv <- round(gsea_ob@result[select_term,"p.adjust"], digits = 6)
  # pv <- round(gsea_ob@result[select_term,"pvalue"], digits = 6)
  gs_info <- utils::getFromNamespace("gsInfo", "enrichplot")
  gsdata <- do.call(
    rbind,
    lapply(select_term, gs_info, object = gsea_ob)
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
    # theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) +
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
    # annotate("text", label = paste0("P value = ", pv), x=xpos,y= ypos -0.05) +
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
    # geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax),color=color, show.legend = FALSE,alpha =0.6,size = 0.4) +
    blank +
    xlab(NULL) +
    theme(axis.ticks = element_blank(), axis.text = element_blank()) +
    theme(
      panel.border = element_rect(fill = NA, linetype = 1, size = 1),
      axis.line = element_blank()
    ) +
    # geom_hline(yintercept = 0) +
    theme(plot.margin = margin(t = -0.5, r = 0.2, b = 0, l = 1.2, unit = "cm"))

  p_gsea_3 <- ggplot(gsdata, aes_(x = ~x, y = ~geneList)) +
    # rasterise(geom_segment(aes_(xend = ~x, yend = 0), color = color, show.legend = FALSE)) +
    rasterise(
      geom_area(color = color, fill = color, show.legend = FALSE),
      dpi = 300
    ) +
    # geom_area(color = color, fill = color, show.legend = FALSE) +
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
  # plot_grid(p_gsea_1, p_gsea_2, nrow= 2, rel_heights = c(8,1), align = "v")

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


#' Plot grouped GSEA results with dot and heatmap views
#'
#' Selects pathways from one GSEA result or from multiple group-level results
#' generated from the same gene-set database. For multiple groups, the union of
#' independently selected top pathways is compared across groups in both a dot
#' plot and an NES heatmap. Both plots and the returned table use the same
#' pathway and group order. When every selected term has the same uppercase
#' database prefix, such as `HALLMARK_`, the prefix is moved to the y-axis title
#' and removed from the displayed term labels.
#'
#' @param gsea_result A `clusterProfiler` GSEA result object, a data frame, or a
#'   list of group-level results from the same gene-set database. Each result
#'   must contain `Description`, `NES`, `pvalue`, and `p.adjust`. List names are
#'   used as group labels; unnamed elements are labelled `Group1`, `Group2`, and
#'   so on with a warning.
#' @param top_n Number of pathways selected independently from each group using
#'   the smallest adjusted P values, or raw P values when adjusted P values are
#'   all missing. Supply one number for every group, or a named numeric vector
#'   to set different values for different groups. Set to `0` or `NULL` to use
#'   only `label_pathways`.
#' @param label_pathways Optional character vector of pathway descriptions to
#'   include. These pathways are placed first and displayed from top to bottom
#'   in exactly the supplied order. Automatically selected top pathways not
#'   already listed are appended below them.
#' @param group_order Optional character vector specifying the left-to-right
#'   order of groups. The input list order is used by default.
#' @param sentence_case Logical. If `TRUE`, replace hyphens and underscores in
#'   pathway labels with spaces, convert all letters to lowercase, and capitalize
#'   only the first letter. Raw pathway identifiers in the returned table remain
#'   unchanged.
#' @param low_color,mid_color,high_color Colours representing negative, zero,
#'   and positive NES values in both plots.
#' @param na_color Fill colour for missing NES values in the heatmap.
#' @param star_size Text size for heatmap significance stars.
#' @param tile_color,tile_size Heatmap tile border colour and line width.
#'
#' @return A named list containing `dotplot`, `heatmap`, and `table`. The table
#'   contains one row for every selected pathway and group combination, including
#'   combinations without a reported GSEA result.
#' @export
plot_gsea_summary <- function(
  gsea_result,
  top_n = 20,
  label_pathways = NULL,
  group_order = NULL,
  sentence_case = FALSE,
  low_color = "#08306B",
  mid_color = "grey90",
  high_color = "#67000D",
  na_color = "grey85",
  star_size = 4,
  tile_color = "grey90",
  tile_size = 0.35
) {
  required_cols <- c("Description", "NES", "pvalue", "p.adjust")
  if (!is.logical(sentence_case) || length(sentence_case) != 1L ||
    is.na(sentence_case)) {
    stop("sentence_case must be TRUE or FALSE.")
  }
  label_pathways <- unique(as.character(label_pathways))
  label_pathways <- label_pathways[
    !is.na(label_pathways) & nzchar(label_pathways)
  ]

  coerce_result <- function(x, group_name, allow_skip = FALSE) {
    result <- tryCatch(
      as.data.frame(x),
      error = function(e) {
        if (allow_skip) {
          warning(
            "Ignoring group '", group_name, "': ", conditionMessage(e),
            call. = FALSE
          )
          return(NULL)
        }
        stop("gsea_result must be coercible to a data frame: ", conditionMessage(e))
      }
    )
    if (is.null(result)) {
      return(NULL)
    }
    missing_cols <- setdiff(required_cols, colnames(result))
    if (length(missing_cols)) {
      if (allow_skip) {
        warning(
          "Ignoring group '", group_name, "' because it is missing: ",
          paste(missing_cols, collapse = ", "),
          call. = FALSE
        )
        return(NULL)
      }
      stop(
        "gsea_result is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    result
  }

  is_result_list <- is.list(gsea_result) && !is.data.frame(gsea_result)
  if (is_result_list) {
    result_list <- gsea_result
    if (!length(result_list)) {
      stop("gsea_result contains no results.")
    }
    group_names <- names(result_list)
    if (is.null(group_names)) {
      group_names <- rep("", length(result_list))
    }
    blank_names <- is.na(group_names) | !nzchar(group_names)
    if (any(blank_names)) {
      warning(
        "gsea_result has unnamed groups; assigning Group1, Group2, and so on.",
        call. = FALSE
      )
      group_names[blank_names] <- paste0("Group", which(blank_names))
    }
    if (anyDuplicated(group_names)) {
      warning(
        "gsea_result contains duplicated group names; making them unique.",
        call. = FALSE
      )
      group_names <- make.unique(group_names)
    }
    names(result_list) <- group_names

    null_results <- vapply(result_list, is.null, logical(1))
    if (any(null_results)) {
      warning(
        "Ignoring NULL GSEA groups: ",
        paste(names(result_list)[null_results], collapse = ", "),
        call. = FALSE
      )
      result_list <- result_list[!null_results]
    }
    if (!length(result_list)) {
      stop("gsea_result contains only NULL results.")
    }
  } else {
    result_list <- list(GSEA = gsea_result)
  }

  group_names <- names(result_list)
  result_list <- lapply(
    group_names,
    function(group_name) {
      coerce_result(
        result_list[[group_name]],
        group_name,
        allow_skip = is_result_list
      )
    }
  )
  names(result_list) <- group_names
  invalid_results <- vapply(result_list, is.null, logical(1))
  if (any(invalid_results)) {
    result_list <- result_list[!invalid_results]
  }
  if (!length(result_list)) {
    stop("gsea_result contains no valid GSEA results.")
  }

  if (!is.null(group_order)) {
    if (!is.character(group_order) || anyNA(group_order) ||
      any(!nzchar(group_order)) || anyDuplicated(group_order)) {
      stop("group_order must contain unique, non-empty group names.")
    }
    unknown_groups <- setdiff(group_order, names(result_list))
    if (length(unknown_groups)) {
      stop(
        "group_order contains unknown groups: ",
        paste(unknown_groups, collapse = ", ")
      )
    }
    omitted_groups <- setdiff(names(result_list), group_order)
    group_order <- c(group_order, omitted_groups)
    result_list <- result_list[group_order]
  }
  group_names <- names(result_list)

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || !length(top_n) || anyNA(top_n) ||
      any(!is.finite(top_n)) || any(top_n < 0)) {
      stop("top_n must be NULL or contain non-negative numbers.")
    }
    if (length(top_n) == 1L) {
      top_by_group <- stats::setNames(
        rep(as.integer(top_n), length(result_list)),
        group_names
      )
    } else if (!is.null(names(top_n)) &&
      all(group_names %in% names(top_n))) {
      top_by_group <- as.integer(top_n[group_names])
      names(top_by_group) <- group_names
    } else if (length(top_n) == length(result_list)) {
      top_by_group <- as.integer(top_n)
      names(top_by_group) <- group_names
    } else {
      stop(
        "top_n must have length 1, match the number of groups, or be named ",
        "for every group."
      )
    }
  } else {
    top_by_group <- stats::setNames(rep(0L, length(result_list)), group_names)
  }

  rank_result <- function(result) {
    rank_col <- if (any(!is.na(result$p.adjust))) "p.adjust" else "pvalue"
    ranked_result <- result[
      !is.na(result[[rank_col]]) & is.finite(result[[rank_col]]), ,
      drop = FALSE
    ]
    ranked_result <- ranked_result[
      order(ranked_result[[rank_col]], -abs(ranked_result$NES)), ,
      drop = FALSE
    ]
    ranked_result <- ranked_result[
      !duplicated(as.character(ranked_result$Description)),
      required_cols,
      drop = FALSE
    ]
    ranked_result
  }

  ranked_list <- lapply(
    group_names,
    function(group_name) {
      rank_result(result_list[[group_name]])
    }
  )
  names(ranked_list) <- group_names
  empty_groups <- vapply(ranked_list, nrow, integer(1)) == 0L
  if (any(empty_groups)) {
    warning(
      "Ignoring groups without finite pathway P values: ",
      paste(group_names[empty_groups], collapse = ", "),
      call. = FALSE
    )
    ranked_list <- ranked_list[!empty_groups]
    group_names <- group_names[!empty_groups]
    top_by_group <- top_by_group[group_names]
  }
  if (!length(ranked_list)) {
    stop("No groups contain finite pathway P values.")
  }

  top_terms <- unlist(
    lapply(group_names, function(group_name) {
      utils::head(
        as.character(ranked_list[[group_name]]$Description),
        top_by_group[[group_name]]
      )
    }),
    use.names = FALSE
  )
  available_terms <- unique(unlist(
    lapply(ranked_list, function(x) as.character(x$Description)),
    use.names = FALSE
  ))
  missing_labels <- setdiff(label_pathways, available_terms)
  if (length(missing_labels)) {
    warning(
      "label_pathways not found in any group: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }
  ordered_labels <- label_pathways[label_pathways %in% available_terms]
  pathway_order <- unique(c(ordered_labels, top_terms))
  if (!length(pathway_order)) {
    stop(
      "No pathways were selected. Use top_n > 0 or provide label_pathways."
    )
  }

  has_database_prefix <- grepl("^[A-Z][A-Z0-9]*_", pathway_order)
  database_prefix <- sub("_.*$", "", pathway_order)
  common_database <- if (
    all(has_database_prefix) && length(unique(database_prefix)) == 1L
  ) {
    unique(database_prefix)
  } else {
    NULL
  }
  if (is.null(common_database)) {
    display_order <- pathway_order
    pathway_axis_title <- "Pathway"
  } else {
    display_order <- sub(
      paste0("^", common_database, "_+"),
      "",
      pathway_order
    )
    pathway_axis_title <- paste(common_database, "pathway")
  }
  if (sentence_case) {
    display_order <- gsub("[-_]+", " ", display_order)
    display_order <- gsub("[[:space:]]+", " ", trimws(display_order))
    display_order <- tolower(display_order)
    display_order <- paste0(
      toupper(substr(display_order, 1L, 1L)),
      substring(display_order, 2L)
    )
  }
  display_term <- stats::setNames(display_order, pathway_order)

  plot_list <- lapply(group_names, function(group_name) {
    group_result <- ranked_list[[group_name]]
    group_result <- group_result[
      as.character(group_result$Description) %in% pathway_order,
      required_cols,
      drop = FALSE
    ]
    if (!nrow(group_result)) {
      return(NULL)
    }
    group_result$Group <- group_name
    group_result
  })
  gsea.plot <- do.call(rbind, plot_list)
  rownames(gsea.plot) <- NULL
  gsea.plot$Description <- as.character(gsea.plot$Description)
  gsea.plot <- dplyr::left_join(
    tidyr::expand_grid(
      Description = pathway_order,
      Group = group_names
    ),
    gsea.plot,
    by = c("Description", "Group")
  )
  gsea.plot$Direction <- ifelse(
    is.na(gsea.plot$NES),
    NA_character_,
    ifelse(gsea.plot$NES < 0, "Down", "Up")
  )
  gsea.plot$PlotP <- ifelse(
    is.na(gsea.plot$p.adjust),
    gsea.plot$pvalue,
    gsea.plot$p.adjust
  )
  gsea.plot$PlotP <- pmax(gsea.plot$PlotP, .Machine$double.xmin)
  gsea.plot$Significance <- dplyr::case_when(
    is.na(gsea.plot$PlotP) ~ "",
    gsea.plot$PlotP < 1e-4 ~ "****",
    gsea.plot$PlotP < 1e-3 ~ "***",
    gsea.plot$PlotP < 1e-2 ~ "**",
    gsea.plot$PlotP < 5e-2 ~ "*",
    TRUE ~ ""
  )
  gsea.plot$Description <- factor(
    gsea.plot$Description,
    levels = rev(pathway_order)
  )
  gsea.plot$DisplayTerm <- factor(
    unname(display_term[as.character(gsea.plot$Description)]),
    levels = rev(display_order)
  )
  gsea.plot$Group <- factor(
    gsea.plot$Group,
    levels = group_names
  )

  dotplot <- ggplot2::ggplot(
    gsea.plot,
    ggplot2::aes(
      x = .data[["Group"]],
      y = .data[["DisplayTerm"]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = -log10(.data[["PlotP"]]),
        color = .data[["NES"]]
      ),
      shape = 16,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = -log10(.data[["PlotP"]])),
      shape = 1,
      color = "black",
      na.rm = TRUE
    ) +
    ggplot2::theme_test() +
    ggplot2::scale_color_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color
    ) +
    ggplot2::labs(
      x = "GSEA",
      y = pathway_axis_title,
      size = "-log10(P)"
    )
  if (length(group_names) == 1L) {
    dotplot <- dotplot + ggplot2::theme(
      axis.text.x = ggplot2::element_blank()
    )
  } else {
    dotplot <- dotplot + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  }

  heatmap <- ggplot2::ggplot(
    gsea.plot,
    ggplot2::aes(
      x = .data[["Group"]],
      y = .data[["DisplayTerm"]],
      fill = .data[["NES"]]
    )
  ) +
    ggplot2::geom_tile(color = tile_color, linewidth = tile_size) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["Significance"]]),
      size = star_size,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      na.value = na_color,
      name = "NES"
    ) +
    ggplot2::labs(x = "GSEA", y = pathway_axis_title) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.ticks = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    )

  list(
    dotplot = dotplot,
    heatmap = heatmap,
    table = gsea.plot
  )
}


#' Draw a customized GO enrichment dot plot
#'
#' Selects terms from an enrichment result, calculates gene counts, caps very
#' small adjusted P values, and plots term significance and size.
#'
#' @param Enricher Enrichment result data frame or an object supporting data
#'   frame-style column access.
#' @param select_term Row indices of terms to display.
#' @param color Name of an RColorBrewer sequential palette.
#' @param Label Label displayed on the plot's categorical axis.
#' @param pmin Minimum adjusted P value used when plotting `-log10` values.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
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


#' Plot enrichment trees for DESeq2 gene sets
#'
#' Splits an annotated DESeq2 result into upregulated, downregulated, and all
#' differential genes, runs human GO and KEGG enrichment, and creates semantic
#' similarity tree plots for each set.
#'
#' @param deg.df Annotated DESeq2 result containing `Entrez`, `padj`, and
#'   `log2FoldChange`.
#' @param label Character label identifying the comparison.
#' @param pv Adjusted P-value threshold.
#' @param lfc Absolute log2 fold-change threshold.
#'
#' @return A list containing gene counts, enrichment objects for each direction,
#'   the comparison label, and a named list of plots.
#' @export
GO_BP_treeplot_DESeq2 <- function(deg.df, label, pv = 0.05, lfc = log2(1.5)) {
  #  deg.df <- HMGA2_sh2.DESeq2.result$result
  #  label <- "shHMGA2"
  up.entrez <- (deg.df %>% filter(padj < pv, log2FoldChange > lfc))$Entrez
  dw.entrez <- (deg.df %>% filter(padj < pv, log2FoldChange < -lfc))$Entrez
  up.entrez <- unique(na.omit(up.entrez))
  dw.entrez <- unique(na.omit(dw.entrez))
  deg.entrez <- unique(c(up.entrez, dw.entrez))
  color <- "YlGnBu"
  up.enricher <- enrich_combind(
    up.entrez,
    species = "human",
    reactome = FALSE,
    hallmark = FALSE
  )
  dw.enricher <- enrich_combind(
    dw.entrez,
    species = "human",
    reactome = FALSE,
    hallmark = FALSE
  )
  deg.enricher <- enrich_combind(
    deg.entrez,
    species = "human",
    reactome = FALSE,
    hallmark = FALSE
  )
  ################## go_bp
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
        "Up-regulated Genes in GO_BP (",
        label,
        " n= ",
        length(up.entrez),
        ")"
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
        "Down-regulated Genes in GO_BP (",
        label,
        " n= ",
        length(dw.entrez),
        ")"
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
        "Diff. Expr. Genes in GO_BP (",
        label,
        " n= ",
        length(deg.entrez),
        ")"
      )
    )
  ################## kegg
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
        "Up-regulated Genes in KEGG pathway (",
        label,
        " n= ",
        length(up.entrez),
        ")"
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
        "Down-regulated Genes in KEGG pathway (",
        label,
        " n= ",
        length(dw.entrez),
        ")"
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
        "Diff. Expr. Genes in KEGG pathway (",
        label,
        " n= ",
        length(deg.entrez),
        ")"
      )
    )

  ###################

  p_combind_up <- enricher_plot(up.enricher) +
    ggtitle(
      label = paste0(
        "Up-regulated Genes Enrichment Analysis(",
        label,
        " n= ",
        length(up.entrez),
        ")"
      )
    )
  p_combind_dw <- enricher_plot(dw.enricher) +
    ggtitle(
      label = paste0(
        "Down-regulated Genes Enrichment Analysis(",
        label,
        " n= ",
        length(dw.entrez),
        ")"
      )
    )
  p_combind_deg <- enricher_plot(deg.enricher) +
    ggtitle(
      label = paste0(
        "Diff. Expr. Genes Enrichment Analysis(",
        label,
        " n= ",
        length(deg.entrez),
        ")"
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


#' Plot enrichment trees for single-cell marker genes
#'
#' Converts human marker symbols to Entrez identifiers, separates upregulated
#' and downregulated markers, performs GO and KEGG enrichment, and creates
#' semantic similarity tree plots.
#'
#' @param deg.df Single-cell marker data frame with gene symbols as row names and
#'   `p_val` and `avg_log2FC` columns.
#' @param label Character label identifying the comparison.
#' @param pv P-value threshold.
#' @param lfc Absolute average log2 fold-change threshold.
#'
#' @return A list containing gene counts, directional enrichment objects, the
#'   comparison label, and a named list of plots.
#' @export
GO_BP_treeplot_scRNAseq <- function(deg.df, label, pv = 0.05, lfc = 0.25) {
  #  deg.df <- XGJ.NEMOBvsNEMOA.degs
  #  label <- "NEMOBvsNEMOA"
  if (!requireNamespace("MAGeCKFlute", quietly = TRUE)) {
    stop("Package 'MAGeCKFlute' is required to convert gene IDs.")
  }
  deg.df$Entrez <- suppressPackageStartupMessages(
    MAGeCKFlute::TransGeneID(
      rownames(deg.df),
      "Symbol",
      "Entrez",
      organism = "hsa"
    )
  )
  up.entrez <- (deg.df %>% filter(p_val < pv, avg_log2FC > lfc))$Entrez
  dw.entrez <- (deg.df %>% filter(p_val < pv, avg_log2FC < -lfc))$Entrez
  up.entrez <- unique(na.omit(up.entrez))
  dw.entrez <- unique(na.omit(dw.entrez))
  color <- "YlGnBu"
  up.enricher <- enrich_combind(
    up.entrez,
    species = "human",
    reactome = FALSE,
    hallmark = FALSE
  )
  dw.enricher <- enrich_combind(
    dw.entrez,
    species = "human",
    reactome = FALSE,
    hallmark = FALSE
  )
  ################## go_bp
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
        "Up-regulated Genes in GO_BP (",
        label,
        " n= ",
        length(up.entrez),
        ")"
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
        "Down-regulated Genes in GO_BP (",
        label,
        " n= ",
        length(dw.entrez),
        ")"
      )
    )
  ################## kegg
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
        "Up-regulated Genes in KEGG pathway (",
        label,
        " n= ",
        length(up.entrez),
        ")"
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
        "Down-regulated Genes in KEGG pathway (",
        label,
        " n= ",
        length(dw.entrez),
        ")"
      )
    )
  ###################

  p_combind_up <- enricher_plot(up.enricher) +
    ggtitle(
      label = paste0(
        "Up-regulated Genes Enrichment Analysis(",
        label,
        " n= ",
        length(up.entrez),
        ")"
      )
    )
  p_combind_dw <- enricher_plot(dw.enricher) +
    ggtitle(
      label = paste0(
        "Down-regulated Genes Enrichment Analysis(",
        label,
        " n= ",
        length(dw.entrez),
        ")"
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


#' Run configurable gene-set enrichment analysis
#'
#' Builds a ranked Entrez gene vector from common bulk or single-cell DEG
#' layouts and runs GSEA against any supported human or mouse MSigDB collection,
#' as well as the native GO and KEGG functions from `clusterProfiler`.
#'
#' @param deg.df Differential-expression data frame.
#' @param species Species name or code for gene annotation and gene sets.
#' @param rank_by Ranking preset (`"log2FC"`, `"stat"`, `"signed_p"`,
#'   `"signed_padj"`, `"log2FC_p"`, or `"log2FC_padj"`) or the name of a
#'   numeric column.
#' @param dbs Character vector of databases to run. In addition to `"GO_BP"`,
#'   `"GO_CC"`, `"GO_MF"`, `"GO_ALL"`, `"KEGG"`, and `"Reactome"`, this accepts
#'   every database name supported by [EnrichMSigDB()], such as `"H"`, `"C1"`,
#'   `"C2_CP_REACTOME"`, `"C5_GO_BP"`, `"C7"`, `"MH"`, `"M2_CP_REACTOME"`,
#'   or `"M5_GO_BP"`. Use `"MSIGDB_ALL"` to run every MSigDB database available
#'   for the selected species. `"Hallmark"` selects `H` for human and `MH` for
#'   mouse.
#' @param gene_id_col Column containing Entrez identifiers. When it is `Entrez`
#'   and absent, a `Symbol` column is converted automatically.
#' @param log2fc_col Optional log2 fold-change column; common names are detected
#'   when `NULL`.
#' @param stat_col Column containing the DESeq2 or equivalent test statistic.
#' @param p_col,padj_col Optional raw and adjusted P-value columns; common names
#'   are detected when `NULL`.
#' @param collapse_dup Method used for duplicated gene IDs: largest absolute
#'   score, mean score, or first occurrence.
#' @param pvalueCutoff P-value cutoff supplied to GSEA functions.
#' @param minGSSize,maxGSSize Minimum and maximum gene-set sizes.
#' @param p_floor Smallest P value used when constructing logarithmic rankings.
#' @param seed Random seed.
#' @param kegg_source Run KEGG using `clusterProfiler` or MSigDB gene sets.
#' @param readable Logical; convert leading-edge Entrez IDs to symbols where
#'   possible.
#' @param verbose Logical; show enrichment progress.
#'
#' @return A named list of requested GSEA result objects plus `geneList` and
#'   `rank_info` describing the constructed ranking.
#' @export
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
  p_floor = 10^-300,
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
  species <- .abel_normalize_species(species)

  if (species == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package not installed: org.Hs.eg.db")
    }
    msig_species <- "Homo sapiens"
    db_species <- "HS"
    kegg_org <- "hsa"
    current_OrgDb <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package not installed: org.Mm.eg.db")
    }
    msig_species <- "Mus musculus"
    db_species <- "MM"
    kegg_org <- "mmu"
    current_OrgDb <- org.Mm.eg.db::org.Mm.eg.db
  }
  msigdb_map <- .abel_msigdb_database_map(species)

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
    msigdbr_args <- names(formals(msigdbr::msigdbr))
    if (species == "mouse" && !"db_species" %in% msigdbr_args) {
      stop(
        "Native Mouse MSigDB requires msigdbr >= 10.0.0. ",
        "Please update the 'msigdbr' package."
      )
    }

    if ("collection" %in% msigdbr_args) {
      args <- list(
        db_species = db_species,
        species = msig_species,
        collection = collection
      )
      if (!is.null(subcollection)) {
        args$subcollection <- subcollection
      }
    } else {
      args <- list(species = msig_species, category = collection)
      if (!is.null(subcollection)) {
        args$subcategory <- subcollection
      }
    }
    m <- do.call(msigdbr::msigdbr, args)

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

  run_msigdb_gsea <- function(database, label = database) {
    database_info <- msigdb_map[[database]]
    collection <- unname(database_info[1])
    subcollection <- unname(database_info[2])
    if (is.na(subcollection)) {
      subcollection <- NULL
    }

    t2g <- safe_run(paste0(label, " MSigDB loading"), function() {
      get_msig_t2g(collection, subcollection)
    })
    if (is.null(t2g)) {
      return(NULL)
    }
    run_t2g_gsea(label, t2g)
  }

  ## 6. run selected databases
  result <- list()

  if (!is.character(dbs) || !length(dbs) || anyNA(dbs) ||
    any(!nzchar(trimws(dbs)))) {
    stop("dbs must contain at least one non-empty database name.")
  }
  dbs <- unlist(lapply(dbs, function(db) {
    if (toupper(trimws(db)) %in% c("MSIGDB_ALL", "ALL_MSIGDB")) {
      names(msigdb_map)
    } else {
      db
    }
  }), use.names = FALSE)

  for (db in dbs) {
    db <- trimws(db)
    db_key <- toupper(gsub("[ .-]", "_", db))

    if (db_key == "HALLMARK") {
      hallmark_db <- if (species == "human") "H" else "MH"
      result[["Hallmark"]] <- run_msigdb_gsea(hallmark_db, "Hallmark")
    } else if (db_key %in% c("H", "MSIGDB_H")) {
      if (species != "human") {
        warning("H is only available for species = 'human'.", call. = FALSE)
      } else {
        output_name <- if (db_key == "H") "H" else "Hallmark"
        result[[output_name]] <- run_msigdb_gsea("H", output_name)
      }
    } else if (db_key %in% c("MH", "MOUSE_HALLMARK", "MSIGDB_MH")) {
      if (species != "mouse") {
        warning("MH is only available for species = 'mouse'.", call. = FALSE)
      } else {
        output_name <- if (db_key == "MH") "MH" else "Hallmark"
        result[[output_name]] <- run_msigdb_gsea("MH", output_name)
      }
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
        if (species == "human") {
          result[["KEGG_msigdb"]] <- run_msigdb_gsea(
            "C2_CP_KEGG_LEGACY",
            "KEGG_msigdb"
          )
        } else {
          warning(
            "Native Mouse MSigDB does not provide KEGG; use ",
            "kegg_source = 'clusterProfiler'.",
            call. = FALSE
          )
        }
      }
    } else if (db_key %in% c("REACTOME", "MSIGDB_REACTOME")) {
      reactome_db <- if (species == "human") {
        "C2_CP_REACTOME"
      } else {
        "M2_CP_REACTOME"
      }
      result[["Reactome"]] <- run_msigdb_gsea(reactome_db, "Reactome")
    } else if (db_key %in% c("C2_CP_KEGG", "MSIGDB_KEGG")) {
      if (species == "human") {
        result[["KEGG_msigdb"]] <- run_msigdb_gsea(
          "C2_CP_KEGG_LEGACY",
          "KEGG_msigdb"
        )
      } else {
        warning(
          "Native Mouse MSigDB does not provide KEGG collections.",
          call. = FALSE
        )
      }
    } else {
      msigdb_database <- .abel_match_msigdb_database(db, species)
      if (!is.na(msigdb_database)) {
        result[[msigdb_database]] <- run_msigdb_gsea(
          msigdb_database,
          msigdb_database
        )
      } else {
        warning(
          "Unknown database for species '",
          species,
          "': ",
          db,
          ". Supported MSigDB databases: ",
          paste(names(msigdb_map), collapse = ", "),
          ". Additional options: GO_BP, GO_CC, GO_MF, GO_ALL, KEGG, Reactome, MSIGDB_ALL.",
          call. = FALSE
        )
      }
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
