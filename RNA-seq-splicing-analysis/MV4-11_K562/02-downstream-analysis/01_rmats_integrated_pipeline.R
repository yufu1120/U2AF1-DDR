#!/usr/bin/env Rscript

############################################################
# Integrated rMATS downstream pipeline
#
# Features
# - Reads rMATS JC/JCEC outputs for all AS types
# - Preserves original coordinate-defining columns
# - Builds analysis-ready master tables
# - Classifies events into:
#     more_skipped_in_group2
#     more_included_in_group2
#     nonsignificant
#     unassigned
# - Builds/loads cached longest-transcript exon reference from GTF
# - Annotates event_display labels from longest transcript
# - Builds sashimi manifest + one-row event files
# - Generates general summary plots for all AS types
# - Extracts 3' acceptor windows for:
#     SE   -> cassette_acceptor
#     A3SS -> long_acceptor, short_acceptor
# - Writes sequence tables, FASTA inputs, A3SS audit table
# - Generates sequence QC plots, nucleotide-frequency plots, and logos
#
# Notes
# - Sequence extraction requires genome FASTA + .fai index
# - rMATS coordinates are handled using the conventions discussed:
#     starts are 0-based
#     ends are 1-based inclusive
# - Final extracted sequences are oriented to transcript direction
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(readr)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(IRanges)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Biostrings)
  library(Rsamtools)
  library(ggseqlogo)
})

options(stringsAsFactors = FALSE)
setwd(".")

# ==========================================================
# Part 2. Global CONFIG
# ==========================================================
CONFIG <- list(
  # ----- paths -----
  rmats_root = ".",
  output_root = "rmats_integrated_output-HDAC8_OE",#change the folder here
  gtf_file = "./genome.gtf",
  genome_fasta = "./genome.fa.gz",   # REQUIRED for sequence extraction

  # ----- cohorts / comparison labels -----
  cohorts = c(
    "MV4CD531" = "1_MV4CD531",
    "MV4HDAC8" = "2_MV4HDAC8"
    # "WTKD"     = "3_WTKD",
    # "K23RKD"   = "4_K23RKD",
    # "K175RKD"  = "5_K175RKD",
    # "MUTKD"    = "6_MUTKD"
  ),
  group1_label = "NIR",
  group2_label = "IR",

  # ----- rMATS settings -----
  all_as_types = c("SE", "RI", "MXE", "A3SS", "A5SS"),
  motif_as_types = c("SE", "A3SS"),
  cutoff_pvalue = 0.05,
  cutoff_fdr = 0.25,
  cutoff_delta_psi = 0.05,
  min_total_counts = 10,
  use_jcec = TRUE,

  # ----- IncLevelDifference direction safety check -----
  # rMATS defines IncLevelDifference as mean(IncLevel1) - mean(IncLevel2).
  # Downstream reporting uses the explicit opposite convention:
  # delta_psi_group2_minus_group1 = mean(IncLevel2) - mean(IncLevel1).
  inclevel_difference_tolerance = 0.002,
  fail_on_inclevel_direction_mismatch = TRUE,

  # ----- annotation settings -----
  exon_match_tolerance = 3,

  # ----- sashimi / event selection filters -----
  plot_filters = list(
    gene_symbols = c("BRAD1", "BRCA2", "ATM", "MDC1", "UIMC1", "RAD50", "RAD9A", "RAD9B", "WRN", "HELLS", "CHAF1A"),
    AS_type = NULL,
    cohort = NULL,
    filter_counts = TRUE,
    filter_PSI = TRUE,
    filter_PValue = TRUE,
    filter_FDR = FALSE,
    event_class = NULL,
    annotation_status = NULL,
    min_abs_dpsi = NULL,
    max_fdr = NULL
  ),

  # ----- sequence-analysis settings -----
  species = "human",
  genome_build = "hg38",
  chr_name_mode = "auto",   # one of: auto, as_is, add_chr, drop_chr
  upstream_nt = 20,
  downstream_nt = 3,
  a3ss_audit_n_per_cohort = 20,
  keep_primary_chromosomes_only = TRUE,

  # ----- plotting settings -----
  point_alpha = 0.45,
  point_size = 0.9,
  min_label_count = 30
)

# ==========================================================
# Sashimi-only branch CONFIG
# ==========================================================
# This branch loads the cached annotated master table and does NOT rerun
# rMATS import, summary plotting, sequence extraction, or motif analysis.
SASHIMI_CONFIG <- list(
  run_name = "DDR_gene_list_20260713",

  # Use either gene_symbols or gene_file. gene_file may contain one gene per
  # line, or a TSV/CSV whose first column contains gene symbols.
  gene_symbols = c("FANCA", "RAD51D", "USP28", "NUCKS1", "IKBKG", "FAAP20"),
  gene_file = NULL,

  annotated_master_rds = file.path(
    CONFIG$output_root,
    "processed", "annotated_master", "all_cohorts.annotated_master.rds"
  ),
  gtf_reference_rds = file.path(
    CONFIG$output_root,
    "reference", "longest_transcript_exons.rds"
  ),
  runs_root = file.path(CONFIG$output_root, "sashimi_runs"),

  # Event filters for the new gene list.
  AS_type = NULL,
  cohort = NULL,
  filter_counts = TRUE,
  filter_PSI = TRUE,
  filter_PValue = TRUE,
  filter_FDR = FALSE,
  event_class = NULL,
  annotation_status = NULL,
  min_abs_dpsi = NULL,
  max_fdr = NULL,

  overwrite_event_inputs = TRUE
)

{
  
  # ==========================================================
  # Part 3. General utility helpers
  # ==========================================================
  ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
  }
  
  safe_numeric <- function(x) suppressWarnings(as.numeric(x))
  
  mean_from_comma_string <- function(x) {
    if (is.null(x) || length(x) == 0) return(rep(NA_real_, length(x)))
    vapply(x, function(one) {
      if (is.na(one) || one == "") return(NA_real_)
      vals <- strsplit(one, ",", fixed = TRUE)[[1]]
      vals <- suppressWarnings(as.numeric(vals))
      mean(vals, na.rm = TRUE)
    }, numeric(1))
  }
  
  n_items_from_comma_string <- function(x) {
    if (is.null(x) || length(x) == 0) return(rep(NA_integer_, length(x)))
    vapply(x, function(one) {
      if (is.na(one) || one == "") return(NA_integer_)
      length(strsplit(one, ",", fixed = TRUE)[[1]])
    }, integer(1))
  }
  
  sanitize_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
  
  write_table_files <- function(df, out_prefix) {
    readr::write_tsv(df, paste0(out_prefix, ".tsv"))
    saveRDS(df, paste0(out_prefix, ".rds"))
    openxlsx::write.xlsx(df, paste0(out_prefix, ".xlsx"), overwrite = TRUE)
  }
  
  make_class_factor <- function(x) {
    factor(
      x,
      levels = c("more_skipped_in_group2", "more_included_in_group2", "nonsignificant", "unassigned")
    )
  }
  
  make_plot_theme <- function() {
    theme_bw(base_size = 11) +
      theme(
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
      )
  }
  
  save_plot <- function(plot_obj, filename, width = 7, height = 5, dpi = 300) {
    ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, units = "in")
  }
  
  class_colors <- c(
    more_skipped_in_group2 = "#990000",
    more_included_in_group2 = "#0b5394",
    nonsignificant = "grey65",
    unassigned = "grey85"
  )
  
  all_as_order <- c("SE", "RI", "MXE", "A3SS", "A5SS")
  site_type_order <- c("cassette_acceptor", "long_acceptor", "short_acceptor")
  
  # ==========================================================
  # Part 4. rMATS reading + analysis-ready master table
  # ==========================================================
  get_rmats_filename <- function(as_type, use_jcec = TRUE) {
    suffix <- if (isTRUE(use_jcec)) "JCEC" else "JC"
    paste0(as_type, ".MATS.", suffix, ".txt")
  }
  
  get_rmats_file <- function(root, cohort_dir, as_type, use_jcec = TRUE) {
    file.path(root, cohort_dir, "output", "final_output", get_rmats_filename(as_type, use_jcec))
  }
  
  add_feature_label <- function(df, as_type) {
    if (as_type == "SE") {
      df %>% mutate(feature_label = paste0(chr, ":", upstreamEE, ",", exonStart_0base, "-", exonEnd, ",", downstreamES))
    } else if (as_type == "RI") {
      df %>% mutate(feature_label = paste0(chr, ":", upstreamEE, ",", riExonStart_0base, "-", riExonEnd, ",", downstreamES))
    } else if (as_type == "MXE") {
      df %>% mutate(feature_label = paste0(chr, ":", upstreamEE, ",", `1stExonStart_0base`, "-", `1stExonEnd`, ",", `2ndExonStart_0base`, "-", `2ndExonEnd`, ",", downstreamES))
    } else if (as_type %in% c("A3SS", "A5SS")) {
      df %>% mutate(feature_label = paste0(chr, ":", longExonStart_0base, "-", longExonEnd, ",", shortES, "-", shortEE, ",", flankingES, "-", flankingEE))
    } else {
      df %>% mutate(feature_label = NA_character_)
    }
  }
  
  required_coordinate_columns <- list(
    SE   = c("chr", "strand", "upstreamES", "upstreamEE", "exonStart_0base", "exonEnd", "downstreamES", "downstreamEE"),
    RI   = c("chr", "strand", "upstreamES", "upstreamEE", "riExonStart_0base", "riExonEnd", "downstreamES", "downstreamEE"),
    MXE  = c("chr", "strand", "upstreamES", "upstreamEE", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "downstreamES", "downstreamEE"),
    A3SS = c("chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"),
    A5SS = c("chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
  )
  
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
  }

  summarize_inclevel_direction_qc <- function(df, config, out_dir = NULL) {
    tolerance <- config$inclevel_difference_tolerance %||% 0.002

    summary_tbl <- df %>%
      group_by(cohort, AS_type) %>%
      summarise(
        n_rows = n(),
        n_checkable = sum(inclevel_direction_check != "uncheckable"),
        n_match_group1_minus_group2 = sum(inclevel_direction_check == "pass_group1_minus_group2"),
        n_match_group2_minus_group1 = sum(inclevel_direction_check == "fail_matches_group2_minus_group1"),
        n_unresolved = sum(inclevel_direction_check == "fail_unresolved"),
        n_uncheckable = sum(inclevel_direction_check == "uncheckable"),
        max_absolute_error = ifelse(
          any(is.finite(inclevel_difference_abs_error)),
          max(inclevel_difference_abs_error[is.finite(inclevel_difference_abs_error)]),
          NA_real_
        ),
        group1_label = dplyr::first(group1_label),
        group2_label = dplyr::first(group2_label),
        raw_difference_definition = paste0(dplyr::first(group1_label), " minus ", dplyr::first(group2_label)),
        canonical_difference_definition = paste0(dplyr::first(group2_label), " minus ", dplyr::first(group1_label)),
        tolerance = tolerance,
        check_status = ifelse(
          n_match_group2_minus_group1 == 0 & n_unresolved == 0,
          "PASS", "FAIL"
        ),
        .groups = "drop"
      )

    mismatch_tbl <- df %>%
      filter(inclevel_direction_check %in% c(
        "fail_matches_group2_minus_group1", "fail_unresolved"
      )) %>%
      select(
        cohort, AS_type, ID, geneSymbol, feature_label,
        IncLevel1, IncLevel2, mean1, mean2,
        IncLevelDifference_raw,
        rmats_delta_psi_group1_minus_group2,
        delta_psi_group2_minus_group1,
        inclevel_difference_abs_error,
        inclevel_opposite_abs_error,
        inclevel_direction_check
      )

    if (!is.null(out_dir)) {
      ensure_dir(out_dir)
      write_table_files(summary_tbl, file.path(out_dir, "inclevel_difference_check"))
      if (nrow(mismatch_tbl) > 0) {
        write_table_files(mismatch_tbl, file.path(out_dir, "inclevel_difference_mismatches"))
      }
    }

    print(summary_tbl, n = Inf)

    if (isTRUE(config$fail_on_inclevel_direction_mismatch) &&
        any(summary_tbl$check_status == "FAIL")) {
      stop(
        "IncLevelDifference direction QC failed. rMATS is expected to use ",
        "mean(IncLevel1) - mean(IncLevel2). Review: ",
        file.path(out_dir %||% "<QC directory>", "inclevel_difference_check.tsv")
      )
    }

    invisible(list(summary = summary_tbl, mismatches = mismatch_tbl))
  }

  upgrade_delta_psi_columns <- function(df, config) {
    required <- c("mean1", "mean2", "IncLevelDifference")
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
      stop("Cached master is missing required columns: ", paste(missing, collapse = ", "))
    }

    if (!"IncLevelDifference_raw" %in% names(df)) {
      df$IncLevelDifference_raw <- safe_numeric(df$IncLevelDifference)
    }
    df$rmats_delta_psi_group1_minus_group2 <- df$mean1 - df$mean2
    df$delta_psi_group2_minus_group1 <- df$mean2 - df$mean1
    df$inclevel_difference_abs_error <- abs(
      df$IncLevelDifference_raw - df$rmats_delta_psi_group1_minus_group2
    )
    df$inclevel_opposite_abs_error <- abs(
      df$IncLevelDifference_raw - df$delta_psi_group2_minus_group1
    )
    tolerance <- config$inclevel_difference_tolerance %||% 0.002
    df$inclevel_direction_check <- case_when(
      !is.finite(df$mean1) | !is.finite(df$mean2) | !is.finite(df$IncLevelDifference_raw) ~ "uncheckable",
      df$inclevel_difference_abs_error <= tolerance ~ "pass_group1_minus_group2",
      df$inclevel_opposite_abs_error <= tolerance ~ "fail_matches_group2_minus_group1",
      TRUE ~ "fail_unresolved"
    )
    df
  }

  read_one_rmats <- function(file, cohort_name, as_type, config) {
    if (!file.exists(file)) stop("rMATS file not found: ", file)

    df <- read.delim(file, sep = "\t", header = TRUE, check.names = FALSE)

    original_names <- names(df)
    if (anyDuplicated(original_names)) {
      dup_names <- unique(original_names[duplicated(original_names)])
      message("[read_one_rmats] Duplicate columns detected in ", basename(file), ": ",
              paste(dup_names, collapse = ", "),
              ". Renaming with make.unique() so the first copy keeps its original name.")
      names(df) <- make.unique(original_names, sep = "__dup")
    }

    df <- as_tibble(df)
    df <- add_feature_label(df, as_type)

    coord_cols <- required_coordinate_columns[[as_type]]
    missing_coord_cols <- setdiff(coord_cols, colnames(df))
    if (length(missing_coord_cols) > 0) stop(
      "Missing required coordinate columns in ", as_type, ": ",
      paste(missing_coord_cols, collapse = ", ")
    )

    required_stat_cols <- c(
      "ID", "geneSymbol", "IncLevel1", "IncLevel2",
      "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
      "PValue", "FDR", "IncLevelDifference"
    )
    missing_stat_cols <- setdiff(required_stat_cols, colnames(df))
    if (length(missing_stat_cols) > 0) stop(
      "Missing expected columns in ", basename(file), ": ",
      paste(missing_stat_cols, collapse = ", ")
    )

    tolerance <- config$inclevel_difference_tolerance %||% 0.002

    df %>% mutate(
      cohort = cohort_name,
      AS_type = as_type,
      feature_label_short = paste0(AS_type, "-", ID, "-", geneSymbol),
      mean1 = mean_from_comma_string(IncLevel1),
      mean2 = mean_from_comma_string(IncLevel2),
      mean_IJC_SAMPLE_1 = mean_from_comma_string(IJC_SAMPLE_1),
      mean_SJC_SAMPLE_1 = mean_from_comma_string(SJC_SAMPLE_1),
      mean_IJC_SAMPLE_2 = mean_from_comma_string(IJC_SAMPLE_2),
      mean_SJC_SAMPLE_2 = mean_from_comma_string(SJC_SAMPLE_2),
      nrep_1 = n_items_from_comma_string(IncLevel1),
      nrep_2 = n_items_from_comma_string(IncLevel2),
      PValue = safe_numeric(PValue),
      FDR = safe_numeric(FDR),
      IncLevelDifference_raw = safe_numeric(IncLevelDifference),
      rmats_delta_psi_group1_minus_group2 = mean1 - mean2,
      delta_psi_group2_minus_group1 = mean2 - mean1,
      inclevel_difference_abs_error = abs(
        IncLevelDifference_raw - rmats_delta_psi_group1_minus_group2
      ),
      inclevel_opposite_abs_error = abs(
        IncLevelDifference_raw - delta_psi_group2_minus_group1
      ),
      inclevel_direction_check = case_when(
        !is.finite(mean1) | !is.finite(mean2) | !is.finite(IncLevelDifference_raw) ~ "uncheckable",
        inclevel_difference_abs_error <= tolerance ~ "pass_group1_minus_group2",
        inclevel_opposite_abs_error <= tolerance ~ "fail_matches_group2_minus_group1",
        TRUE ~ "fail_unresolved"
      ),
      # Keep the historical column name for compatibility, but never use it
      # without the explicit definition columns above.
      IncLevelDifference = IncLevelDifference_raw
    )
  }
  
  build_analysis_table_one_cohort <- function(cohort_name, cohort_dir, config) {
    dfs <- lapply(config$all_as_types, function(as_type) {
      file <- get_rmats_file(config$rmats_root, cohort_dir, as_type, config$use_jcec)
      read_one_rmats(file, cohort_name = cohort_name, as_type = as_type, config = config)
    })
    
    bind_rows(dfs) %>%
      mutate(
        cohort = as.character(cohort),
        AS_type = factor(AS_type, levels = all_as_order),
        group1_label = config$group1_label,
        group2_label = config$group2_label,
        filter_counts = case_when(
          (mean_IJC_SAMPLE_1 * nrep_1 >= config$min_total_counts & mean_SJC_SAMPLE_1 * nrep_1 >= config$min_total_counts) ~ TRUE,
          (mean_IJC_SAMPLE_2 * nrep_2 >= config$min_total_counts & mean_SJC_SAMPLE_2 * nrep_2 >= config$min_total_counts) ~ TRUE,
          TRUE ~ FALSE
        ),
        filter_PSI = abs(delta_psi_group2_minus_group1) >= config$cutoff_delta_psi,
        filter_PValue = PValue <= config$cutoff_pvalue,
        filter_FDR = FDR <= config$cutoff_fdr
      )
  }
  
  build_analysis_tables_all <- function(config) {
    out_dir <- ensure_dir(file.path(config$output_root, "processed", "analysis_ready"))
    qc_dir <- ensure_dir(file.path(config$output_root, "processed", "qc"))

    all_tables <- imap(config$cohorts, function(cohort_dir, cohort_name) {
      build_analysis_table_one_cohort(cohort_name, cohort_dir, config)
    })
    all_df <- bind_rows(all_tables)

    # Mandatory checkpoint before tables are classified or used downstream.
    summarize_inclevel_direction_qc(all_df, config, out_dir = qc_dir)

    walk(names(config$cohorts), function(cohort_name) {
      one <- all_df %>% filter(cohort == cohort_name)
      prefix <- file.path(out_dir, sanitize_filename(paste0(cohort_name, ".analysis_ready")))
      write_table_files(one, prefix)
    })

    all_df
  }
  
  classify_events <- function(df, config) {
    df %>%
      mutate(
        is_significant_for_motif = filter_counts & filter_PValue & filter_PSI,
        is_background_for_motif = filter_counts & (FDR > config$cutoff_fdr) & (abs(delta_psi_group2_minus_group1) < config$cutoff_delta_psi),
        event_class = case_when(
          is_significant_for_motif & delta_psi_group2_minus_group1 < 0 ~ "more_skipped_in_group2",
          is_significant_for_motif & delta_psi_group2_minus_group1 > 0 ~ "more_included_in_group2",
          is_background_for_motif ~ "nonsignificant",
          TRUE ~ "unassigned"
        ),
        event_class = make_class_factor(event_class)
      )
  }
  
  # ==========================================================
  # Part 5. GTF reference + longest-transcript annotation
  # ==========================================================
  build_or_load_gtf_reference <- function(gtf_file, output_root, force_rebuild = FALSE) {
    if (!file.exists(gtf_file)) {
      stop("GTF not found: ", gtf_file)
    }
    
    ref_dir <- ensure_dir(file.path(output_root, "reference"))
    ref_rds <- file.path(ref_dir, "longest_transcript_exons.rds")
    
    if (file.exists(ref_rds) && !isTRUE(force_rebuild)) {
      message("[GTF reference] Loading cached reference: ", ref_rds)
      ref_obj <- readRDS(ref_rds)
      return(ref_obj)
    }
    
    message("[GTF reference] Building reference from GTF: ", gtf_file)
    
    flatten_gtf_col <- function(x, n) {
      if (is.null(x)) return(rep(NA_character_, n))
      
      if (isS4(x)) {
        x <- as.vector(x)
      }
      
      if (is.list(x)) {
        x <- vapply(x, function(one) {
          if (length(one) == 0 || all(is.na(one))) return(NA_character_)
          as.character(one[[1]])
        }, character(1))
      } else {
        x <- as.character(x)
      }
      
      if (length(x) != n) {
        if (length(x) == 1) {
          x <- rep(x, n)
        } else {
          stop("GTF column length mismatch during flattening.")
        }
      }
      
      x
    }
    
    normalize_gtf_df <- function(gtf) {
      mcols_df <- as.data.frame(mcols(gtf), optional = TRUE, stringsAsFactors = FALSE)
      n <- length(gtf)
      
      tibble(
        seqnames = as.character(as.vector(seqnames(gtf))),
        start = as.integer(start(gtf)),
        end = as.integer(end(gtf)),
        strand = as.character(as.vector(strand(gtf))),
        type = if ("type" %in% names(mcols_df)) flatten_gtf_col(mcols_df$type, n) else rep(NA_character_, n),
        gene_id = if ("gene_id" %in% names(mcols_df)) flatten_gtf_col(mcols_df$gene_id, n) else rep(NA_character_, n),
        gene_name = if ("gene_name" %in% names(mcols_df)) flatten_gtf_col(mcols_df$gene_name, n) else rep(NA_character_, n),
        transcript_id = if ("transcript_id" %in% names(mcols_df)) flatten_gtf_col(mcols_df$transcript_id, n) else rep(NA_character_, n),
        transcript_name = if ("transcript_name" %in% names(mcols_df)) flatten_gtf_col(mcols_df$transcript_name, n) else rep(NA_character_, n)
      ) %>%
        dplyr::mutate(
          seqnames = as.character(seqnames),
          strand = as.character(strand),
          type = as.character(type),
          gene_id = as.character(gene_id),
          gene_name = as.character(gene_name),
          transcript_id = as.character(transcript_id),
          transcript_name = as.character(transcript_name),
          start = as.integer(start),
          end = as.integer(end)
        )
    }
    
    load_longest_transcript_exons <- function(gtf_file) {
      gtf <- rtracklayer::import(gtf_file)
      gtf_df <- normalize_gtf_df(gtf)
      
      exon_df <- gtf_df %>%
        dplyr::filter(type == "exon", !is.na(transcript_id), !is.na(gene_name)) %>%
        dplyr::distinct(seqnames, start, end, strand, gene_id, gene_name, transcript_id, transcript_name)
      
      if (nrow(exon_df) == 0) {
        stop("No exon records found in GTF: ", gtf_file)
      }
      
      tx_lengths <- exon_df %>%
        dplyr::group_by(gene_name, transcript_id, transcript_name, gene_id, strand, seqnames) %>%
        dplyr::summarise(tx_len = sum(end - start + 1L), .groups = "drop") %>%
        dplyr::arrange(gene_name, dplyr::desc(tx_len), transcript_id)
      
      longest_tx <- tx_lengths %>%
        dplyr::group_by(gene_name) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      
      longest_exons <- exon_df %>%
        dplyr::inner_join(
          longest_tx %>% dplyr::select(gene_name, transcript_id, tx_len),
          by = c("gene_name", "transcript_id")
        ) %>%
        dplyr::group_by(gene_name, transcript_id) %>%
        dplyr::arrange(start, end, .by_group = TRUE) %>%
        dplyr::mutate(
          exon_number_plus = dplyr::row_number(),
          exon_number_minus = rev(dplyr::row_number()),
          exon_number = dplyr::case_when(
            strand == "+" ~ exon_number_plus,
            strand == "-" ~ exon_number_minus,
            TRUE ~ exon_number_plus
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(gene_name, gene_id, transcript_id, transcript_name, tx_len, seqnames, strand, start, end, exon_number)
      
      longest_exons
    }
    
    longest_exons <- load_longest_transcript_exons(gtf_file)
    
    ref_obj <- list(
      gtf_file = normalizePath(gtf_file, winslash = "/", mustWork = FALSE),
      build_time = as.character(Sys.time()),
      longest_exons = longest_exons
    )
    
    saveRDS(ref_obj, ref_rds)
    message("[GTF reference] Saved reference to: ", ref_rds)
    
    ref_obj
  }
  
  closest_exon_match <- function(exons_df, chr, strand, target_start, target_end, tol = 3) {
    exons_sub <- exons_df %>% dplyr::filter(seqnames == chr, strand == strand)
    
    if (nrow(exons_sub) == 0) {
      return(list(exon_number = NA_integer_, match_type = "no_match", annotation_status = "no_match", overlap_bp = NA_integer_))
    }
    
    exons_sub <- exons_sub %>%
      dplyr::mutate(
        exact = (start == target_start & end == target_end),
        near_boundary = (abs(start - target_start) <= tol & abs(end - target_end) <= tol),
        overlap_bp = pmax(0L, pmin(end, target_end) - pmax(start, target_start) + 1L),
        width = end - start + 1L,
        target_width = target_end - target_start + 1L
      )
    
    exact_hits <- exons_sub %>% dplyr::filter(exact)
    if (nrow(exact_hits) >= 1) {
      chosen <- exact_hits %>% dplyr::arrange(exon_number) %>% dplyr::slice(1)
      return(list(
        exon_number = chosen$exon_number[[1]],
        match_type = "exact",
        annotation_status = ifelse(nrow(exact_hits) == 1, "matched", "ambiguous"),
        overlap_bp = chosen$overlap_bp[[1]]
      ))
    }
    
    near_hits <- exons_sub %>% dplyr::filter(near_boundary)
    if (nrow(near_hits) >= 1) {
      chosen <- near_hits %>%
        dplyr::mutate(dist_sum = abs(start - target_start) + abs(end - target_end)) %>%
        dplyr::arrange(dist_sum, exon_number) %>%
        dplyr::slice(1)
      return(list(
        exon_number = chosen$exon_number[[1]],
        match_type = "tolerant_boundary",
        annotation_status = ifelse(nrow(near_hits) == 1, "matched", "ambiguous"),
        overlap_bp = chosen$overlap_bp[[1]]
      ))
    }
    
    overlap_hits <- exons_sub %>% dplyr::filter(overlap_bp > 0)
    if (nrow(overlap_hits) >= 1) {
      max_overlap <- max(overlap_hits$overlap_bp, na.rm = TRUE)
      n_best <- overlap_hits %>% dplyr::filter(overlap_bp == max_overlap) %>% nrow()
      chosen <- overlap_hits %>% dplyr::arrange(dplyr::desc(overlap_bp), exon_number) %>% dplyr::slice(1)
      return(list(
        exon_number = chosen$exon_number[[1]],
        match_type = "max_overlap",
        annotation_status = ifelse(n_best == 1, "matched", "ambiguous"),
        overlap_bp = chosen$overlap_bp[[1]]
      ))
    }
    
    list(exon_number = NA_integer_, match_type = "no_match", annotation_status = "no_match", overlap_bp = 0L)
  }
  
  get_flanking_exon_numbers_for_introns <- function(exons_df, chr, strand, intron_start, intron_end, tol = 3) {
    exons_sub <- exons_df %>% dplyr::filter(seqnames == chr, strand == strand)
    
    if (nrow(exons_sub) < 2) {
      return(list(left_exon = NA_integer_, right_exon = NA_integer_, match_type = "no_match", annotation_status = "no_match"))
    }
    
    intron_df <- exons_sub %>%
      dplyr::arrange(start, end) %>%
      dplyr::mutate(
        next_start = dplyr::lead(start),
        next_end = dplyr::lead(end),
        next_exon_number = dplyr::lead(exon_number),
        intron_start_tx = end + 1L,
        intron_end_tx = next_start - 1L
      ) %>%
      dplyr::filter(!is.na(next_start), intron_start_tx <= intron_end_tx) %>%
      dplyr::mutate(
        exact = (intron_start_tx == intron_start & intron_end_tx == intron_end),
        near_boundary = (abs(intron_start_tx - intron_start) <= tol & abs(intron_end_tx - intron_end) <= tol),
        overlap_bp = pmax(0L, pmin(intron_end_tx, intron_end) - pmax(intron_start_tx, intron_start) + 1L)
      )
    
    if (nrow(intron_df) == 0) {
      return(list(left_exon = NA_integer_, right_exon = NA_integer_, match_type = "no_match", annotation_status = "no_match"))
    }
    
    exact_hits <- intron_df %>% dplyr::filter(exact)
    if (nrow(exact_hits) >= 1) {
      chosen <- exact_hits %>% dplyr::slice(1)
      return(list(
        left_exon = chosen$exon_number[[1]],
        right_exon = chosen$next_exon_number[[1]],
        match_type = "exact",
        annotation_status = ifelse(nrow(exact_hits) == 1, "matched", "ambiguous")
      ))
    }
    
    near_hits <- intron_df %>% dplyr::filter(near_boundary)
    if (nrow(near_hits) >= 1) {
      chosen <- near_hits %>%
        dplyr::mutate(dist_sum = abs(intron_start_tx - intron_start) + abs(intron_end_tx - intron_end)) %>%
        dplyr::arrange(dist_sum, exon_number) %>%
        dplyr::slice(1)
      return(list(
        left_exon = chosen$exon_number[[1]],
        right_exon = chosen$next_exon_number[[1]],
        match_type = "tolerant_boundary",
        annotation_status = ifelse(nrow(near_hits) == 1, "matched", "ambiguous")
      ))
    }
    
    overlap_hits <- intron_df %>% dplyr::filter(overlap_bp > 0)
    if (nrow(overlap_hits) >= 1) {
      max_overlap <- max(overlap_hits$overlap_bp, na.rm = TRUE)
      n_best <- overlap_hits %>% dplyr::filter(overlap_bp == max_overlap) %>% nrow()
      chosen <- overlap_hits %>% dplyr::arrange(dplyr::desc(overlap_bp), exon_number) %>% dplyr::slice(1)
      return(list(
        left_exon = chosen$exon_number[[1]],
        right_exon = chosen$next_exon_number[[1]],
        match_type = "max_overlap",
        annotation_status = ifelse(n_best == 1, "matched", "ambiguous")
      ))
    }
    
    list(left_exon = NA_integer_, right_exon = NA_integer_, match_type = "no_match", annotation_status = "no_match")
  }
  
  normalize_chr <- function(chr) {
    chr <- as.character(chr)
    chr <- gsub("^chr", "", chr)
    chr
  }
  
  annotate_one_event <- function(row, exon_models, tol = 3) {
    gene <- as.character(row$geneSymbol)
    chr_raw <- as.character(row$chr)
    strand <- as.character(row$strand)
    as_type <- as.character(row$AS_type)
    
    chr <- normalize_chr(chr_raw)
    
    gene_exons <- exon_models %>%
      mutate(seqnames = normalize_chr(seqnames)) %>%
      filter(gene_name == gene, seqnames == chr)
    
    if (nrow(gene_exons) == 0) {
      return(tibble(
        transcript_id_used = NA_character_,
        transcript_name_used = NA_character_,
        annotation_basis = "longest_transcript",
        annotation_status = "no_gene_model",
        match_type = "no_match",
        event_display = paste0(gene, "-unannotated_", as_type)
      ))
    }
    
    tx_id <- gene_exons$transcript_id[[1]]
    tx_name <- gene_exons$transcript_name[[1]]
    
    if (as_type == "SE") {
      hit <- closest_exon_match(gene_exons, chr, strand, as.integer(row$exonStart_0base) + 1L, as.integer(row$exonEnd), tol)
      event_display <- if (is.na(hit$exon_number)) paste0(gene, "-unannotated_SE") else paste0(gene, "-exon", hit$exon_number)
      return(tibble(
        transcript_id_used = tx_id,
        transcript_name_used = tx_name,
        annotation_basis = "longest_transcript",
        annotation_status = hit$annotation_status,
        match_type = hit$match_type,
        event_display = event_display
      ))
    }
    
    if (as_type == "RI") {
      hit <- get_flanking_exon_numbers_for_introns(gene_exons, chr, strand, as.integer(row$riExonStart_0base) + 1L, as.integer(row$riExonEnd), tol)
      event_display <- if (is.na(hit$left_exon) || is.na(hit$right_exon)) paste0(gene, "-unannotated_RI") else paste0(gene, "-intron", hit$left_exon, "_", hit$right_exon)
      return(tibble(
        transcript_id_used = tx_id,
        transcript_name_used = tx_name,
        annotation_basis = "longest_transcript",
        annotation_status = hit$annotation_status,
        match_type = hit$match_type,
        event_display = event_display
      ))
    }
    
    if (as_type == "MXE") {
      hit1 <- closest_exon_match(gene_exons, chr, strand, as.integer(row$`1stExonStart_0base`) + 1L, as.integer(row$`1stExonEnd`), tol)
      hit2 <- closest_exon_match(gene_exons, chr, strand, as.integer(row$`2ndExonStart_0base`) + 1L, as.integer(row$`2ndExonEnd`), tol)
      status <- ifelse(any(c(hit1$annotation_status, hit2$annotation_status) == "no_match"),
                       "partial_or_no_match",
                       ifelse(any(c(hit1$annotation_status, hit2$annotation_status) == "ambiguous"), "ambiguous", "matched"))
      event_display <- if (is.na(hit1$exon_number) || is.na(hit2$exon_number)) paste0(gene, "-unannotated_MXE") else paste0(gene, "-exon", hit1$exon_number, "_exon", hit2$exon_number)
      return(tibble(
        transcript_id_used = tx_id,
        transcript_name_used = tx_name,
        annotation_basis = "longest_transcript",
        annotation_status = status,
        match_type = paste(hit1$match_type, hit2$match_type, sep = ";"),
        event_display = event_display
      ))
    }
    
    if (as_type %in% c("A3SS", "A5SS")) {
      hit <- closest_exon_match(gene_exons, chr, strand, as.integer(row$longExonStart_0base) + 1L, as.integer(row$longExonEnd), tol)
      event_display <- if (is.na(hit$exon_number)) paste0(gene, "-unannotated_", as_type) else paste0(gene, "-exon", hit$exon_number, "-", as_type)
      return(tibble(
        transcript_id_used = tx_id,
        transcript_name_used = tx_name,
        annotation_basis = "longest_transcript",
        annotation_status = hit$annotation_status,
        match_type = hit$match_type,
        event_display = event_display
      ))
    }
    
    tibble(
      transcript_id_used = tx_id,
      transcript_name_used = tx_name,
      annotation_basis = "longest_transcript",
      annotation_status = "unsupported_AS_type",
      match_type = "no_match",
      event_display = paste0(gene, "-unannotated_", as_type)
    )
  }
  
  annotate_event_display <- function(
    df,
    gtf_ref,
    tol = 3,
    progress_step = 1000,
    annotate_only_significant = TRUE,
    gene_symbols = NULL
  ) {
    exon_models <- gtf_ref$longest_exons
    
    df_out <- df %>%
      mutate(
        transcript_id_used = NA_character_,
        transcript_name_used = NA_character_,
        annotation_basis = "not_aiming_annotated",
        annotation_status = "not_aiming_annotated",
        match_type = "not_aiming_annotated",
        event_display = feature_label_short
      )
    
    idx <- seq_len(nrow(df_out))
    
    if (isTRUE(annotate_only_significant)) {
      idx <- idx[df_out$is_significant_for_motif[idx] == TRUE]
    }
    
    if (!is.null(gene_symbols)) {
      idx <- idx[df_out$geneSymbol[idx] %in% gene_symbols]
    }
    
    n <- length(idx)
    
    if (n == 0) {
      message("[annotate_event_display] No rows selected for annotation.")
      return(df_out)
    }
    
    message("[annotate_event_display] Annotating ", n, " selected rows out of ", nrow(df_out), " total rows.")
    
    anno_list <- vector("list", n)
    
    for (j in seq_along(idx)) {
      i <- idx[[j]]
      
      if (j %% progress_step == 0 || j == 1 || j == n) {
        message(sprintf(
          "[annotate_event_display] %d / %d selected rows (%.2f%%)",
          j, n, 100 * j / n
        ))
      }
      
      anno_list[[j]] <- annotate_one_event(
        df_out[i, ],
        exon_models,
        tol = tol
      )
    }
    
    anno_tbl <- bind_rows(anno_list)
    
    anno_cols <- c(
      "transcript_id_used",
      "transcript_name_used",
      "annotation_basis",
      "annotation_status",
      "match_type",
      "event_display"
    )
    
    df_out[idx, anno_cols] <- anno_tbl[, anno_cols]
    
    df_out
  }
  
  # ==========================================================
  # Part 6. General filtering + sashimi manifest branch
  # ==========================================================
  filter_events_general <- function(df, gene_symbols = NULL, AS_type = NULL, cohort = NULL,
                                    filter_counts = NULL, filter_PSI = NULL, filter_PValue = NULL,
                                    filter_FDR = NULL, event_class = NULL, annotation_status = NULL,
                                    min_abs_dpsi = NULL, max_fdr = NULL) {
    out <- df
    if (!is.null(gene_symbols)) out <- out %>% filter(geneSymbol %in% gene_symbols)
    if (!is.null(AS_type)) out <- out %>% filter(as.character(.data$AS_type) %in% AS_type)
    if (!is.null(cohort)) out <- out %>% filter(.data$cohort %in% cohort)
    if (!is.null(filter_counts)) out <- out %>% filter(.data$filter_counts == filter_counts)
    if (!is.null(filter_PSI)) out <- out %>% filter(.data$filter_PSI == filter_PSI)
    if (!is.null(filter_PValue)) out <- out %>% filter(.data$filter_PValue == filter_PValue)
    if (!is.null(filter_FDR)) out <- out %>% filter(.data$filter_FDR == filter_FDR)
    if (!is.null(event_class)) out <- out %>% filter(as.character(.data$event_class) %in% event_class)
    if (!is.null(annotation_status)) out <- out %>% filter(.data$annotation_status %in% annotation_status)
    if (!is.null(min_abs_dpsi)) out <- out %>% filter(abs(.data$delta_psi_group2_minus_group1) >= min_abs_dpsi)
    if (!is.null(max_fdr)) out <- out %>% filter(.data$FDR <= max_fdr)
    out
  }
  
  add_plot_ids <- function(df) {
    df %>% arrange(cohort, AS_type, geneSymbol, event_display, feature_label) %>%
      mutate(plot_id = sprintf("plot_%05d", seq_len(n())))
  }
  
  get_rmats_cols_for_type <- function(df, as_type) {
    id_dup_col <- grep("^ID__dup", colnames(df), value = TRUE)
    
    common_tail <- c(
      id_dup_col,
      "IJC_SAMPLE_1", "SJC_SAMPLE_1",
      "IJC_SAMPLE_2", "SJC_SAMPLE_2",
      "IncFormLen", "SkipFormLen",
      "PValue", "FDR",
      "IncLevel1", "IncLevel2",
      "IncLevelDifference"
    )
    
    cols <- switch(
      as_type,
      
      SE = c(
        "ID", "GeneID", "geneSymbol", "chr", "strand",
        "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
        common_tail
      ),
      
      RI = c(
        "ID", "GeneID", "geneSymbol", "chr", "strand",
        "riExonStart_0base", "riExonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
        common_tail
      ),
      
      MXE = c(
        "ID", "GeneID", "geneSymbol", "chr", "strand",
        "1stExonStart_0base", "1stExonEnd",
        "2ndExonStart_0base", "2ndExonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
        common_tail
      ),
      
      A3SS = c(
        "ID", "GeneID", "geneSymbol", "chr", "strand",
        "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE",
        "flankingES", "flankingEE",
        common_tail
      ),
      
      A5SS = c(
        "ID", "GeneID", "geneSymbol", "chr", "strand",
        "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE",
        "flankingES", "flankingEE",
        common_tail
      ),
      
      stop("Unsupported AS type: ", as_type)
    )
    
    cols <- unique(cols)
    
    missing <- setdiff(cols, colnames(df))
    if (length(missing) > 0) {
      stop(
        "Missing required rMATS columns for ", as_type, ": ",
        paste(missing, collapse = ", ")
      )
    }
    
    cols
  }
  
  write_one_row_rmats_inputs <- function(manifest_df, out_dir) {
    input_dir <- ensure_dir(file.path(out_dir, "event_inputs"))
    out_paths <- vector("character", nrow(manifest_df))
    
    for (i in seq_len(nrow(manifest_df))) {
      row <- manifest_df[i, ]
      as_type <- as.character(row$AS_type)
      
      keep <- get_rmats_cols_for_type(manifest_df, as_type)
      one <- row[, keep, drop = FALSE]
      
      # Restore original rMATS duplicate ID header expected by rmats2sashimiplot
      names(one) <- sub("^ID__dup[0-9]*$", "ID", names(one))
      
      file_path <- file.path(
        input_dir,
        as_type,
        paste0(row$plot_id, ".", as_type, ".txt")
      )
      
      ensure_dir(dirname(file_path))
      
      readr::write_tsv(
        one,
        file_path,
        na = "NA",
        quote = "none"
      )
      
      out_paths[i] <- file_path
    }
    
    manifest_df$rmats_event_file <- out_paths
    manifest_df
  }
  
  fmt_num <- function(x, digits = 3) {
    ifelse(is.na(x), "NA", formatC(x, format = "e", digits = digits))
  }
  
  build_sashimi_manifest <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sashimi_manifest"))
    
    filt <- do.call(filter_events_general, c(list(df = df), config$plot_filters))
    filt <- add_plot_ids(filt)
    filt <- write_one_row_rmats_inputs(filt, out_dir)
    
    filt <- filt %>% mutate(
      plot_subdir = file.path(cohort, as.character(AS_type), plot_id),
      expected_pdf_dir = file.path("output", plot_subdir),
      ppt_title = event_display,
      ppt_subtitle = paste0(
        "cohort=", cohort,
        " | AS_type=", AS_type,
        " | feature_label=", feature_label,
        " | PValue=", fmt_num(PValue, 2),
        " | FDR=", fmt_num(FDR, 2),
        " | dPSI_", group2_label, "_minus_", group1_label, "=",
        signif(delta_psi_group2_minus_group1, 3)
      )
    )
    
    write_table_files(filt, file.path(out_dir, "sashimi_manifest"))
    
    filt_sig_gene <- filt %>%
      filter(
        is_significant_for_motif == TRUE,
        geneSymbol %in% config$plot_filters$gene_symbols
      )
    
    readr::write_tsv(
      filt_sig_gene,
      file.path(out_dir, "sashimi_manifest_filter.tsv")
    )
    
    openxlsx::write.xlsx(
      filt_sig_gene,
      file.path(out_dir, "sashimi_manifest_filter.xlsx"),
      overwrite = TRUE
    )
    
    saveRDS(
      filt_sig_gene,
      file.path(out_dir, "sashimi_manifest_filter.rds")
    )
    
    filt
  }

  
  # ==========================================================
  # Part 7. Sequence-analysis module
  # ==========================================================
  # ---------- 7A. Genome FASTA helpers ----------
  open_fasta <- function(fasta_path) {
    if (!file.exists(fasta_path)) stop("Genome FASTA not found: ", fasta_path)
    fai_path <- paste0(fasta_path, ".fai")
    if (!file.exists(fai_path)) {
      stop("FASTA index not found: ", fai_path,
           "\nCreate it first, e.g. Rsamtools::indexFa('", fasta_path, "')")
    }
    fa <- FaFile(fasta_path)
    open(fa)
    fa
  }
  
  close_fasta_safe <- function(fa) {
    try(close(fa), silent = TRUE)
  }
  
  resolve_chr_name <- function(chr, fasta_seqnames, mode = "auto") {
    if (is.na(chr) || is.null(chr) || length(chr) != 1) return(NA_character_)
    chr <- as.character(chr)
    
    if (mode == "as_is") return(if (chr %in% fasta_seqnames) chr else NA_character_)
    
    if (mode == "add_chr") {
      cand <- if (startsWith(chr, "chr")) chr else paste0("chr", chr)
      return(if (cand %in% fasta_seqnames) cand else NA_character_)
    }
    
    if (mode == "drop_chr") {
      cand <- sub("^chr", "", chr)
      return(if (cand %in% fasta_seqnames) cand else NA_character_)
    }
    
    candidates <- unique(c(
      chr,
      if (!startsWith(chr, "chr")) paste0("chr", chr) else sub("^chr", "", chr),
      if (chr %in% c("MT", "M")) c("chrM", "MT", "M") else NULL,
      if (chr == "chrM") c("chrM", "MT", "M") else NULL
    ))
    hit <- candidates[candidates %in% fasta_seqnames]
    if (length(hit) > 0) return(hit[[1]])
    NA_character_
  }
  
  get_seq_safe <- function(fa, chr, start, end, chr_name_mode = "auto") {
    if (is.na(chr) || is.na(start) || is.na(end)) {
      return(list(
        ok = FALSE,
        seq = NA_character_,
        reason = "missing_coordinate",
        chr_used = NA_character_
      ))
    }
    
    if (start < 1 || end < 1 || end < start) {
      return(list(
        ok = FALSE,
        seq = NA_character_,
        reason = "invalid_window",
        chr_used = NA_character_
      ))
    }
    
    si <- seqinfo(fa)
    fasta_seqnames <- names(si)
    
    chr_used <- resolve_chr_name(chr, fasta_seqnames, mode = chr_name_mode)
    
    if (is.na(chr_used)) {
      return(list(
        ok = FALSE,
        seq = NA_character_,
        reason = "chr_not_in_fasta",
        chr_used = NA_character_
      ))
    }
    
    chr_len <- seqlengths(si)[chr_used]
    
    if (is.na(chr_len)) {
      return(list(
        ok = FALSE,
        seq = NA_character_,
        reason = "missing_chr_length",
        chr_used = chr_used
      ))
    }
    
    if (end > chr_len) {
      return(list(
        ok = FALSE,
        seq = NA_character_,
        reason = "window_out_of_bounds",
        chr_used = chr_used
      ))
    }
    
    gr <- GenomicRanges::GRanges(
      seqnames = chr_used,
      ranges = IRanges::IRanges(start = start, end = end)
    )
    
    seq <- tryCatch(
      {
        as.character(Biostrings::getSeq(fa, gr))
      },
      error = function(e) {
        NA_character_
      }
    )
    
    if (is.na(seq)) {
      return(list(
        ok = FALSE,
        seq = NA_character_,
        reason = "getSeq_error_or_invalid_DNA_letters",
        chr_used = chr_used
      ))
    }
    
    list(
      ok = TRUE,
      seq = seq,
      reason = NA_character_,
      chr_used = chr_used
    )
  }
  
  orient_acceptor_window <- function(raw_seq, strand) {
    if (is.na(raw_seq)) return(NA_character_)
    if (strand == "+") return(raw_seq)
    if (strand == "-") return(as.character(reverseComplement(DNAString(raw_seq))))
    NA_character_
  }
  
  # ---------- 7B. 3' site coordinate/window helpers ----------
  compute_window <- function(acceptor_anchor, strand, upstream_nt, downstream_nt) {
    if (is.na(acceptor_anchor) || is.na(strand)) {
      return(list(start = NA_integer_, end = NA_integer_, acceptor_boundary_genomic = NA_integer_, rc = NA))
    }
    
    if (strand == "+") {
      list(
        start = as.integer(acceptor_anchor - upstream_nt + 1),
        end = as.integer(acceptor_anchor + downstream_nt),
        acceptor_boundary_genomic = as.integer(acceptor_anchor + 1),
        rc = FALSE
      )
    } else if (strand == "-") {
      list(
        start = as.integer(acceptor_anchor - downstream_nt + 1),
        end = as.integer(acceptor_anchor + upstream_nt),
        acceptor_boundary_genomic = as.integer(acceptor_anchor),
        rc = TRUE
      )
    } else {
      list(start = NA_integer_, end = NA_integer_, acceptor_boundary_genomic = NA_integer_, rc = NA)
    }
  }
  
  extract_one_site <- function(fa, chr, strand, acceptor_anchor, site_type, row_meta, config) {
    win <- compute_window(
      acceptor_anchor = acceptor_anchor,
      strand = strand,
      upstream_nt = config$upstream_nt,
      downstream_nt = config$downstream_nt
    )
    
    raw <- get_seq_safe(
      fa = fa,
      chr = chr,
      start = win$start,
      end = win$end,
      chr_name_mode = config$chr_name_mode
    )
    
    final_seq <- if (isTRUE(raw$ok)) orient_acceptor_window(raw$seq, strand) else NA_character_
    
    tibble(
      ID = row_meta$ID,
      geneSymbol = row_meta$geneSymbol,
      cohort = row_meta$cohort,
      AS_type = row_meta$AS_type,
      feature_label = row_meta$feature_label,
      feature_label_short = row_meta$feature_label_short,
      chr = chr,
      chr_used_in_fasta = raw$chr_used,
      strand = strand,
      event_class = row_meta$event_class,
      site_type = site_type,
      species = config$species,
      genome_build = config$genome_build,
      group1_label = row_meta$group1_label,
      group2_label = row_meta$group2_label,
      IncLevelDifference = row_meta$IncLevelDifference,
      PValue = row_meta$PValue,
      FDR = row_meta$FDR,
      mean1 = row_meta$mean1,
      mean2 = row_meta$mean2,
      raw_window_start_genomic = win$start,
      raw_window_end_genomic = win$end,
      acceptor_boundary_genomic = win$acceptor_boundary_genomic,
      is_reverse_complemented = win$rc,
      raw_genomic_sequence = if (isTRUE(raw$ok)) raw$seq else NA_character_,
      final_oriented_sequence = final_seq,
      sequence_length = ifelse(is.na(final_seq), NA_integer_, nchar(final_seq)),
      extraction_ok = isTRUE(raw$ok) & !is.na(final_seq),
      extraction_reason = ifelse(isTRUE(raw$ok), NA_character_, raw$reason)
    )
  }
  
  # ---------- 7C. Event-specific extraction ----------
  extract_se_acceptor <- function(row_df, fa, config) {
    row_df <- as.list(row_df)
    strand <- row_df$strand
    chr <- row_df$chr
    
    acceptor_anchor <- if (strand == "+") {
      as.integer(row_df$exonStart_0base)
    } else if (strand == "-") {
      as.integer(row_df$exonEnd)
    } else {
      NA_integer_
    }
    
    extract_one_site(
      fa = fa,
      chr = chr,
      strand = strand,
      acceptor_anchor = acceptor_anchor,
      site_type = "cassette_acceptor",
      row_meta = row_df,
      config = config
    )
  }
  
  extract_a3ss_acceptors <- function(row_df, fa, config) {
    row_df <- as.list(row_df)
    strand <- row_df$strand
    chr <- row_df$chr
    
    if (strand == "+") {
      long_anchor <- as.integer(row_df$longExonStart_0base)
      short_anchor <- as.integer(row_df$shortES)
    } else if (strand == "-") {
      long_anchor <- as.integer(row_df$longExonEnd)
      short_anchor <- as.integer(row_df$shortEE)
    } else {
      long_anchor <- NA_integer_
      short_anchor <- NA_integer_
    }
    
    bind_rows(
      extract_one_site(fa, chr, strand, long_anchor, "long_acceptor", row_df, config),
      extract_one_site(fa, chr, strand, short_anchor, "short_acceptor", row_df, config)
    )
  }
  
  extract_sequences_all <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sequence_analysis", "sequence_tables"))
    
    fa <- open_fasta(config$genome_fasta)
    on.exit(close_fasta_safe(fa), add = TRUE)
    
    motif_df <- df %>%
      filter(AS_type %in% config$motif_as_types, event_class != "unassigned")
    
    if (isTRUE(config$keep_primary_chromosomes_only)) {
      primary_chr <- c(
        paste0("chr", 1:22), "chrX", "chrY", "chrM",
        as.character(1:22), "X", "Y", "MT", "M"
      )
      
      motif_df <- motif_df %>%
        filter(chr %in% primary_chr)
    }
    
    if (nrow(motif_df) == 0) {
      warning("No SE/A3SS events with assigned classes found for sequence extraction.")
      return(tibble())
    }
    
    seq_rows <- lapply(seq_len(nrow(motif_df)), function(i) {
      row_df <- motif_df[i, ]
      as_type <- as.character(row_df$AS_type)
      
      if (as_type == "SE") {
        extract_se_acceptor(row_df, fa, config)
      } else if (as_type == "A3SS") {
        extract_a3ss_acceptors(row_df, fa, config)
      } else {
        NULL
      }
    })
    
    seq_tbl <- bind_rows(seq_rows) %>%
      mutate(
        AS_type = factor(AS_type, levels = all_as_order),
        event_class = make_class_factor(as.character(event_class)),
        site_type = factor(site_type, levels = site_type_order)
      )
    
    write_table_files(seq_tbl, file.path(out_dir, "all_cohorts.acceptor_sequences"))
    
    seq_tbl
  }
  
  # ---------- 7D. Sequence exports ----------
  export_a3ss_audit <- function(df, seq_df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sequence_analysis", "audit_tables"))
    
    a3ss_events <- df %>%
      filter(AS_type == "A3SS", event_class != "unassigned") %>%
      mutate(.row_id = row_number())
    
    if (nrow(a3ss_events) == 0) return(invisible(NULL))
    
    sampled_events <- a3ss_events %>%
      group_by(cohort) %>%
      slice_head(n = config$a3ss_audit_n_per_cohort) %>%
      ungroup() %>%
      select(
        ID, geneSymbol, cohort, strand, event_class,
        longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE,
        feature_label, feature_label_short
      )
    
    sampled_seq <- seq_df %>%
      filter(AS_type == "A3SS") %>%
      semi_join(sampled_events, by = c("ID", "cohort")) %>%
      select(
        ID, cohort, site_type,
        raw_window_start_genomic, raw_window_end_genomic,
        acceptor_boundary_genomic, is_reverse_complemented,
        raw_genomic_sequence, final_oriented_sequence,
        extraction_ok, extraction_reason
      )
    
    audit <- sampled_events %>%
      left_join(sampled_seq, by = c("ID", "cohort")) %>%
      arrange(cohort, ID, site_type)
    
    write_table_files(audit, file.path(out_dir, "A3SS_acceptor_site_audit"))
    invisible(audit)
  }
  
  write_fasta_group <- function(df, file) {
    if (nrow(df) == 0) return(invisible(NULL))
    df <- df %>% filter(extraction_ok)
    if (nrow(df) == 0) return(invisible(NULL))
    
    headers <- paste0(
      ">",
      df$cohort, "|", df$AS_type, "|", df$site_type, "|",
      df$event_class, "|", df$ID, "|", df$geneSymbol
    )
    lines <- as.vector(rbind(headers, df$final_oriented_sequence))
    writeLines(lines, con = file)
  }
  
  export_fasta_sets <- function(seq_df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sequence_analysis", "fasta_inputs"))
    
    split_key <- seq_df %>%
      filter(event_class != "unassigned") %>%
      mutate(group_key = case_when(
        AS_type == "SE" ~ paste(cohort, AS_type, event_class, sep = "__"),
        AS_type == "A3SS" ~ paste(cohort, AS_type, site_type, event_class, sep = "__"),
        TRUE ~ paste(cohort, AS_type, site_type, event_class, sep = "__")
      ))
    
    split_list <- split(split_key, split_key$group_key, drop = TRUE)
    walk2(split_list, names(split_list), function(one_df, nm) {
      out_file <- file.path(out_dir, paste0(sanitize_filename(nm), ".fa"))
      write_fasta_group(one_df, out_file)
    })
  }
  
  # ---------- 7E. Sequence-focused plots ----------
  plot_sequence_extraction_summary <- function(seq_df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sequence_analysis", "extraction_summary"))
    if (nrow(seq_df) == 0) return(invisible(NULL))
    
    tab <- seq_df %>%
      count(cohort, AS_type, site_type, event_class, extraction_ok) %>%
      mutate(status = ifelse(extraction_ok, "extracted", "failed"))
    
    p <- ggplot(tab, aes(x = event_class, y = n, fill = status)) +
      geom_col(position = "stack") +
      facet_grid(AS_type ~ cohort + site_type, scales = "free_y") +
      labs(title = "Sequence extraction summary", x = "Event class", y = "Rows") +
      make_plot_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_plot(p, file.path(out_dir, "sequence_extraction_summary.png"), width = 12, height = 6)
  }
  
  make_nt_frequency_table <- function(seq_df) {
    seq_use <- seq_df %>%
      filter(extraction_ok, !is.na(final_oriented_sequence)) %>%
      mutate(sequence_length = nchar(final_oriented_sequence))
    
    if (nrow(seq_use) == 0) return(tibble())
    
    all_rows <- lapply(seq_len(nrow(seq_use)), function(i) {
      row <- seq_use[i, ]
      chars <- strsplit(row$final_oriented_sequence, split = "")[[1]]
      tibble(
        cohort = row$cohort,
        AS_type = as.character(row$AS_type),
        site_type = as.character(row$site_type),
        event_class = as.character(row$event_class),
        position = seq_along(chars),
        nt = chars
      )
    })
    
    bind_rows(all_rows) %>%
      count(cohort, AS_type, site_type, event_class, position, nt) %>%
      group_by(cohort, AS_type, site_type, event_class, position) %>%
      mutate(freq = n / sum(n)) %>%
      ungroup()
  }
  
  plot_nt_frequency <- function(seq_df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sequence_analysis", "nucleotide_frequency"))
    nt_tab <- make_nt_frequency_table(seq_df)
    if (nrow(nt_tab) == 0) return(invisible(NULL))
    
    se_tab <- nt_tab %>% filter(AS_type == "SE", site_type == "cassette_acceptor")
    if (nrow(se_tab) > 0) {
      by_cohort <- split(se_tab, se_tab$cohort)
      walk2(by_cohort, names(by_cohort), function(sub, cohort_name) {
        p <- ggplot(sub, aes(x = position, y = freq, color = nt)) +
          geom_line() +
          facet_wrap(~ event_class, nrow = 1) +
          labs(title = paste0(cohort_name, ": SE cassette acceptor nucleotide frequency"),
               x = "Position (1..23)", y = "Frequency") +
          make_plot_theme()
        save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".SE.cassette_acceptor.nt_frequency.png")), width = 10, height = 4)
      })
    }
    
    a3_tab <- nt_tab %>% filter(AS_type == "A3SS", site_type %in% c("long_acceptor", "short_acceptor"))
    if (nrow(a3_tab) > 0) {
      by_cohort_site <- split(a3_tab, list(a3_tab$cohort, a3_tab$site_type), drop = TRUE)
      walk2(by_cohort_site, names(by_cohort_site), function(sub, nm) {
        cohort_name <- unique(sub$cohort)
        site_name <- unique(sub$site_type)
        p <- ggplot(sub, aes(x = position, y = freq, color = nt)) +
          geom_line() +
          facet_wrap(~ event_class, nrow = 1) +
          labs(title = paste0(cohort_name, ": A3SS ", site_name, " nucleotide frequency"),
               x = "Position (1..23)", y = "Frequency") +
          make_plot_theme()
        save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".A3SS.", site_name, ".nt_frequency.png")), width = 10, height = 4)
      })
    }
  }
  
  plot_sequence_logos <- function(seq_df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "sequence_analysis", "logos"))
    
    seq_use <- seq_df %>%
      filter(extraction_ok, !is.na(final_oriented_sequence),
             sequence_length == (config$upstream_nt + config$downstream_nt))
    
    if (nrow(seq_use) == 0) return(invisible(NULL))
    
    se_groups <- seq_use %>%
      filter(AS_type == "SE", site_type == "cassette_acceptor") %>%
      split(list(.$cohort, .$event_class), drop = TRUE)
    
    walk2(se_groups, names(se_groups), function(sub, nm) {
      if (nrow(sub) == 0) return(NULL)
      cohort_name <- unique(sub$cohort)
      class_name <- unique(as.character(sub$event_class))
      p <- ggseqlogo(sub$final_oriented_sequence, method = "bits") +
        ggtitle(paste0(cohort_name, ": SE cassette acceptor - ", class_name)) +
        theme_bw(base_size = 11)
      save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".SE.cassette_acceptor.", class_name, ".logo.png")), width = 8, height = 3)
    })
    
    a3_groups <- seq_use %>%
      filter(AS_type == "A3SS") %>%
      split(list(.$cohort, .$site_type, .$event_class), drop = TRUE)
    
    walk2(a3_groups, names(a3_groups), function(sub, nm) {
      if (nrow(sub) == 0) return(NULL)
      cohort_name <- unique(sub$cohort)
      site_name <- unique(as.character(sub$site_type))
      class_name <- unique(as.character(sub$event_class))
      p <- ggseqlogo(sub$final_oriented_sequence, method = "bits") +
        ggtitle(paste0(cohort_name, ": A3SS ", site_name, " - ", class_name)) +
        theme_bw(base_size = 11)
      save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".A3SS.", site_name, ".", class_name, ".logo.png")), width = 8, height = 3)
    })
  }
  
  # ==========================================================
  # Part 8. General summary plots for all AS types
  # ==========================================================
  make_event_count_table <- function(df, config) {
    df %>%
      filter(event_class != "unassigned") %>%
      mutate(
        AS_type = factor(AS_type, levels = all_as_order),
        event_class = make_class_factor(as.character(event_class))
      ) %>%
      count(cohort, AS_type, event_class, .drop = FALSE) %>%
      mutate(n_plot = ifelse(n < config$min_label_count, NA_integer_, n))
  }
  
  make_filter_funnel_table <- function(df) {
    all_types <- levels(factor(df$AS_type, levels = all_as_order))
    
    by_type <- bind_rows(lapply(all_types, function(one_type) {
      sub <- df %>% filter(AS_type == one_type)
      tibble(
        AS_type = one_type,
        stage = c("total_events", "pass_count_filter", "pass_deltaPSI_filter", "significant_motif_set", "nonsignificant_background_set"),
        n = c(
          nrow(sub),
          sum(sub$filter_counts, na.rm = TRUE),
          sum(sub$filter_counts & sub$filter_PSI, na.rm = TRUE),
          sum(sub$is_significant_for_motif, na.rm = TRUE),
          sum(sub$is_background_for_motif, na.rm = TRUE)
        )
      )
    }))
    
    overall <- tibble(
      AS_type = "ALL",
      stage = c("total_events", "pass_count_filter", "pass_deltaPSI_filter", "significant_motif_set", "nonsignificant_background_set"),
      n = c(
        nrow(df),
        sum(df$filter_counts, na.rm = TRUE),
        sum(df$filter_counts & df$filter_PSI, na.rm = TRUE),
        sum(df$is_significant_for_motif, na.rm = TRUE),
        sum(df$is_background_for_motif, na.rm = TRUE)
      )
    )
    
    bind_rows(overall, by_type) %>%
      mutate(
        AS_type = factor(AS_type, levels = c("ALL", all_as_order)),
        stage = factor(stage, levels = c("total_events", "pass_count_filter", "pass_deltaPSI_filter", "significant_motif_set", "nonsignificant_background_set"))
      )
  }
  
  plot_general_overview_all_as <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "general_overview"))
    tab <- make_event_count_table(df, config)
    
    p <- ggplot(tab, aes(x = AS_type, y = n, fill = event_class)) +
      geom_col(position = "stack") +
      facet_wrap(~ cohort, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = class_colors, drop = FALSE) +
      labs(title = "General splicing overview across all AS types", x = "AS type", y = "Events") +
      make_plot_theme()
    
    save_plot(p, file.path(out_dir, "general_splicing_overview_all_AS_types.png"), width = 10, height = 4)
    p
  }
  
  plot_event_counts_across_all_as_types <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "general_overview"))
    tab <- make_event_count_table(df, config)
    
    p <- ggplot(tab, aes(x = cohort, y = n, fill = event_class)) +
      geom_col(position = "stack") +
      geom_text(aes(label = n_plot), position = position_stack(vjust = 0.5), color = "white", size = 3, na.rm = TRUE) +
      facet_wrap(~ AS_type, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = class_colors, drop = FALSE) +
      labs(title = "Event count plot across all AS types", x = "", y = "Events") +
      make_plot_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_plot(p, file.path(out_dir, "event_count_across_all_AS_types.png"), width = 12, height = 4)
    p
  }
  
  plot_per_as_type_event_counts <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "per_AS_type", "event_counts"))
    tabs <- make_event_count_table(df, config)
    
    lapply(all_as_order, function(one_type) {
      sub <- tabs %>% filter(AS_type == one_type)
      p <- ggplot(sub, aes(x = cohort, y = n, fill = event_class)) +
        geom_col(position = "stack") +
        geom_text(aes(label = n_plot), position = position_stack(vjust = 0.5), color = "white", size = 3, na.rm = TRUE) +
        scale_fill_manual(values = class_colors, drop = FALSE) +
        labs(title = paste0(one_type, ": event counts by class"), x = "", y = "Events") +
        make_plot_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      save_plot(p, file.path(out_dir, paste0(one_type, ".event_counts.png")), width = 6, height = 4)
    })
  }
  
  plot_filter_funnel <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "filter_funnel"))
    
    by_cohort <- split(df, df$cohort)
    lapply(names(by_cohort), function(cohort_name) {
      sub <- by_cohort[[cohort_name]]
      tab <- make_filter_funnel_table(sub)
      
      p_all <- ggplot(tab %>% filter(AS_type == "ALL"), aes(x = stage, y = n, group = 1)) +
        geom_line() +
        geom_point(size = 2) +
        labs(title = paste0(cohort_name, ": filter funnel (all AS types)"), x = "Stage", y = "Events") +
        make_plot_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      save_plot(p_all, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".filter_funnel.ALL.png")), width = 7, height = 4)
      
      p_by_type <- ggplot(tab %>% filter(AS_type != "ALL"), aes(x = stage, y = n, group = AS_type, color = AS_type)) +
        geom_line() +
        geom_point(size = 1.8) +
        facet_wrap(~ AS_type, nrow = 1, scales = "free_y") +
        labs(title = paste0(cohort_name, ": filter funnel by AS type"), x = "Stage", y = "Events") +
        make_plot_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      
      save_plot(p_by_type, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".filter_funnel.by_AS_type.png")), width = 12, height = 4)
    })
  }
  
  plot_dpsi_density_by_class <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "psi_density"))
    
    df2 <- df %>%
      filter(event_class != "unassigned") %>%
      mutate(
        event_class = make_class_factor(as.character(event_class)),
        AS_type = factor(AS_type, levels = all_as_order)
      )
    
    by_cohort <- split(df2, df2$cohort)
    lapply(names(by_cohort), function(cohort_name) {
      sub <- by_cohort[[cohort_name]]
      p <- ggplot(sub, aes(x = IncLevelDifference, fill = event_class, color = event_class)) +
        geom_density(alpha = 0.20) +
        facet_wrap(~ AS_type, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = class_colors, drop = FALSE) +
        scale_color_manual(values = class_colors, drop = FALSE) +
        labs(title = paste0(cohort_name, ": dPSI density by class"), x = "IncLevelDifference (group1 - group2)", y = "Density") +
        make_plot_theme()
      
      save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".dPSI_density.by_class.png")), width = 12, height = 4)
    })
  }
  
  plot_psi_scatter_all_types <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "psi_scatter"))
    
    df2 <- df %>%
      filter(event_class != "unassigned") %>%
      mutate(
        event_class = make_class_factor(as.character(event_class)),
        AS_type = factor(AS_type, levels = all_as_order)
      )
    
    by_cohort <- split(df2, df2$cohort)
    lapply(names(by_cohort), function(cohort_name) {
      sub <- by_cohort[[cohort_name]]
      
      p_all <- ggplot(sub, aes(x = mean1, y = mean2)) +
        geom_point(color = "grey80", size = config$point_size, alpha = config$point_alpha) +
        geom_point(data = sub %>% filter(event_class == "more_skipped_in_group2"), aes(color = event_class), size = config$point_size, alpha = config$point_alpha) +
        geom_point(data = sub %>% filter(event_class == "more_included_in_group2"), aes(color = event_class), size = config$point_size, alpha = config$point_alpha) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        facet_wrap(~ AS_type, nrow = 1) +
        scale_color_manual(values = class_colors, drop = FALSE) +
        labs(title = paste0(cohort_name, ": PSI scatter across all AS types"),
             x = paste0("mean PSI (", config$group1_label, ")"),
             y = paste0("mean PSI (", config$group2_label, ")")) +
        make_plot_theme()
      
      save_plot(p_all, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".PSI_scatter.all_AS_types.png")), width = 12, height = 4)
      
      walk(all_as_order, function(one_type) {
        sub_type <- sub %>% filter(AS_type == one_type)
        if (nrow(sub_type) == 0) return(NULL)
        p <- ggplot(sub_type, aes(x = mean1, y = mean2)) +
          geom_point(color = "grey80", size = config$point_size, alpha = config$point_alpha) +
          geom_point(data = sub_type %>% filter(event_class == "more_skipped_in_group2"), aes(color = event_class), size = config$point_size, alpha = config$point_alpha) +
          geom_point(data = sub_type %>% filter(event_class == "more_included_in_group2"), aes(color = event_class), size = config$point_size, alpha = config$point_alpha) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
          scale_color_manual(values = class_colors, drop = FALSE) +
          labs(title = paste0(cohort_name, ": ", one_type, " PSI scatter"),
               x = paste0("mean PSI (", config$group1_label, ")"),
               y = paste0("mean PSI (", config$group2_label, ")")) +
          make_plot_theme()
        save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".", one_type, ".PSI_scatter.png")), width = 5, height = 5)
      })
    })
  }
  
  plot_dpsi_vs_fdr_all_types <- function(df, config) {
    out_dir <- ensure_dir(file.path(config$output_root, "summary_plots", "dpsi_vs_fdr"))
    
    df2 <- df %>%
      mutate(
        neglog10FDR = -log10(pmax(FDR, 1e-300)),
        event_class = make_class_factor(as.character(event_class)),
        AS_type = factor(AS_type, levels = all_as_order)
      )
    
    by_cohort <- split(df2, df2$cohort)
    lapply(names(by_cohort), function(cohort_name) {
      sub <- by_cohort[[cohort_name]]
      
      p_all <- ggplot(sub, aes(x = IncLevelDifference, y = neglog10FDR, color = event_class)) +
        geom_point(size = config$point_size, alpha = config$point_alpha) +
        geom_vline(xintercept = c(-config$cutoff_delta_psi, config$cutoff_delta_psi), linetype = "dashed") +
        geom_hline(yintercept = -log10(config$cutoff_fdr), linetype = "dashed") +
        facet_wrap(~ AS_type, nrow = 1) +
        scale_color_manual(values = class_colors, drop = FALSE) +
        labs(title = paste0(cohort_name, ": dPSI vs -log10(FDR)"),
             x = "IncLevelDifference (group1 - group2)", y = "-log10(FDR)") +
        make_plot_theme()
      
      save_plot(p_all, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".dPSI_vs_FDR.all_AS_types.png")), width = 12, height = 4)
      
      walk(all_as_order, function(one_type) {
        sub_type <- sub %>% filter(AS_type == one_type)
        if (nrow(sub_type) == 0) return(NULL)
        p <- ggplot(sub_type, aes(x = IncLevelDifference, y = neglog10FDR, color = event_class)) +
          geom_point(size = config$point_size, alpha = config$point_alpha) +
          geom_vline(xintercept = c(-config$cutoff_delta_psi, config$cutoff_delta_psi), linetype = "dashed") +
          geom_hline(yintercept = -log10(config$cutoff_fdr), linetype = "dashed") +
          scale_color_manual(values = class_colors, drop = FALSE) +
          labs(title = paste0(cohort_name, ": ", one_type, " dPSI vs -log10(FDR)"),
               x = "IncLevelDifference", y = "-log10(FDR)") +
          make_plot_theme()
        save_plot(p, file.path(out_dir, paste0(sanitize_filename(cohort_name), ".", one_type, ".dPSI_vs_FDR.png")), width = 5, height = 5)
      })
    })
  }
  
  # ==========================================================
  # Part 8B. Sashimi-only branch from cached master
  # ==========================================================
  read_gene_list <- function(gene_symbols = NULL, gene_file = NULL) {
    genes <- character()

    if (!is.null(gene_symbols)) {
      genes <- c(genes, as.character(gene_symbols))
    }

    if (!is.null(gene_file)) {
      if (!file.exists(gene_file)) stop("Gene-list file not found: ", gene_file)
      raw <- readLines(gene_file, warn = FALSE)
      raw <- raw[nzchar(trimws(raw)) & !grepl("^#", trimws(raw))]
      from_file <- vapply(strsplit(raw, "[\\t,]"), function(x) trimws(x[[1]]), character(1))
      # Remove a likely header while retaining ordinary gene symbols.
      from_file <- from_file[!tolower(from_file) %in% c("gene", "genes", "genesymbol", "gene_symbol")]
      genes <- c(genes, from_file)
    }

    genes <- unique(trimws(genes))
    genes <- genes[nzchar(genes) & !is.na(genes)]
    if (length(genes) == 0) stop("No valid genes were supplied.")
    genes
  }

  make_event_key <- function(df) {
    sanitize_filename(paste(
      df$cohort, as.character(df$AS_type), df$geneSymbol,
      df$feature_label, sep = "__"
    ))
  }

  build_sashimi_run_manifest <- function(df, genes, config, sashimi_config, run_dir) {
    manifest_dir <- ensure_dir(file.path(run_dir, "manifest"))
    event_root <- ensure_dir(file.path(run_dir, "event_inputs"))
    plots_root <- ensure_dir(file.path(run_dir, "plots"))

    requested_summary <- tibble(geneSymbol = genes) %>%
      left_join(
        df %>% count(geneSymbol, name = "n_events_before_filtering"),
        by = "geneSymbol"
      ) %>%
      mutate(
        n_events_before_filtering = replace_na(n_events_before_filtering, 0L),
        found_in_master = n_events_before_filtering > 0
      )

    selected <- do.call(
      filter_events_general,
      list(
        df = df,
        gene_symbols = genes,
        AS_type = sashimi_config$AS_type,
        cohort = sashimi_config$cohort,
        filter_counts = sashimi_config$filter_counts,
        filter_PSI = sashimi_config$filter_PSI,
        filter_PValue = sashimi_config$filter_PValue,
        filter_FDR = sashimi_config$filter_FDR,
        event_class = sashimi_config$event_class,
        annotation_status = NULL,
        min_abs_dpsi = sashimi_config$min_abs_dpsi,
        max_fdr = sashimi_config$max_fdr
      )
    )

    if (nrow(selected) == 0) {
      write_table_files(requested_summary, file.path(manifest_dir, "requested_gene_summary"))
      stop("No events passed the sashimi-only branch filters.")
    }

    gtf_ref_file <- sashimi_config$gtf_reference_rds
    if (!file.exists(gtf_ref_file)) stop("Cached GTF reference not found: ", gtf_ref_file)
    gtf_ref <- readRDS(gtf_ref_file)

    # Annotate only the filtered rows for the requested genes.
    selected <- annotate_event_display(
      selected,
      gtf_ref = gtf_ref,
      tol = config$exon_match_tolerance,
      progress_step = 500,
      annotate_only_significant = FALSE,
      gene_symbols = genes
    )

    if (!is.null(sashimi_config$annotation_status)) {
      selected <- selected %>%
        filter(annotation_status %in% sashimi_config$annotation_status)
    }
    if (nrow(selected) == 0) stop("No events remained after annotation-status filtering.")

    selected <- selected %>%
      arrange(cohort, AS_type, geneSymbol, event_display, feature_label) %>%
      mutate(
        run_name = sashimi_config$run_name,
        event_key = sanitize_filename(paste(
          cohort, as.character(AS_type), geneSymbol, feature_label, sep = "__"
        )),
        plot_id = sprintf("plot_%05d", row_number())
      )

    # Reuse the original one-row rMATS writer, then move its output contract
    # into this run-specific directory.
    selected <- write_one_row_rmats_inputs(selected, run_dir)
    selected <- selected %>%
      mutate(
        rmats_event_file = normalizePath(rmats_event_file, winslash = "/", mustWork = TRUE),
        plot_output_dir = file.path(
          normalizePath(plots_root, winslash = "/", mustWork = TRUE),
          cohort, as.character(AS_type), plot_id
        ),
        expected_pdf_dir = plot_output_dir,
        ppt_title = event_display,
        ppt_subtitle = paste0(
          cohort, " | ", AS_type, " | ", feature_label,
          " | PValue=", fmt_num(PValue, 2),
          " | FDR=", fmt_num(FDR, 2),
          " | dPSI_", group2_label, "_minus_", group1_label, "=",
          signif(delta_psi_group2_minus_group1, 3)
        ),
        plot_status = "pending"
      )

    after_counts <- selected %>% count(geneSymbol, name = "n_events_selected")
    requested_summary <- requested_summary %>%
      left_join(after_counts, by = "geneSymbol") %>%
      mutate(n_events_selected = replace_na(n_events_selected, 0L))

    write_table_files(requested_summary, file.path(manifest_dir, "requested_gene_summary"))
    write_table_files(
      requested_summary %>% filter(!found_in_master),
      file.path(manifest_dir, "missing_genes")
    )
    write_table_files(selected, file.path(manifest_dir, "sashimi_manifest"))

    selected
  }

  run_sashimi_only <- function(config = CONFIG, sashimi_config = SASHIMI_CONFIG) {
    genes <- read_gene_list(
      gene_symbols = sashimi_config$gene_symbols,
      gene_file = sashimi_config$gene_file
    )

    master_file <- sashimi_config$annotated_master_rds
    if (!file.exists(master_file)) stop("Cached annotated master not found: ", master_file)

    run_dir <- ensure_dir(file.path(
      sashimi_config$runs_root,
      sanitize_filename(sashimi_config$run_name)
    ))
    walk(c("manifest", "event_inputs", "plots", "reports", "logs", "ppt"),
         ~ ensure_dir(file.path(run_dir, .x)))

    message("Sashimi-only Step 1/4: Load cached annotated master")
    df <- readRDS(master_file)
    df <- upgrade_delta_psi_columns(df, config)

    # Revalidate the cached table before any sashimi filtering.
    message("Sashimi-only Step 2/4: Validate IncLevelDifference direction")
    summarize_inclevel_direction_qc(
      df,
      config,
      out_dir = file.path(run_dir, "reports", "inclevel_direction_qc")
    )

    message("Sashimi-only Step 3/4: Filter, annotate, and write manifest")
    manifest <- build_sashimi_run_manifest(
      df = df,
      genes = genes,
      config = config,
      sashimi_config = sashimi_config,
      run_dir = run_dir
    )

    config_lines <- c(
      paste0("run_name: ", sashimi_config$run_name),
      paste0("master: ", normalizePath(master_file, winslash = "/", mustWork = TRUE)),
      paste0("genes: ", paste(genes, collapse = ",")),
      paste0("group1: ", config$group1_label),
      paste0("group2: ", config$group2_label),
      paste0("raw_rmats_difference: ", config$group1_label, " minus ", config$group2_label),
      paste0("reported_difference: ", config$group2_label, " minus ", config$group1_label),
      paste0("events_selected: ", nrow(manifest))
    )
    writeLines(config_lines, file.path(run_dir, "run_config.txt"))

    message("Sashimi-only Step 4/4: Manifest ready")
    message("Manifest: ", file.path(run_dir, "manifest", "sashimi_manifest.tsv"))
    message("Next: run 2-run_sashimi_from_manifest.sh using this manifest.")

    invisible(list(run_dir = run_dir, genes = genes, manifest = manifest))
  }

  # ==========================================================
  # Part 9. Output tree + run metadata
  # ==========================================================
  initialize_output_tree <- function(config) {
    dirs <- c(
      file.path(config$output_root, "processed", "analysis_ready"),
      file.path(config$output_root, "processed", "annotated_master"),
      file.path(config$output_root, "processed", "qc"),
      file.path(config$output_root, "reference"),
      file.path(config$output_root, "sashimi_manifest"),
      file.path(config$output_root, "summary_plots", "general_overview"),
      file.path(config$output_root, "summary_plots", "per_AS_type", "event_counts"),
      file.path(config$output_root, "summary_plots", "filter_funnel"),
      file.path(config$output_root, "summary_plots", "psi_density"),
      file.path(config$output_root, "summary_plots", "psi_scatter"),
      file.path(config$output_root, "summary_plots", "dpsi_vs_fdr"),
      file.path(config$output_root, "sequence_analysis", "sequence_tables"),
      file.path(config$output_root, "sequence_analysis", "extraction_summary"),
      file.path(config$output_root, "sequence_analysis", "logos"),
      file.path(config$output_root, "sequence_analysis", "nucleotide_frequency"),
      file.path(config$output_root, "sequence_analysis", "audit_tables"),
      file.path(config$output_root, "sequence_analysis", "fasta_inputs")
    )
    walk(dirs, ensure_dir)
  }
  
  write_run_metadata <- function(config) {
    out_file <- file.path(config$output_root, "run_metadata.txt")
    lines <- c(
      paste0("rmats_root: ", config$rmats_root),
      paste0("output_root: ", config$output_root),
      paste0("gtf_file: ", config$gtf_file),
      paste0("genome_fasta: ", config$genome_fasta),
      paste0("species: ", config$species),
      paste0("genome_build: ", config$genome_build),
      paste0("chr_name_mode: ", config$chr_name_mode),
      paste0("group1_label: ", config$group1_label),
      paste0("group2_label: ", config$group2_label),
      paste0("rmats_raw_difference: ", config$group1_label, " minus ", config$group2_label),
      paste0("canonical_downstream_difference: ", config$group2_label, " minus ", config$group1_label),
      paste0("inclevel_difference_tolerance: ", config$inclevel_difference_tolerance),
      paste0("cutoff_pvalue: ", config$cutoff_pvalue),
      paste0("cutoff_fdr: ", config$cutoff_fdr),
      paste0("cutoff_delta_psi: ", config$cutoff_delta_psi),
      paste0("min_total_counts: ", config$min_total_counts),
      paste0("upstream_nt: ", config$upstream_nt),
      paste0("downstream_nt: ", config$downstream_nt),
      paste0("cohorts: ", paste(names(config$cohorts), config$cohorts, sep = "=", collapse = "; "))
    )
    writeLines(lines, con = out_file)
  }
  
  # ==========================================================
  # Part 10. Unified driver
  # ==========================================================
  
  config = CONFIG
  
  run_pipeline_integrated <- function(config = CONFIG) {
    initialize_output_tree(config)
    write_run_metadata(config)
    
    message("Step 1/8: Build analysis-ready table")
    df <- build_analysis_tables_all(config)
    
    message("Step 2/8: Classify events")
    df <- classify_events(df, config)
    
    message("Step 3/8: Build or load cached GTF reference")
    gtf_ref <- build_or_load_gtf_reference(
      gtf_file = config$gtf_file,
      output_root = config$output_root,
      force_rebuild = FALSE
    )
    
    message("Step 4/8: Annotate event_display from longest transcript")
    df <- annotate_event_display(
      df,
      gtf_ref = gtf_ref,
      tol = config$exon_match_tolerance,
      progress_step = 1000,
      annotate_only_significant = TRUE,
      gene_symbols = config$plot_filters$gene_symbols
    )
    
    annotated_dir <- ensure_dir(file.path(config$output_root, "processed", "annotated_master"))
    write_table_files(df, file.path(annotated_dir, "all_cohorts.annotated_master"))
    
    # 
    # df <- readRDS("./rmats_integrated_output-U2AF1/processed/annotated_master/all_cohorts.annotated_master.rds")
    
    message("Step 5/8: General summary plots")
    plot_general_overview_all_as(df, config)
    plot_event_counts_across_all_as_types(df, config)
    plot_per_as_type_event_counts(df, config)
    plot_filter_funnel(df, config)
    plot_dpsi_density_by_class(df, config)
    plot_psi_scatter_all_types(df, config)
    plot_dpsi_vs_fdr_all_types(df, config)
    
    message("Step 6/8: Build sashimi manifest and one-row event files")
    manifest <- build_sashimi_manifest(df, config)
    
    message("Step 7/8: Sequence extraction for SE + A3SS")
    seq_df <- extract_sequences_all(df, config)
    
    message("Step 8/8: Sequence outputs and motif plots")
    export_a3ss_audit(df, seq_df, config)
    export_fasta_sets(seq_df, config)
    plot_sequence_extraction_summary(seq_df, config)
    plot_nt_frequency(seq_df, config)
    plot_sequence_logos(seq_df, config)
    
    message("Done.")
    invisible(list(
      annotated_master = df,
      sashimi_manifest = manifest,
      gtf_reference = gtf_ref,
      extracted_sequences = seq_df
    ))
  }
  
}
# ==========================================================
# Part 11. Optional run block
# ==========================================================
# Choose ONE mode and uncomment it.
#
# A) Full downstream pipeline (raw rMATS tables changed):
# result <- run_pipeline_integrated(CONFIG)
#
# B) Sashimi-only branch (new gene list; cached master already exists):
sashimi_result <- run_sashimi_only(CONFIG, SASHIMI_CONFIG)
#
# Common sequence-extraction fixes:
# - If FASTA uses names like '1' but rMATS uses 'chr1', keep chr_name_mode='auto'
# - If needed, force chr_name_mode='add_chr' or 'drop_chr'
# - Create FASTA index first:
#     Rsamtools::indexFa(CONFIG$genome_fasta)
