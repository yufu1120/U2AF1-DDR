# make_relative_frequency_logos_only.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(ggseqlogo)
})

CONFIG <- list(
  output_root = "./../rmats_integrated_output",
  upstream_nt = 20,
  downstream_nt = 3
)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

sanitize_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

save_plot <- function(plot_obj, filename, width = 8, height = 3, dpi = 300) {
  ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, units = "in")
}

plot_sequence_logos_relative_frequency_only <- function(config = CONFIG) {
  seq_rds <- file.path(
    config$output_root,
    "sequence_analysis",
    "sequence_tables",
    "all_cohorts.acceptor_sequences.rds"
  )
  
  if (!file.exists(seq_rds)) {
    stop("Cannot find cached sequence table: ", seq_rds)
  }
  
  seq_df <- readRDS(seq_rds)
  
  out_dir <- ensure_dir(file.path(
    config$output_root,
    "sequence_analysis",
    "logos_relative_frequency"
  ))
  
  seq_use <- seq_df %>%
    filter(
      extraction_ok,
      !is.na(final_oriented_sequence),
      sequence_length == (config$upstream_nt + config$downstream_nt)
    )
  
  if (nrow(seq_use) == 0) {
    stop("No usable extracted sequences found.")
  }
  
  # SE cassette acceptor
  se_groups <- seq_use %>%
    filter(AS_type == "SE", site_type == "cassette_acceptor") %>%
    split(list(.$cohort, .$event_class), drop = TRUE)
  
  walk2(se_groups, names(se_groups), function(sub, nm) {
    if (nrow(sub) == 0) return(NULL)
    
    cohort_name <- unique(sub$cohort)
    class_name <- unique(as.character(sub$event_class))
    
    p <- ggseqlogo(sub$final_oriented_sequence, method = "prob") +
      ggtitle(paste0(
        cohort_name,
        ": SE cassette acceptor - ",
        class_name,
        " [relative frequency]"
      )) +
      ylab("Relative frequency") +
      theme_bw(base_size = 11)
    
    save_plot(
      p,
      file.path(
        out_dir,
        paste0(
          sanitize_filename(cohort_name),
          ".SE.cassette_acceptor.",
          class_name,
          ".relative_frequency_logo.png"
        )
      )
    )
  })
  
  # A3SS long/short acceptors
  a3_groups <- seq_use %>%
    filter(AS_type == "A3SS") %>%
    split(list(.$cohort, .$site_type, .$event_class), drop = TRUE)
  
  walk2(a3_groups, names(a3_groups), function(sub, nm) {
    if (nrow(sub) == 0) return(NULL)
    
    cohort_name <- unique(sub$cohort)
    site_name <- unique(as.character(sub$site_type))
    class_name <- unique(as.character(sub$event_class))
    
    p <- ggseqlogo(sub$final_oriented_sequence, method = "prob") +
      ggtitle(paste0(
        cohort_name,
        ": A3SS ",
        site_name,
        " - ",
        class_name,
        " [relative frequency]"
      )) +
      ylab("Relative frequency") +
      theme_bw(base_size = 11)
    
    save_plot(
      p,
      file.path(
        out_dir,
        paste0(
          sanitize_filename(cohort_name),
          ".A3SS.",
          site_name,
          ".",
          class_name,
          ".relative_frequency_logo.png"
        )
      )
    )
  })
  
  message("Done. Relative-frequency logos saved to: ", out_dir)
}

plot_sequence_logos_relative_frequency_only(CONFIG)
