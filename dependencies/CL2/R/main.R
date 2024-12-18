


#' @title Analyze a single CyTOF file
#' @description Read and parse file, cluster and label main cell types, detect
#' debris and doublets, gate downstream T and B cell subpopulations.
#' @param file FCS file name.
#' @param dir_in A path to read fcs files from.
#' @param dir_out A path to write output.
#' @param cols Columns to use in analysis.
#' @returns Nothing; all output written to file.
#' @export
analyze_cytof_file <- function(file, dir_in, dir_out, cols) {

  dir_in <- paste0(dir_in, "/")
  dir_out <- paste0(dir_out, "/")

  set.seed(0)
  path <- paste0(dir_in, file)
  fn <- file %>%
    str_remove("_Normalized.fcs") %>%
    str_remove("_normalized.fcs") %>%
    str_remove("_Processed.fcs") %>%
    str_remove(".fcs")
  print("##################")
  print(paste("Starting processing for file", fn, "..."))

  ### Read FCS file and pregate on Bead/Gaussian/LiveDead channels ###
  ff <- read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE)
  df_raw <- ff %>%
    read_data() %>%
    as_tibble()
  pregating <- pregate_data(df_raw, fn, dir_out, plot=TRUE)
  event_type <- pregating$event_type
  channels_remove <- "Time|Bead|Live|Center|Offset|Residual|Width"
  channels_keep <- names(df_raw)[!grepl(channels_remove, names(df_raw))]

  ### Filter and scale data ###
  df <- df_raw %>%
    filter(event_type == "") %>%
    select(all_of(channels_keep))

  # df <- df[sample(nrow(df), min(nrow(df), 1e5)),]
  names(df) <- str_replace(names(df), "CD8a/CD64", "CD8a")

  x_full <- df %>% as.matrix() %>% scale_data()
  x <- x_full[,cols]

  ### Build UMAP on a subset of cells ###
  n_sel <- min(nrow(df), 3e4)
  sel_umap <- sample(nrow(df), n_sel)
  df_um <- get_umap(df, x, sel_umap)

  ### Classify cells using pretrained logistic regression model ###
  pred_probs <- predict_cell_type(x_full, defs_major, return_probs = TRUE)
  cell_type <- colnames(pred_probs)[unname(apply(pred_probs, 1, which.max))]
  confidence <- apply(pred_probs, 1, max)
  uncertain_cells <- which(confidence < 0.9 & !grepl("Debris|Neutr", cell_type))
  cell_type[uncertain_cells] <- "Uncertain"

  ### Detect doublets ###
  print(paste("Doublet detection for file", fn, "..."))
  cleanet_res <- detect_doublets(df, cols, cell_type, dir_out, fn)
  cell_type <- cleanet_res$cell_type

  write_stats_file(pregating$df_stats, cell_type, dir_out, fn)

  ### Visualize main cell types ###
  df_viz <- df_um %>% mutate(ct = cell_type[sel_umap])

  p_markers <- plot_umap_markers(df_viz, fn)
  ggsave(p_markers, filename = paste0(dir_out, "umap_markers/", fn, ".png"),
         width=16, height=10)

  p_major <- plot_umap_major_only(df_viz, fn)
  ggsave(p_major, filename = paste0(dir_out, "umap_major/", fn, ".png"),
         width=9, height=7)

  p_dna <- plot_dna_cd45(df_viz, fn)
  ggsave(p_dna, filename = paste0(dir_out, "DNA_CD45/", fn, ".png"),
         width=9, height=7)

  ### Refine T cell subsets ###
  tcells <- which(cell_type=="T cell")
  n_tcells <- length(tcells)

  if(n_tcells > 20) {
    cols_tcell <- c("CD4", "CD8a", "TCRgd")
    x_tcell <- x_full[tcells,cols_tcell]
    pred_tcell <- predict_cell_type(x_tcell, defs_tcell)
    cell_type[tcells] <- pred_tcell
  }

  ### Refine monocytes and mdc ###
  mono <- which(cell_type=="Myeloid")
  n_mono <- length(mono)

  if (n_mono > 20) {
    cols_mono <- c("CD11c", "CD14", "CD38", "CD123", "CD294", "HLA-DR", "CD45RA")
    x_mono <- x_full[mono,cols_mono]
    pred_mono <- predict_cell_type(x_mono, defs_myel)
    cell_type[mono] <- pred_mono
  }

  event_type[which(event_type=="")] <- cell_type
  cell_idx <- which(!grepl("Debris|_|Bead|Offset|Residual|Width|Center|Dead", event_type))
  ff_clean <- ff[cell_idx,]

  write.FCS(ff_clean, paste0(dir_out, "fcs_clean/", fn, ".fcs"))
  df_file <- tibble(event_type=event_type)
  write_csv(df_file, file=paste0(dir_out, "files_labeled/", fn, ".csv"), progress=FALSE)

  backgate_major(df, cell_type, dir_out, fn)

  print(paste("Gating file", fn, "..."))
  gate_detailed_phenos(df, x_full, cell_type, dir_out, fn)

  ### Compute density estimates to be used later for QC ###
  valid_cell_types <- c("Neutrophil", "Eosinophil", "Basophil", "B cell",
                        "Myeloid", "T cell CD4", "T cell CD8", "T cell gd")
  channels <- setdiff(names(df), c("DNA1", "DNA2", "Event_length"))
  df_kdes <- estimate_distributions(cell_type, df, fn, channels, valid_cell_types)
  write_csv(df_kdes, file=paste0(dir_out, "/kdes_for_qc/", fn, ".csv"), progress=FALSE)

  print(paste("Finished file", fn, "!"))
  return(NULL)
}


#' @title Create output directory structure
#' @description Create subdirectories for saving figures and csv output.
#' @param dir_out Base output path.
#' @returns Nothing.
#' @export
create_dirs <- function(dir_out) {
  dir.create(paste0(dir_out, "/files_labeled"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/fcs_clean"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/feat_major"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/feat_adaptive"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/kdes_for_qc"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/cleanup_gates"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/cleanup_stats"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/umap_major"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/umap_markers"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/DNA_CD45"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/doublet_csv"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/doublet_fig"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/gating"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/thresholds"), showWarnings = FALSE)

  for (gate_group in unique(gate_hierarchy_full$GateGroup)) {
    dir.create(paste0(dir_out, "/gating/", gate_group), showWarnings = FALSE)
  }

  dir.create(paste0(dir_out, "/backgating"), showWarnings = FALSE)
  for (n in c("T cell CD4 Naive", "T cell CD8 Naive", "CD45 CD66b",
              "CD4 CD8a", "CD3 CD19", "CD11c CD14", "CD3 CD56",
              "CD16 CD66b", "CD123 CD294", "CD3 TCRgd", "CD14 CD38")) {
    dir.create(paste0(dir_out, "/backgating/", n), showWarnings = FALSE)
  }
}



