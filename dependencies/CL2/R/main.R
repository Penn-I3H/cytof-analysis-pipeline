

#' @title Analyze a single CyTOF file
#' @description Read and parse file, cluster and label main cell types, detect
#' debris and doublets, gate downstream T and B cell subpopulations.
#' @param file FCS file name.
#' @param dir_in A path to read fcs files from.
#' @param dir_out A path to write output.
#' @param cols Columns to use in analysis.
#' @param cofactor Parameter for asinh transformation.
#' @returns Nothing; all output written to file.
#' @export
analyze_cytof_file <- function(file, dir_in, dir_out, cols, cofactor=5) {

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
    read_data(transform = FALSE) %>%
    as_tibble()
  pregating <- pregate_data(df_raw, fn, dir_out, plot=TRUE)
  event_type <- pregating$event_type

  df <- filter_bead_gaussian_live(df_raw, event_type)
  names(df) <- str_replace(names(df), "CD8a/CD64", "CD8a") ## for extension panel
  # df <- df[sample(nrow(df), min(nrow(df), 5e4)),]

  data_transf <- df %>%
    mutate(across(!matches("length"), ~asinh(.x/cofactor))) %>%
    as.matrix()
  data_scaled <- scale_data(data_transf)

  ### Build UMAP on a subset of cells ###
  n_sel <- min(nrow(df), 3e4)
  sel_umap <- sample(nrow(df), n_sel)
  df_um <- get_umap(df, data_scaled, cols, sel_umap)

  ### Classify cells using pretrained logistic regression model ###
  cell_type <- classify_cells_logit(data_scaled, defs_major)

  ### Detect doublets ###
  print(paste("Doublet detection for file", fn, "..."))
  cleanet_res <- detect_doublets(df, cols, cell_type, dir_out, fn, cofactor=cofactor)
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

  ### Refine cell subsets ###
  cell_type <- refine_subsets(data_scaled, cell_type, "T cell", defs_tcell)
  cell_type <- refine_subsets(data_scaled, cell_type, "T cell CD4", defs_cd4_mem)
  cell_type <- refine_subsets(data_scaled, cell_type, "T cell CD8", defs_cd8_mem)
  cell_type <- refine_subsets(data_scaled, cell_type, "Myeloid", defs_myel)
  cell_type <- refine_subsets(data_scaled, cell_type, "Neutrophil", defs_neut)

  backgate_major(df, cell_type, dir_out, fn)

  print(paste("Gating file", fn, "..."))
  gate_detailed_phenos(df, data_scaled, data_transf, cell_type, dir_out, fn)

  event_type[which(event_type=="")] <- cell_type
  ff_clean <- clean_fcs_file(ff, event_type)
  write.FCS(ff_clean, paste0(dir_out, "fcs_clean/", fn, "_cleaned.fcs"))

  df_file <- tibble(event_type=event_type)
  write_csv(df_file, file=paste0(dir_out, "files_labeled/", fn, ".csv"), progress=FALSE)

  ### Compute density estimates to be used later for QC ###
  df_kdes <- estimate_distributions(cell_type, data_transf, fn)
  write_csv(df_kdes, file=paste0(dir_out, "/kdes_for_qc/kdes_", fn, ".csv"), progress=FALSE)

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
  for (n in c("CD45RA CD27", "CD45 CD66b",
              "CD4 CD8a", "CD3 CD19", "CD11c CD14", "CD3 CD56",
              "CD16 CD66b", "CD123 CD294", "CD3 TCRgd", "CD14 CD38")) {
    dir.create(paste0(dir_out, "/backgating/", n), showWarnings = FALSE)
  }
}


#### de novo AML : inv16, t(8:21), M0, M1, M2, M4, M5
#### anything with MLD, Secondary

# IS AML different for HD?
# Is AML with MLD different from AML?
# Does this change with therapy?
# TP53, FLT3, IDH1/2: compare molecular subtypes after comparison of MLD






