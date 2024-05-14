

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
    str_remove("_Processed.fcs")
  message(fn)

  df <- path %>%
    read_data() %>%
    as_tibble() %>%
    pregate_data(fn, dir_out)

  df <- df[sample(nrow(df), min(nrow(df),1e5)),]
  x_full <- df %>% as.matrix() %>% scale_data()
  x <- x_full[,cols]
  x <- df %>% select(all_of(cols)) %>% as.matrix() %>% scale_data()

  n_sel <- min(nrow(df), 3e4)
  sel_umap <- sample(nrow(df), n_sel)
  df_um <- get_umap(df, x, sel_umap)

  ### clustering and annotating clusters
  clustering_raw <- run_fastpg(x, resolution = 1, n_threads = 1)
  clustering <- merge_clusters_c(df, cols, clustering_raw, min_bhatt = 0.2)

  centroids <- get_medians(df, clustering)
  mat <- data.matrix(df) %>% scale_data()
  centroids_sc <- get_medians(mat, clustering)
  labels <- label_clusters_score_opt(centroids_sc, defs_major) %>% make.unique()


  rownames(centroids) <- rownames(centroids_sc) <- levels(clustering) <- labels

  cl_tmp <- clustering
  cell_type <- as.character(cl_tmp) %>% str_split(fixed(".")) %>% sapply("[",1)

  ######### new doublet detection #########

  cell_type <- detect_doublets_apr(df, cols, cell_type)#, thresh = 4)
  ##########################################

  mat <- as.matrix(df) %>% scale_data()
  ypred <- label_clusters_score_opt(mat, defs_major, mdipa_main = FALSE, return_ypred = TRUE)
  pred <- colnames(ypred)[unname(apply(ypred, 1, which.max))]
  conf <- apply(ypred, 1, max)

  cell_type <- update_pred(cell_type, pred, conf)


  ### visualization
  df_clust <- df_um %>%
    mutate(clust = cl_tmp[sel_umap]) %>%
    mutate(ct = cell_type[sel_umap]) %>%
    mutate(event = case_when(grepl("agg", ct) ~ "doublet",
                             grepl("debris", ct) ~ "debris",
                             TRUE ~ "single cell"))
  # write_csv(df_clust, file=paste0(dir_out, "tmp/df_clust_", fn, ".csv"), progress=FALSE)

  p_major <- plot_umap_major(df_clust, fn)
  ggsave(p_major, filename = paste0(dir_out, "umap_major/", fn, ".png"),
         width=12, height=10)

  # p_debris <- plot_umap_debris(df_clust, fn)
  # ggsave(p_debris, filename = paste0(dir_out, "umap_debris/", fn, ".png"),
  #        width=16, height=10)

  p_dna <- plot_dna_cd45(df_clust, fn)
  ggsave(p_dna, filename = paste0(dir_out, "DNA_CD45/", fn, ".png"),
         width=9, height=7)

  #### monocytes and DC

  is_myel <- cell_type %in% c("myeloid", "pdc", "basophil")

  if (length(which(is_myel)) > 20) {
    df_mono <- df %>% filter(is_myel)

    cols_mono <- c("CD45", "CD123", "CD294", "CD11c", "CD16", "CD14",
                   "CD66b", "CD38", "HLA-DR", "CD45RA", "CD45RO")
    cols_mono_strict <- c("CD11c", "CD14", "CD38", "CD123", "CD294",
                          "HLA-DR", "CD45RA", "CD45RO")

    x_mono <- x_full[which(is_myel),cols_mono]

    set.seed(0)
    sel_mono <- sample(nrow(x_mono), min(nrow(x_mono), 3e4))
    df_um_mono <- get_umap(df_mono, x_mono, sel_mono, "CD14", "CD294")

    clustering_mono_raw <- run_fastpg(x_mono, resolution = 1, n_threads = 1)
    clustering_mono <- merge_clusters_c(df_mono, cols=cols_mono_strict,
                                        clustering_mono_raw, min_bhatt = 0.35)

    mono_cent <- get_medians(x_mono, clustering_mono)
    labels <- label_clusters_score_opt(mono_cent, defs_myel)
    rownames(mono_cent) <- levels(clustering_mono) <- labels

    df_um_mono <- df_um_mono %>%
      mutate(clust = clustering_mono[sel_mono])
    # write_csv(df_um_mono, file=paste0(dir_out, "tmp/df_mono_", fn, ".csv"), progress=FALSE)

    p_mono <- plot_umap_mono(df_um_mono, fn)
    ggsave(p_mono, filename=paste0(dir_out, "umap_mono/", fn, ".png"),
           width=12, height=10)

    cell_type <- update_clustering(cell_type, clustering_mono, which(is_myel))
  } else {
    p_mono <- NULL
  }

  #### t cells

  is_tcell <- cell_type=="tcell"

  if(length(which(is_tcell)) > 20) {
    df_tcell <- df %>% filter(is_tcell)
    cols_tcell <- c("CD3", "CD45", "CD4", "CD8a", "TCRgd")
    cols_tcell_main <- c("CD4", "CD8a", "TCRgd")

    x_tcell <- x_full[which(is_tcell),cols_tcell]

    set.seed(0)
    sel_tcell <- sample(nrow(x_tcell), min(nrow(x_tcell), 3e4))
    df_um_tcell <- get_umap(df_tcell, x_tcell, sel_tcell, "CD4", "TCRgd")

    clustering_tcell_raw <- run_fastpg(x_tcell, resolution = 0.5, n_threads = 1)
    clustering_tcell <- merge_clusters_c(df_tcell, cols_tcell_main, clustering_tcell_raw,
                                         min_bhatt = 0.5)

    tcell_cent <- get_medians(x_tcell, clustering_tcell)
    labels <- label_clusters_score_opt(tcell_cent, defs_tcell, mdipa_main=FALSE)
    rownames(tcell_cent) <- levels(clustering_tcell) <- labels

    df_um_tcell <- df_um_tcell %>%
      mutate(clust = clustering_tcell[sel_tcell])
    # write_csv(df_um_tcell, file=paste0(dir_out, "tmp/df_tcell_", fn, ".csv"), progress=FALSE)

    p_tcell <- plot_umap_tcell(df_um_tcell, fn)
    ggsave(p_tcell, filename=paste0(dir_out, "umap_tcell/", fn, ".png"), width=12, height=10)

    cell_type <- update_clustering(cell_type, clustering_tcell, which(is_tcell))
  } else {
    p_tcell <- NULL
  }

  df_file <- tibble(cell_type=cell_type)
  write_csv(df_file, file=paste0(dir_out, "files_labeled/", fn, ".csv"), progress=FALSE)

  ### gating secondary cell types
  gating_stuff(df, cell_type, dir_out, fn)

  return(NULL)

  # ### density estimates by cell type
  # ### to be used later for QC
  # valid_cell_types <- c("neutrophil", "eosinophil", "basophil", "bcell", "pdc",
  #                       "monocyte_classical", "tcell_cd4", "tcell_cd8", "tcell_gd")
  # df_kdes <- estimate_distributions(cell_type, df, fn, cols, valid_cell_types)
  # write_csv(df_kdes, file=paste0(dir_out, "/kdes_for_qc/", fn, ".csv"), progress=FALSE)

}


#' @title Create output directory structure
#' @description Create subdirectories for saving figures and csv output.
#' @param dir_out Base output path.
#' @returns Nothing.
#' @export
create_dirs <- function(dir_out) {
  dir.create(paste0(dir_out, "/files_labeled"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/feat_major"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/feat_adaptive"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/kdes_for_qc"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/cleanup_gates"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/cleanup_stats"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/umap_major"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/umap_tcell"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/umap_mono"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/DNA_CD45"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/gating"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/bcell_mem"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/cd4_func"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/mait_nkt"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/nk_late_early"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/tcell_act"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/tcell_em_cm"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/tcell_mem"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/tfh"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/gating/treg"), showWarnings = FALSE)

  dir.create(paste0(dir_out, "/gating/kdes_thresh"), showWarnings = FALSE)
}



