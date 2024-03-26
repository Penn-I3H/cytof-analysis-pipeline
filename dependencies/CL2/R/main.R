

#' @title Analyze a single CyTOF file
#' @description Read and parse file, cluster and label main cell types, detect
#' debris and doublets, gate downstream T and B cell subpopulations.
#' @param file FCS file name.
#' @param dir_in A path to read fcs files from.
#' @param dir_out A path to write output.
#' @param cols Columns to use in analysis.
#' @param defs Model definitions for main cell types.
#' @param defs_tcell Model definitions for T cell subtypes.
#' @param defs_myeloid Model definitions for monocytes and DCs.
#' @returns Nothing; all output written to file.
#' @export
analyze_cytof_file <- function(file, dir_in, dir_out, cols,
                               defs, defs_tcell, defs_myeloid) {
  set.seed(0)
  path <- paste0(dir_in, file)
  fn <- file %>%
    str_remove("_Normalized.fcs") %>%
    str_remove("_normalized.fcs") %>%
    str_remove("_Processed.fcs")
  message(fn)

  df <- path %>%
    read_data() %>%
    as_tibble()

  df <- df[sample(nrow(df), min(nrow(df),1e5)),]
  x_full <- df %>% as.matrix() %>% scale_data()
  x <- x_full[,cols]
  x <- df %>% select(all_of(cols)) %>% as.matrix() %>% scale_data()

  n_sel <- min(nrow(df), 1e4)
  sel_umap <- sample(nrow(df), n_sel)
  df_um <- get_umap(df, x, sel_umap)

  ### clustering and annotating clusters
  clustering_raw <- run_fastpg(x, resolution = 1, n_threads = 1)
  clustering <- merge_clusters_c(df, cols, clustering_raw, min_bhatt = 0.2)

  gr <- detect_doublets(df, cols, clustering)

  centroids <- get_medians(df, clustering)
  mat <- data.matrix(df) %>% scale_data()
  centroids_sc <- get_medians(mat, clustering)

  labels <- get_labels_with_doublets_mar(clustering, df, gr, defs, centroids_sc)

  rownames(centroids) <- rownames(centroids_sc) <- levels(clustering) <- labels

  cl_tmp <- separate_doublets_feb(clustering, df, cols)
  cell_type <- as.character(cl_tmp) %>% str_split(fixed(".")) %>% sapply("[",1)

  mat <- as.matrix(df) %>% scale_data()
  ypred <- label_clusters_score_opt(mat, defs, mdipa_main = FALSE, return_ypred = TRUE)
  pred <- colnames(ypred)[unname(apply(ypred, 1, which.max))]
  conf <- apply(ypred, 1, max)

  cell_type <- update_pred(cell_type, pred, conf)


  ### visualization
  df_clust <- df_um %>%
    mutate(clust = cl_tmp[sel_umap]) %>%
    mutate(ct = cell_type[sel_umap]) %>%
    mutate(event = case_when(grepl("agg", clust) ~ "doublet",
                             grepl("debris", clust) ~ "debris",
                             TRUE ~ "single cell"))

  p <- plot_umap_major(df_clust, fn)
  ggsave(p, filename = paste0(dir_out, "umap_major/", fn, ".png"),
         width=12, height=10)

  p <- plot_dna_cd45(df_clust, fn)
  ggsave(p, filename = paste0(dir_out, "DNA_CD45/", fn, ".png"),
         width=9, height=7)

  df_file <- tibble(cell_type=cell_type,
                    pred = if_else(grepl("agg", cell_type), cell_type, pred))
  write_csv(df_file, file=paste0(dir_out, "/files_labeled/", fn, ".csv"), progress=FALSE)

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
    sel_mono <- sample(nrow(x_mono), min(nrow(x_mono), 1e4))
    df_um_mono <- get_umap(df_mono, x_mono, sel_mono, "CD14", "CD294")

    clustering_mono_raw <- run_fastpg(x_mono, resolution = 1, n_threads = 1)
    clustering_mono <- merge_clusters_c(df_mono, cols=cols_mono_strict,
                                        clustering_mono_raw, min_bhatt = 0.35)

    mono_cent <- get_medians(x_mono, clustering_mono)
    labels <- label_clusters_score_opt(mono_cent, defs_myel)
    rownames(mono_cent) <- levels(clustering_mono) <- labels

    df_um_mono <- df_um_mono %>%
      mutate(clust = clustering_mono[sel_mono])

    p <- plot_umap_mono(df_um_mono, fn)
    ggsave(p, filename=paste0(dir_out, "umap_myel/", fn, ".png"),
           width=12, height=10)

    cell_type <- update_clustering(cell_type, clustering_mono, which(is_myel))
  }


  #### t cells

  is_tcell <- cell_type=="tcell"

  if(length(which(is_tcell)) > 20) {
    df_tcell <- df %>% filter(is_tcell)
    cols_tcell <- c("CD3", "CD45", "CD4", "CD8a", "TCRgd")
    cols_tcell_main <- c("CD4", "CD8a", "TCRgd")

    x_tcell <- x_full[which(is_tcell),cols_tcell]

    set.seed(0)
    sel_tcell <- sample(nrow(x_tcell), min(nrow(x_tcell), 1e4))
    df_um_tcell <- get_umap(df_tcell, x_tcell, sel_tcell, "CD4", "TCRgd")

    clustering_tcell_raw <- run_fastpg(x_tcell, resolution = 0.5, n_threads = 1)
    clustering_tcell <- merge_clusters_c(df_tcell, cols_tcell_main, clustering_tcell_raw,
                                         min_bhatt = 0.5)

    tcell_cent <- get_medians(x_tcell, clustering_tcell)
    labels <- label_clusters_score_opt(tcell_cent, defs_tcell, mdipa_main=FALSE)
    rownames(tcell_cent) <- levels(clustering_tcell) <- labels

    df_um_tcell <- df_um_tcell %>%
      mutate(clust = clustering_tcell[sel_tcell])

    p <- plot_umap_tcell(df_um_tcell, fn)
    ggsave(p, filename=paste0(dir_out, "umap_tcell/", fn, ".png"), width=12, height=10)

    cell_type <- update_clustering(cell_type, clustering_tcell, which(is_tcell))
  }

  gating_stuff(df, cell_type, dir_out, fn)
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

  dir.create(paste0(dir_out, "/umap_major"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/umap_tcell"), showWarnings = FALSE)
  dir.create(paste0(dir_out, "/umap_myel"), showWarnings = FALSE)
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



