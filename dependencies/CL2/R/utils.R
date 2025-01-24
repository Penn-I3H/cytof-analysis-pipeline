
read_data <- function(ff, transform=TRUE) {
  pdata <- ff@parameters@data

  # remove unused channels
  non_metal <- grep("Time|Center|Offset|Width|Residual|length", pdata$name)
  markers <- setdiff(grep("_", pdata$desc), non_metal)

  if (transform) {
    data <- cbind(asinh(ff@exprs[,markers]/5),
                  ff@exprs[,non_metal])
    data[,c("Center", "Offset", "Width", "Residual")] <- asinh(data[,c("Center", "Offset", "Width", "Residual")]/5)
  } else {
    data <- ff@exprs[,c(markers, non_metal)]
  }

  # use protein rather than isotope for channel names
  nice_names <- pdata$desc[markers] %>%
    str_split("_") %>%
    sapply(function(x) x[2])

  bead_ch <- which(nice_names=="Bead")
  nice_names[bead_ch] <- pdata$desc[markers[bead_ch]]
  colnames(data)[seq_along(markers)] <- nice_names

  return(data)
}


plot_cleanup_gate <- function(df, cutoffs, marker, perc) {
  t0 <- min(df$Time)
  t1 <- max(df$Time)
  tdiff <- t1-t0
  t0 <- t0 - 0.05*tdiff
  t1 <- t1 + 0.05*tdiff
  m0 <- cutoffs[[marker]][1]
  m1 <- cutoffs[[marker]][2]

  df_gate <- tibble(Time = c(t0, t0, t1, t1, t0),
                 Marker = c(m0, m1, m1, m0, m0))
  df_text <- tibble(label = paste0(round(100*perc,2), "%"),
                    t1=t1, m1=m1)

  plot_bin2d(df, channels = c("Time", marker), scales = c("linear", "asinh")) +
    geom_path(data=df_gate, aes(x=Time, y=Marker), color="red") +
    geom_label(data=df_text, aes(x=(t1+t0)/2, y=1.05*m1, label=label), alpha=0.7) +
    theme(axis.text.x = element_text(angle=45))
}


apply_cleanup_gate <- function(df, event_type, cutoffs, marker, label) {
  filt <- which(event_type=="")
  dat <- df[[marker]][filt]

  idx_out <- which(dat < cutoffs[[marker]][1] | dat > cutoffs[[marker]][2])
  event_type[filt[idx_out]] <- label
  return(event_type)
}

asinh_trans <- trans_new(
  name = "asinh",
  transform = function(x) asinh(x / 5),
  inverse = function(x) sinh(x) * 5
)


asinh_breaks_major <- c(0,10^seq(1,5))
asinh_breaks_minor <- as.numeric(outer(seq(2,9),10^seq(1,4)))


set_plot_scale <- function(p, scales) {

  if (scales[1]=="asinh")
    p <- p + scale_x_continuous(transform = asinh_trans,
                                breaks = asinh_breaks_major,
                                labels = asinh_breaks_major,
                                minor_breaks = asinh_breaks_minor)

  if (scales[2]=="asinh")
    p <- p + scale_y_continuous(transform = asinh_trans,
                                breaks = asinh_breaks_major,
                                labels = asinh_breaks_major,
                                minor_breaks = asinh_breaks_minor)
  return(p)
}


plot_bin2d <- function(df, channels, scales=c("linear", "linear")) {
  p <- ggplot(df, aes(x=.data[[channels[1]]], y=.data[[channels[2]]])) +
    geom_bin_2d(bins=100) +
    scale_fill_gradientn(colours = hsv(h = seq(0.6667,0, length.out = 11))) +
    guides(fill="none") +
    theme_bw(base_size=14)

  return(set_plot_scale(p, scales))
}


pregate_data <- function(df, fn, dir_out, plot=FALSE) {

  cutoffs <- list("140Ce_Bead" = c(0,500),
                  "Residual" = c(20,400),
                  "Center" = c(300,4000),
                  "Offset" = c(30,400),
                  "Width" = c(30,300),
                  "Live" = c(0,80))
  markers <- names(cutoffs)
  labels <- c("Bead", "Residual", "Center", "Offset", "Width", "Dead")

  event_type <- character(nrow(df))

  for (i in seq_along(labels)) {
    event_type <- apply_cleanup_gate(df, event_type, cutoffs, markers[i], labels[i])
  }

  df1 <- df[which(!event_type %in% labels[seq(1)]),]
  df2 <- df[which(!event_type %in% labels[seq(2)]),]
  df3 <- df[which(!event_type %in% labels[seq(3)]),]
  df4 <- df[which(!event_type %in% labels[seq(4)]),]
  df5 <- df[which(!event_type %in% labels[seq(5)]),]
  df6 <- df[which(!event_type %in% labels[seq(6)]),]

  if(plot) {
    p1 <- plot_cleanup_gate(df, cutoffs, markers[1], nrow(df1)/nrow(df))
    p2 <- plot_cleanup_gate(df1, cutoffs, markers[2], nrow(df2)/nrow(df1))
    p3 <- plot_cleanup_gate(df2, cutoffs, markers[3], nrow(df3)/nrow(df2))
    p4 <- plot_cleanup_gate(df3, cutoffs, markers[4], nrow(df4)/nrow(df3))
    p5 <- plot_cleanup_gate(df4, cutoffs, markers[5], nrow(df5)/nrow(df4))
    p6 <- plot_cleanup_gate(df5, cutoffs, markers[6], nrow(df6)/nrow(df5))
    p <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol=3)
    ggsave(p, filename = paste0(dir_out, "cleanup_gates/", fn, ".png"),
           width=9, height=6)
  }

  df_stats <- tibble(file = fn,
                     n_events = nrow(df),
                     n_bead_gate = nrow(df1),
                     n_residual_gate = nrow(df2),
                     n_center_gate = nrow(df3),
                     n_offset_gate = nrow(df4),
                     n_width_gate = nrow(df5),
                     n_live_gate = nrow(df6))

  return(list(event_type=event_type, df_stats=df_stats))
}


write_stats_file <- function(df_stats, cell_type, dir_out, fn) {
  file_stats <- paste0(dir_out, "cleanup_stats/cleanup_stats_", fn, ".csv")
  df_stats <- df_stats %>%
    mutate(n_not_debris = length(which(cell_type!="Debris")),
           n_cd45_single = length(which(!grepl("Debris|Doublet|_", cell_type))))

  write_csv(df_stats, file = file_stats, progress=FALSE)
}


scale_data <- function(mat) {
  mat_sc <- apply(mat, 2, function(col) {
    col / quantile(col,0.999)
  })
  return(mat_sc)
}


get_umap <- function(df, data_scaled, cols, sel_umap, dir1="DNA1", dir2="CD3") {
  x <- data_scaled[,cols]

  um0 <- umap(x[sel_umap,])
  v <- x[sel_umap,dir1]
  cos1 <- sum(v*um0[,1]) / ( sqrt(sum(v*v)) * sqrt(sum(um0[,1]*um0[,1])) )
  cos2 <- sum(v*um0[,2]) / ( sqrt(sum(v*v)) * sqrt(sum(um0[,2]*um0[,2])) )
  ang <- acos(cos1) * sign(cos2)
  rot <- matrix( c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow=2)
  um <- um0 %*% rot
  if (cor(x[sel_umap,dir2], um[,2]) < 0)
    um[,2] <- -um[,2]

  df_um <- df[sel_umap,] %>% mutate(umap1 = um[,1], umap2 = um[,2])
  return(df_um)
}


plot_umap_channel <- function(df, channel, base_size=8) {
  ggplot(df, aes(x=umap1, y=umap2, color=.data[[channel]])) +
    geom_point(size=0.5, shape=1, alpha=0.5) +
    scale_color_viridis_c(trans = asinh_trans,
                          breaks = asinh_breaks_major,
                          labels = asinh_breaks_major) +
    theme_bw(base_size=base_size) +
    theme(axis.title = element_blank())
}


plot_umap_ct <- function(df, channel, fn, pal=NULL, base_size=8) {
  p <- ggplot(df, aes(x=umap1, y=umap2, color=.data[[channel]])) +
    geom_point(size=0.5, shape=1, alpha=0.5) +
    ggtitle(fn) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape=19))) +
    theme_bw(base_size=base_size) +
    theme(axis.title = element_blank())

  if (!is.null(pal))
    p <- p + scale_color_manual(values=pal, name="Cell type")

  return(p)
}


plot_umap_major_only <- function(df, fn) {
  pal_paired <- brewer.pal(12, "Paired")
  names(pal_paired) <- c("Neutrophil", "Eosinophil",
                         "B cell", "Plasmablast",
                         "Basophil", "pDC",
                         "T cell", "NK cell",
                         "Doublet", "Myeloid",
                         "Debris", "Uncertain")
  pal_paired["Doublet"] <- "gray30"

  plot_umap_ct(df %>% mutate(ct = if_else(grepl("_", ct), "Doublet", ct)),
               "ct", fn, pal_paired, base_size=16)
}


plot_umap_markers <- function(df, fn) {

  p1 <- plot_umap_channel(df, "CD45", base_size=12)
  p2 <- plot_umap_channel(df, "CD3", base_size=12)
  p3 <- plot_umap_channel(df, "CD19", base_size=12)
  p4 <- plot_umap_channel(df, "CD66b", base_size=12)
  p5 <- plot_umap_channel(df, "CD123", base_size=12)
  p6 <- plot_umap_channel(df, "CD294", base_size=12)
  p7 <- plot_umap_channel(df, "CD56", base_size=12)
  p8 <- plot_umap_channel(df, "CD11c", base_size=12)
  p9 <- plot_umap_channel(df, "DNA1", base_size=12)

  p <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3) +
    plot_annotation(title=fn)

  return(p)
}


plot_dna_cd45 <- function(df, fn) {
  df <- df %>%
    filter(ct!="Uncertain") %>%
    mutate(event = case_when(grepl("_", ct) ~ "Doublet",
                             grepl("Debris", ct) ~ "Debris",
                             TRUE ~ "Single cell"))

  dark2 <- brewer.pal(8, "Dark2")
  pal_dark2 <- c("Debris"="#FFFF99", "Single cell"=dark2[1], "Doublet"="gray30")

  p <- ggplot(df, aes(x=DNA1, y=CD45, color=event)) +
    geom_point(size=0.5, shape=1, alpha=0.5) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape=19))) +
    scale_color_manual(values=pal_dark2, name="Event type") +
    ggtitle(fn) +
    theme_bw(base_size=16)

  return(set_plot_scale(p, c("asinh", "asinh")))
}


classify_cells_logit <- function(data_scaled, defs_major) {
  pred_probs <- predict_cell_type(data_scaled, defs_major, return_probs = TRUE)
  uniq <- unique(colnames(pred_probs))
  pred_probs <- lapply(uniq, function(ct) {
    apply(pred_probs[,which(colnames(pred_probs)==ct),drop=FALSE], 1, sum)
  }) %>% do.call(what=cbind)
  colnames(pred_probs) <- uniq
  # return(pred_probs)
  cell_type <- colnames(pred_probs)[unname(apply(pred_probs, 1, which.max))]
  confidence <- apply(pred_probs, 1, function(row) {
    sor <- sort(row, decreasing=TRUE)
    tmp <- intersect( names(sor)[seq(2)], c("Debris", "Neutrophil") )
    if (length(tmp)==2)
      return(sor[1]-sor[2]+0.1)
    return(sor[1]-sor[2])
  })
  uncertain_cells <- which(confidence < 0.1 & cell_type %in% c("pDC", "Plasmablast") |
                             confidence < 0.05)

  cell_type[uncertain_cells] <- "Uncertain"
  return(cell_type)
}


predict_cell_type <- function(data, defs, return_probs=FALSE) {
  coeff <- defs %>% select(-Phenotype)

  coeff <- as.matrix(coeff)
  rownames(coeff) <- defs$Phenotype
  mat <- coeff[,-1,drop=FALSE]
  intercept <- coeff[,1]

  probs_not_norm <- t(mat %*% t(data[,colnames(mat),drop=FALSE]) + intercept)

  probs <- exp(probs_not_norm) %>%
    apply(1, function(row) row/sum(row)) %>% t()

  if (return_probs)
    return(probs)

  lab <- colnames(probs)[unname(apply(probs, 1, which.max))]
  return(lab)
}


kde_single_mat <- function(mat, name, range=c(-3,6)) {
  channels <- colnames(mat)

  lapply(channels, function(channel) {
    kde <- bkde(mat[,channel], range.x = range)
    tabular <- tibble(protein_expression = kde$x,
                      density_estimate = kde$y/max(kde$y),
                      channel = channel,
                      sample = name)
  }) %>% do.call(what=rbind)
}


plot_kdes <- function(df_kdes, sel_vals, sel_col="celltype") {
  df_plot <- df_kdes %>%
    filter(.data[[sel_col]] %in% sel_vals)

  ggplot(df_plot, aes(x=intensity, y=density, color=.data[[sel_col]])) +
    geom_path() +
    facet_wrap(~channel, scales="free") +
    theme_bw()
}


get_kdes <- function(x, cols, cell_type, unq) {

  df_kdes <- lapply(cols, function(ch) {
    dat <- x[,ch]
    m <- min(dat)
    M <- max(dat)
    lapply(unq, function(ct) {
      dat_ct <- dat[which(cell_type==ct)]

      if (length(dat_ct) < 3)
        return(NULL)

      if (length(dat_ct) < 10)
        dat <- jitter(rep(dat_ct,10))
      else {
        if (sd(dat_ct)==0)
          dat_ct <- jitter(dat_ct)
      }

      kde <- bkde(dat_ct, range.x = c(m,M), gridsize=101L)
      return(tibble(celltype=ct,
                    channel=ch,
                    intensity=kde$x,
                    density=kde$y / sum(kde$y)))
    }) %>% do.call(what=rbind)
  }) %>% do.call(what=rbind) %>%
    mutate(density = pmax(density,0))

  return(df_kdes)
}


detect_doublets <- function(df, cols, cell_type, dir_out, fn, thresh=5, cofactor=cofactor) {

  cleanet_res <- cleanet(df, setdiff(cols, "Event_length"),
                         cofactor=cofactor, is_debris=cell_type=="Debris")

  singlet_clas <- cell_type[which(cleanet_res$status!="Doublet")]
  doublet_clas <- classify_doublets(cleanet_res, singlet_clas)
  df_exp_obs <- compare_doublets_exp_obs(doublet_clas, singlet_clas, cleanet_res)

  cell_type[which(cleanet_res$status=="Doublet")] <- doublet_clas
  cleanet_res$cell_type <- cell_type

  write_csv(df_exp_obs %>% arrange(-Observed),
            file=paste0(dir_out, "/doublet_csv/doublet_", fn, ".csv"))

  p <- plot_expected_observed(df_exp_obs, fn)
  ggsave(p, filename=paste0(dir_out, "/doublet_fig/doublet_", fn, ".png"),
         width=6.75, height=7)

  return(cleanet_res)
}


plot_expected_observed <- function(df_exp_obs, fn) {
  df_text <- df_exp_obs %>%
    filter(Observed > 0.03) %>%
    mutate(hjust = if_else(Expected > 0.1, 1, 0))

  ggplot(df_exp_obs, aes(x=Expected, y=Observed)) +
    geom_point() +
    geom_text(data=df_text, aes(label=Type, hjust=hjust), vjust=1) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    ggtitle(paste0("Doublet distribution for ", fn)) +
    theme_bw(base_size=16) +
    theme(plot.title = element_text(size=16))
}


backgate_major <- function(df, cell_type, dir_out, fn) {
  channels <- c("CD4", "CD8a")
  cells <- which(grepl("T cell", cell_type) & !grepl("Doublet|_", cell_type))
  cts <- c("T cell CD4", "T cell CD8", "T cell DN", "T cell gd")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD3", "TCRgd")
  cells <- which(grepl("T cell", cell_type) & !grepl("Doublet|_", cell_type))
  cts <- c("T cell CD4", "T cell CD8", "T cell DN", "T cell gd")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD3", "CD19")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("T cell", "B cell")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD45", "CD66b")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("Neutrophil CD16hi", "Neutrophil CD16lo", "Eosinophil")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD3", "CD56")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("T cell", "NK cell")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD11c", "CD14")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("Mono|mDC", "Neutrophil CD16hi")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD14", "CD38")
  cells <- which(grepl("Mono|mDC", cell_type) & !grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("Monocyte Classical", "Monocyte Nonclassical", "mDC")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD16", "CD66b")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("Eosinophil", "Neutrophil CD16hi", "Neutrophil CD16lo")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD123", "CD294")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain", cell_type))
  cts <- c("Eosinophil", "Basophil", "pDC")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD38", "CD27")
  cells <- which(!grepl("Debris|Doublet|_|Uncertain|Neutrophil|T cell|Mono", cell_type))
  cts <- c("Plasmablast", "B cell")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD45RA", "CD27")
  cells <- which(grepl("T cell CD4", cell_type) & !grepl("Doublet|_", cell_type))
  cts <- c("T cell CD4 Naive", "T cell CD4 Mem")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)

  channels <- c("CD45RA", "CD27")
  cells <- which(grepl("T cell CD8", cell_type) & !grepl("Doublet|_", cell_type))
  cts <- c("T cell CD8 Naive", "T cell CD8 Mem")
  backgate_one_plot(df, cell_type, channels, cells, cts, dir_out, fn)
}


backgate_one_plot <- function(df, cell_type, channels, cells, cts, dir_out, fn) {
  label_array <- cell_type[cells]
  pattern <- paste(cts, collapse="|")
  label_array[which(!grepl(pattern, label_array))] <- "Other"

  for (ct in cts) {
    label_array[grep(ct, label_array)] <- ct
  }

  label_name <- paste(channels, collapse=" ")

  backgate(df[cells,channels], label_array, label_name,
           cts, channels, dir_out, fn)
}


backgate <- function(df, label_array, label_name, cts, channels, dir_out, fn) {
  sel <- sample(nrow(df), min(nrow(df), 1e4))
  df_plot <- df[sel,]

  pal <- gg_palette(length(cts)+1)
  names(pal) <- c("Other", cts)

  p <- ggplot(df_plot, aes(x=.data[[channels[1]]],
                           y=.data[[channels[2]]],
                           color=label_array[sel])) +
    geom_point(size=0.3, alpha=0.3) +
    ggtitle(fn, subtitle = label_name) +
    scale_color_manual(name="Cell type", values=pal) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape=19))) +
    theme_bw()
  p <- set_plot_scale(p, c("asinh", "asinh"))
  ggsave(p, filename=paste0(dir_out, "backgating/", label_name, "/",
                            paste(cts, collapse=" "), "_", fn, ".png"),
         width=7, height=6)
}


gg_palette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[seq(n)]
}


get_thresh_opt <- function(df, cell_type, unq, thresh, channel, pop1, pop2,
                           w = c(1,1), search_func = search_opt) {

  pattern <- paste0(pop1,"|",pop2)

  if (!any(grepl(pattern, unq)))
    return(thresh[channel])

  df_pop <- df %>% select(channel) %>%
    mutate(ct=cell_type) %>% filter(grepl(pattern, ct) & !grepl("_|Doublet",ct))
  class1_idx <- grep(pop1, df_pop$ct)
  class2_idx <- grep(pop2, df_pop$ct)

  if (length(class1_idx) < 5 | length(class2_idx) < 5)
    return(thresh[channel])

  opt <-  search_func(x=asinh(df_pop[[channel]]/5),
                     class1_idx = class1_idx,
                     class2_idx = class2_idx,
                     w=w)
  return(opt)
}


get_thresh_default <- function(df, mat) {
  thresh <- setNames(rep(0, ncol(mat)), names(df))
  thresh_quant <- apply(mat,2,function(col) quantile(col, 0.999))
  ch_020 <- c("CD25")
  # ch_025 <- c("CD16")
  ch_050 <- c("CD194", "CD161")
  thresh[ch_020] <- 0.20 * thresh_quant[ch_020]
  # thresh[ch_025] <- 0.25 * thresh_quant[ch_025]
  thresh[ch_050] <- 0.50 * thresh_quant[ch_050]
  thresh[which(thresh<=0)] <- 0.35 * thresh_quant[which(thresh<=0)]
  return(thresh)
}


update_thresh_nk <- function(df, cell_type, unq, thresh) {
  thresh["CD57"] <- get_thresh_opt(df, cell_type, unq, thresh, "CD57",
                                   "T cell CD4 Naive|T cell CD8 Naive", "NK cell")
  thresh["CD56"] <- get_thresh_opt(df, cell_type, unq, thresh, "CD56",
                                   "T cell CD4 Naive|T cell CD8 Naive", "NK cell")
  thresh["CD57"] <- max(thresh["CD57"], asinh(15/5))
  return(thresh)
}


update_thresh_cd4_func <- function(df, data_scaled, cell_type, unq, thresh) {
  ct_func <- refine_subsets(data_scaled, cell_type, "T cell CD4 Mem",
                            defs_cd4_func, only_pred=TRUE)
  df_func <- df %>% filter(cell_type=="T cell CD4 Mem")

  thresh["CD183"] <- get_thresh_opt(df_func, ct_func, unique(ct_func), thresh, "CD183",
                                    "T cell CD4 Other", "T cell CD4 Th1")
  thresh["CD185"] <- get_thresh_opt(df_func, ct_func, unique(ct_func), thresh, "CD185",
                                    "T cell CD4 Other", "T cell CD4 Tfh")

  thresh["CD196"] <- get_thresh_opt(df, cell_type, unq, thresh, "CD196",
                                    "T cell CD4 Naive|T cell CD8 Naive", "T cell CD4 Mem")

  thresh["CD185"] <- max(thresh["CD185"], asinh(10/5))
  thresh["CD196"] <- max(thresh["CD196"], asinh(5/5))
  thresh["CD196"] <- min(thresh["CD196"], 0.5*max(asinh(df_func[["CD196"]]/5)))

  sel_cd194 <- grep("T cell CD4 Naive|T cell CD4 Mem", cell_type)
  df_cd194 <- df[sel_cd194,] %>%
    mutate(ct = case_when(cell_type[sel_cd194]=="T cell CD4 Naive" ~ "T cell CD4 Naive",
                          CD185 < 5*sinh(thresh["CD185"]) &
                            CD183 < 5*sinh(thresh["CD183"]) ~ "Candidates",
                          TRUE ~ "Other"))
  thresh["CD194"] <- get_thresh_opt(df_cd194, df_cd194$ct, unique(df_cd194$ct),
                                    thresh, "CD194",
                                    "T cell CD4 Naive", "Candidates")
  return(thresh)
}


gate_detailed_phenos <- function(df, data_scaled, data_transf, cell_type, dir_out, fn) {
  unq <- unique(cell_type[which(!grepl("Debris|Doublet|_|Uncertain",cell_type))])
  df_kdes <- get_kdes(data_transf, names(df), cell_type, unq)

  thresh <- get_thresh_default(df, data_transf)
  thresh <- update_thresh_nk(df, cell_type, unq, thresh)

  if ("T cell CD4 Mem" %in% unq) {
    thresh <- update_thresh_cd4_func(df, data_scaled, cell_type, unq, thresh)
  }

  p <- get_thresholds_plot(df_kdes, thresh)
  ggsave(p, width=12, height=10,
         filename=paste0(dir_out, "gating/thresholds/", fn, ".png"))

  df_feat_maj <- get_feat_major(cell_type, unq, fn)
  write_csv(df_feat_maj,
            file=paste0(dir_out, "/feat_major/feat_major_", fn, ".csv"))

  unq_no_mem <- unq %>% str_remove(" Naive") %>% str_remove(" Mem")
  cts <- intersect(unq_no_mem, gate_hierarchy_full$Parent)

  df_feat_adaptive <- lapply(cts, function(ct) {
    if (!ct %in% unq_no_mem) {
      message(paste(ct, "missing!"))
      return(NULL)
    }

    apply_gate_hierarchy(df, cell_type, thresh, dir_out, fn, ct)
  }) %>% do.call(what=rbind) %>%
    pivot_wider(names_from="feature", values_from="frac") %>%
    mutate(file=fn) %>%
    relocate(file)

  df_feat_gated <- df_feat_adaptive %>% select(!contains("out of"))
  write_csv(df_feat_gated,
            file=paste0(dir_out, "/feat_gated/feat_gated_", fn, ".csv"))

  df_feat_comb <- df_feat_adaptive %>% select(matches("out of|file"))
  write_csv(df_feat_comb,
            file=paste0(dir_out, "/feat_combinations/feat_combinations_", fn, ".csv"))
}


get_thresholds_plot <- function(df_kdes, thresh) {
  ct_plot <- c("Neutrophil CD16hi", "T cell CD4 Naive", "T cell CD4 Mem",
               "T cell CD8 Naive", "T cell CD8 Mem",
               "B cell", "Monocyte Classical", "NK cell")
  ch_plot <- c("CD45RA", "CD27", "CD197", "CD4",
               "CD185", "CD183", "CD196",
               "CD38", "HLA-DR", "IgD", "CD56", "CD57", "CD161",
               "CD25", "CD127")

  tmp <- tibble(channel=names(thresh), intensity=unname(thresh)) %>%
    filter(channel %in% ch_plot)
  p <- plot_kdes(df_kdes %>% filter(channel %in% ch_plot), ct_plot) +
    geom_vline(data=tmp, mapping=aes(xintercept=intensity), linetype="dashed")
  return(p)
}


get_feat_major <- function(cell_type, unq, fn) {
  n_tot <- length(which(!grepl("Debris|Doublet|_", cell_type)))

  perc <- sapply(unq, function(ct) length(which(cell_type==ct))/n_tot)
  df_feat_maj <- tibble(ct=unq, perc=perc) %>%
    pivot_wider(names_from="ct", values_from="perc") %>%
    mutate(`T cell CD4` = rowSums(across(matches("T cell CD4"))),
           `T cell CD8` = rowSums(across(matches("T cell CD8"))),
           .keep="unused") %>%
    mutate(`T cell` = rowSums(across(matches("T cell"))),
           `Monocyte` = rowSums(across(matches("Monocyte"))),
           `Neutrophil` = rowSums(across(matches("Neutrophil"))),
           file = fn) %>%
    mutate(across(contains("T cell "), ~.x/`T cell`)) %>%
    mutate(across(contains("Monocyte "), ~.x/`Monocyte`)) %>%
    mutate(across(contains("Neutrophil "), ~.x/`Neutrophil`)) %>%
    relocate(file)
  return(df_feat_maj)
}


refine_subsets <- function(data_scaled, cell_type, ct, defs, only_pred=FALSE) {
  cells <- which(cell_type==ct)
  n_cells <- length(cells)

  if(n_cells > 0) { # > 20
    cols_cells <- intersect(names(defs), colnames(data_scaled))
    x_cells <- data_scaled[cells,cols_cells,drop=FALSE]
    pred_cells <- predict_cell_type(x_cells, defs)
    cell_type[cells] <- pred_cells
  }

  if (only_pred)
    return(pred_cells)

  return(cell_type)
}


update_thresh_adhoc <- function(df, cell_type, unq, thresh, ct) {
  if (ct == "T cell CD8") {
    thresh["CD45RA"] <- get_thresh_opt(df, cell_type, unq, thresh, "CD45RA",
                                       "T cell CD8 Mem", "T cell CD8 Naive")
  }

  if (ct == "T cell CD4") {
    thresh["CD45RA"] <- get_thresh_opt(df, cell_type, unq, thresh, "CD45RA",
                                       "T cell CD4 Mem", "T cell CD4 Naive")
  }

  if (ct == "B cell") {
    thresh["CD38"] <- thresh["CD38"]/2
  }

  if (ct == "NK cell") {
    x56 <- df %>%
      filter(cell_type=="NK cell" & CD57 > 5*sinh(thresh["CD57"])) %>%
      pull(CD56)
    thresh["CD56"] <- asinh(quantile(x56,0.995)/5)

    x16 <- df %>%
      filter(cell_type=="NK cell") %>%
      pull(CD16)
    thresh["CD16"] <- 0.25*asinh(quantile(x16,0.999)/5)
  }

  return(thresh)
}


apply_gate_hierarchy <- function(df, cell_type, thresh,
                                 dir_out, fn, ct="T cell CD4") {
  unq <- unique(cell_type)

  thresh <- update_thresh_adhoc(df, cell_type, unq, thresh, ct)
  thresh <- 5 * sinh(thresh)

  ##### parse hierarchy
  gate_hierarchy <- gate_hierarchy_full %>%
    filter(grepl(ct, Parent))
  gr <- graph_from_edgelist(gate_hierarchy %>%
                              select(Parent, GateName) %>%
                              as.matrix())

  bf <- bfs(gr, root=ct, unreachable = FALSE)
  ord_gate <- V(gr)$name[bf$order[!is.na(bf$order)]]
  group_by_gate <- setNames(c(ct, gate_hierarchy$GateGroup),
                            c(ct,gate_hierarchy$GateName))
  ord_group <- unique(unname(group_by_gate[ord_gate[-1]]))

  ##### apply gates
  df_ct <- df %>% mutate(label=cell_type) %>%
    filter(grepl(ct, label) & !grepl("Debris|_", label))

  ct_list <- lapply(ord_group, function(gn) rep("Not in parent", nrow(df_ct)))
  names(ct_list) <- ord_group
  df_gated <- as_tibble(ct_list)

  for (group_name in ord_group) {
    group_context <- get_group_context(group_name, gate_hierarchy)
    parent_group <- unname(group_by_gate[group_context$parent_name])

    if (parent_group==ct) {
      parent_idx <- seq(nrow(df_ct))
      df_parent <- df_ct
    } else {
      parent_idx <- which(df_gated[[parent_group]]==group_context$parent_name)
      df_parent <- df_ct[parent_idx,]
    }

    df_gated[[group_name]][parent_idx] <- "Ungated"

    p <- plot_bin2d(df_parent,
                    channels = c(group_context$x_col, group_context$y_col),
                    scales = c("asinh", "asinh")) +
      ggtitle(fn, subtitle=paste(group_name, "out of", parent_group))

    for (gate_name in group_context$gate_names) {
      context <- get_context(gate_name, gate_hierarchy)
      gate <- get_gate_conditional(df_parent, context, thresh, gate_name)

      filt <- get_filter_conditional(df_parent, gate, gate_name)
      df_gated[[group_name]][parent_idx][which(filt)] <- gate_name

      p <- update_gate_plot(p, filt, gate, gate_name)
    }

    ggsave(p, filename=paste0(dir_out, "gating/", group_name, "/", fn, ".png"),
           width=6, height=6)
  }

  df_gated <- simplify_gates(df_gated, gate_hierarchy, ord_group)
  df_feat <- tabulate_gate_groups(df_gated, ct)
  return(df_feat)
}


get_gate_conditional <- function(df_parent, context, thresh, gate_name) {
  if (grepl("Act", gate_name)) {
    gate <- get_act_gate(df_parent, c(context$x_col, context$y_col), thresh)
  } else if (grepl("Treg", gate_name)) {
    gate <- get_treg_gate(df_parent, c(context$x_col, context$y_col), thresh)
  } else {
    gate <- get_rectangle_gate(df_parent, thresh, context)
  }
  return(gate)
}


get_filter_conditional <- function(df_parent, gate, gate_name) {
  pattern <- "CD4 Naive|CD8 Naive"
  if (grepl(pattern, gate_name)) {
    filt <- df_parent$label==gate_name
  } else {
    filt <- apply_gate(df_parent, gate) & !grepl(pattern, df_parent$label)
  }
  return(filt)
}


update_gate_plot <- function(p, filt, gate, gate_name) {
  rat <- length(which(filt))/length(filt)
  label <- paste0(gate_name, "\n", round(100*rat,2), "%")
  p <- p + geom_path(data=gate, color="red") +
    annotate("label", x=mean(gate[[1]]), y=mean(gate[[2]]),
             label=label, alpha=0.5)
  return(p)
}


simplify_gates <- function(df_gated, gate_hierarchy, ord_group) {
  nonrep_gates <- gate_hierarchy %>% filter(!Report) %>% pull(GateName)
  if (length(nonrep_gates)==0)
    return(df_gated)

  target_grps <- gate_hierarchy %>% filter(Parent %in% nonrep_gates) %>%
    pull(GateGroup) %>% unique()
  parent_grps <- gate_hierarchy %>% filter(GateName %in% nonrep_gates) %>%
    pull(GateGroup)

  ord <- match(ord_group, target_grps)
  ord <- rev(ord[which(!is.na(ord))])

  target_grps <- target_grps[ord]
  parent_grps <- parent_grps[ord]

  for (i in seq_along(target_grps)) {
    grp <- target_grps[i]
    prt <- parent_grps[i]
    idx <- which(df_gated[[grp]] != "Not in parent")
    df_gated[[prt]][idx] <- df_gated[[grp]][idx]
    df_gated[[grp]] <- NULL
  }

  return(df_gated)
}


search_opt <- function(x, class1_idx, class2_idx, w, size=20) {
  candidates <- seq(from=min(x), to=max(x), length.out=size)[seq(2,size-1)]

  scores <- sapply(candidates, function(cand) {
    tn <- length(which(x[class1_idx] < cand))
    fp <- length(which(x[class1_idx] >= cand))

    fn <- length(which(x[class2_idx] < cand))
    tp <- length(which(x[class2_idx] >= cand))

    tnr <- tn / (tn + fp)
    tpr <- tp / (tp + fn)
    return((w[1]*tnr+w[2]*tpr)/(2*sum(w)))
  })

  return(candidates[which.max(scores)])
}


tabulate_gate_groups <- function(df_gated, ct) {
  df_single <- lapply(df_gated, tabulate_group_single) %>% do.call(what=rbind)

  grps <- names(df_gated)
  if (length(grps)<2)
    return(df_single)

  grp_pairs <- combn(grps, 2, simplify=FALSE)

  df_pairs <- lapply(grp_pairs, tabulate_group_pair, df_gated=df_gated, ct=ct) %>%
    do.call(what=rbind)

  return(rbind(df_single, df_pairs))
}


tabulate_group_pair <- function(df_gated, pair, ct) {
  tab <- table(df_gated[[pair[1]]], df_gated[[pair[2]]])

  perc_mat <- apply(tab, 1, function(row) row/sum(row))
  feat_row <- cross_perc(perc_mat, ct)

  perc_mat <- apply(tab, 2, function(col) col/sum(col))
  feat_col <- cross_perc(perc_mat, ct)

  return(rbind(feat_row, feat_col))
}


cross_perc <- function(perc_mat, ct) {
  lapply(setdiff(rownames(perc_mat), "Ungated"), function(gate_name) {
    res <- perc_mat[gate_name,setdiff(colnames(perc_mat), "Ungated"),drop=FALSE]
    tmp <- tibble(feature = paste0(gate_name, " out of", str_remove(colnames(res), ct)),
                  frac = as.numeric(res))
    return(tmp)
  }) %>% do.call(what=rbind)
}


filter_bead_gaussian_live <- function(df_raw, event_type) {
  channels_remove <- "Time|Bead|Live|Center|Offset|Residual|Width"
  channels_keep <- names(df_raw)[!grepl(channels_remove, names(df_raw))]

  df <- df_raw %>%
    filter(event_type == "") %>%
    select(all_of(channels_keep))
  return(df)
}


tabulate_group_single <- function(col) {
  tab <- table(col)
  in_parent <- setdiff(names(tab), "Not in parent")
  to_report <- setdiff(in_parent, "Ungated")
  denom <- sum(tab[in_parent])

  return(tibble(feature=to_report, frac=as.numeric(tab[to_report])/denom))
}


get_treg_gate <- function(df, params, thresh) {

  t1 <- thresh[params[1]]
  t2 <- thresh[params[2]]
  M <- max(df[[params[1]]])
  delta_x <- min(asinh(t2/5), asinh(M/5)-asinh(t1/5))

  gate <- tibble(var1 = c(t1, 5*sinh(asinh(t1/5)+delta_x), M, M, t1, t1),
                 var2 = c(t2, 5*sinh(asinh(t2/5)+delta_x), 5*sinh(asinh(t2/5)+delta_x), 0, 0, t2))
  names(gate) <- params
  return(gate)
}


get_act_gate <- function(df, params, thresh) {
  t1 <- asinh(thresh[params[1]]/5)
  t2 <- asinh(thresh[params[2]]/5)
  M1 <- max(df[[params[1]]])
  M2 <- max(df[[params[2]]])

  gate <- tibble(var1 = c(M1, 5*sinh(0.9*t1), 5*sinh(0.9*t1), 5*sinh(1.25*t1), M1, M1),
                 var2 = c(M2, M2, 5*sinh(1.25*t2), 5*sinh(0.75*t2), 5*sinh(0.75*t2), M2))
  names(gate) <- params
  return(gate)
}


estimate_distributions <- function(cell_type, data_transf, fn) {

  cell_types <- c("Neutrophil CD16hi", "Eosinophil", "Basophil", "B cell",
                  "Monocyte Classical", "T cell gd", "NK cell",
                  "T cell CD4 Naive", "T cell CD4 Mem",
                  "T cell CD8 Naive", "T cell CD8 Mem")
  channels <- setdiff(colnames(data_transf), c("DNA1", "DNA2", "Event_length"))

  tab <- table(cell_type)
  to_use <- intersect(cell_types, names(which(tab > 50)))
  to_use <- union(to_use, "all")

  df_kdes <- lapply(to_use, function(ct) {

    if (ct == "all") {
      cells <- which(!grepl("Debris|Doublet|_", cell_type))
    } else {
      cells <- which(cell_type == ct)
    }

    lapply(channels, function(ch) {
      x <- data_transf[cells,ch]
      kde <- bkde(x, range.x=c(-1,8), gridsize = 101L)

      return(tibble(file=fn, cell_type=ct, channel=ch,
                    expression = kde$x,
                    density = round(kde$y / max(kde$y), 4)))
    }) %>% do.call(what=rbind)
  }) %>% do.call(what=rbind)

  return(df_kdes)
}


get_hull <- function(data, context, gate_col) {
  df <- data %>%
    filter(.data[[gate_col]] %in% context$this_descendants) %>%
    select(all_of(c(context$x_col, context$y_col)))
  idx <- chull(df)
  hull <- df[c(idx, idx[1]),]
  return(hull)
}


get_group_context <- function(group_name, gate_hierarchy) {
  tmp <- gate_hierarchy %>% filter(GateGroup==group_name)
  return(list(parent_name = tmp$Parent[1],
              x_col = tmp$Channel1[1],
              y_col = tmp$Channel2[1],
              gate_names = tmp$GateName))
}


get_context <- function(gate_name, gate_hierarchy) {
  # build tree as an igraph object
  gr <- graph_from_edgelist(gate_hierarchy %>%
                              select(GateName, Parent) %>%
                              filter(!is.na(Parent)) %>%
                              as.matrix())

  # get all gate names downstream of current gate:
  # those for which there exists a path to the current node
  d <- distances(gr, to=gate_name, mode="out")
  this_descendants <- names(which(!is.infinite(d[,1])))

  # pull parent name and get everything downstream of it, too
  parent_name <- gate_hierarchy %>% filter(GateName==gate_name) %>% pull(Parent)
  d <- distances(gr, to=parent_name, mode="out")
  parent_descendants <- names(which(!is.infinite(d[,1])))

  # channels to gate on
  x_col <- gate_hierarchy %>% filter(GateName==gate_name) %>% pull(Channel1)
  y_col <- gate_hierarchy %>% filter(GateName==gate_name) %>% pull(Channel2)

  # positive or negative for channels
  x_loc <- gate_hierarchy %>% filter(GateName==gate_name) %>% pull(Location1)
  y_loc <- gate_hierarchy %>% filter(GateName==gate_name) %>% pull(Location2)

  context <- list(this_name = gate_name,
                  this_descendants = this_descendants,
                  parent_name = parent_name,
                  parent_descendants = parent_descendants,
                  x_col = x_col,
                  y_col = y_col,
                  x_loc = x_loc,
                  y_loc = y_loc)
  return(context)
}


rectangle_gate <- function(df, thresholds, x_col, y_col, x_loc, y_loc) {
  x_lim <- if_else(x_loc=="hi", max(df[[x_col]]), min(df[[x_col]]))
  y_lim <- if_else(y_loc=="hi", max(df[[y_col]]), min(df[[y_col]]))

  x_th <- thresholds[x_col]
  y_th <- thresholds[y_col]

  gate <- tibble(x = c(x_th, x_th, x_lim, x_lim, x_th),
                 y = c(y_th, y_lim, y_lim, y_th, y_th))
  names(gate) <- c(x_col, y_col)
  return(gate)
}


half_plane_gate <- function(df, thresholds, x_col, y_col, x_loc) {
  x_lim <- if_else(x_loc=="hi", max(df[[x_col]]), min(df[[x_col]]))
  x_th <- thresholds[x_col]

  y_min <- min(df[[y_col]])
  y_max <- max(df[[y_col]])

  gate <- tibble(x = c(x_th, x_th, x_lim, x_lim, x_th),
                 y = c(y_min, y_max, y_max, y_min, y_min))
  names(gate) <- c(x_col, y_col)
  return(gate)
}


get_rectangle_gate <- function(df, thresholds, context) {
  if (is.na(context$y_loc)) {
    gate <- half_plane_gate(df, thresholds,
                            context$x_col, context$y_col,
                            context$x_loc)
  } else {
    gate <-   rectangle_gate(df, thresholds,
                             context$x_col, context$y_col,
                             context$x_loc, context$y_loc)
  }

  return(gate)
}


apply_gate <- function(df, gate) {
  x_col <- names(gate)[1]
  y_col <- names(gate)[2]

  pip <- point.in.polygon(df %>% pull(x_col),
                          df %>% pull(y_col),
                          gate %>% pull(x_col),
                          gate %>% pull(y_col))
  return(pip > 0)
}


clean_fcs_file <- function(ff, event_type) {
  pattern <- "Debris|Doublet|_|Bead|Offset|Residual|Width|Center|Dead"
  cell_idx <- which(!grepl(pattern, event_type))
  ff_clean <- ff[cell_idx,]
  return(ff_clean)
}


get_isotope_mass <- function(pdata) {
  idx <- which(grepl("_", pdata$desc) & !grepl("Bead|Live|DNA", pdata$desc))
  mark <- pdata$desc[idx] %>% str_split("_") %>% sapply("[",2)
  iso <- pdata$desc[idx] %>% str_split("_") %>% sapply("[",1)
  mass <- as.numeric(gsub("([0-9]+).*$", "\\1", iso))
  return(setNames(mass,mark))
}


check_correlations <- function(df, pdata, data_scaled, cell_type, fn, dir_out,
                               cutoff=0.6) {

  cell_types <- c("Neutrophil CD16hi", "B cell",
                  "Monocyte Classical", "NK cell",
                  "T cell CD4 Naive", "T cell CD4 Mem",
                  "T cell CD8 Naive", "T cell CD8 Mem")
  channels <- setdiff(names(df), c("DNA1", "DNA2", "Event_length"))
  channel_pairs <- combn(channels,2)
  channel_pairs <- channel_pairs[,which(!(grepl("CD45",channel_pairs[1,]) &
                                          grepl("CD45", channel_pairs[2,])))]
  n <- ncol(channel_pairs)

  df_cor <- lapply(cell_types, function(ct) {
    cells <- which(cell_type==ct)

    if (length(cells) < 100)
      return(NULL)

    data_cells <- data_scaled[cells,]

    cors <- lapply(seq(n), function(pair_index) {
      ch1 <- channel_pairs[1,pair_index]
      ch2 <- channel_pairs[2,pair_index]
      r <- cor(data_cells[,ch1], data_cells[,ch2])
      return(r)
    }) %>% do.call(what=c)
    return(tibble(cell_type=ct,
                  ch1=channel_pairs[1,],
                  ch2=channel_pairs[2,],
                  r=cors))
  }) %>% do.call(what=rbind) %>%
    arrange(-r, ch1, ch2) %>%
    mutate(file=fn) %>%
    relocate(file)

  df_pairs <- df_cor %>%
    filter(r>cutoff) %>%
    select(ch1, ch2) %>%
    unique()

  df_plot <- df %>% mutate(`Cell Type`=cell_type) %>%
    filter(`Cell Type` %in% cell_types)

  iso_mass <- get_isotope_mass(pdata)

  for (pair_index in seq(nrow(df_pairs))) {
    plot_correlation_pair(df_pairs, df_cor, df_plot, pair_index,
                          iso_mass, dir_out, fn, cutoff)
  }

  return(df_cor)
}


plot_correlation_pair <- function(df_pairs, df_cor, df_plot, pair_index,
                                  iso_mass, dir_out, fn, cutoff) {
  chan1 <- df_pairs$ch1[pair_index]
  chan2 <- df_pairs$ch2[pair_index]

  df_annot <- df_cor %>% mutate(`Cell Type`=cell_type) %>%
    filter(ch1==chan1 & ch2==chan2) %>%
    mutate(label = paste0("r=",round(r,2)))
  m1 <- asinh(min(df_plot[[chan1]])/5)
  m2 <- asinh(min(df_plot[[chan2]])/5)
  M1 <- asinh(max(df_plot[[chan1]])/5)
  M2 <- asinh(max(df_plot[[chan2]])/5)

  p <- plot_bin2d(df_plot, channels = c(chan1, chan2),
                  scales = c("asinh", "asinh")) +
    geom_label(data=df_annot, x=0.5*(m1+M1), y=0.9*M2+0.1*m2,
               aes(label=label, color=r>cutoff), alpha=0.5) +
    scale_color_manual(values=setNames(c("red", "black"), c("TRUE", "FALSE"))) +
    guides(color="none") +
    facet_wrap(~`Cell Type`, ncol=4) +
    ggtitle(label=fn, subtitle=paste0(chan1, " vs ", chan2, "; mass diff: |",
                                      iso_mass[chan1], "-", iso_mass[chan2],
                                      "|=", abs(iso_mass[chan1]-iso_mass[chan2])))
  filename <- paste0(dir_out, "correlations_figs/", fn, "_",
                     chan1, "_", chan2, ".png")
  ggsave(p, filename=filename, width=12, height=6)
}

