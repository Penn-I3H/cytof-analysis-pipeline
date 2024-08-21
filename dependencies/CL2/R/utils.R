
read_data <- function(ff) {
  pdata <- ff@parameters@data

  # remove unused channels
  non_metal <- grep("Time|Center|Offset|Width|Residual|length", pdata$name)
  markers <- setdiff(grep("_", pdata$desc), non_metal)

  data <- cbind(asinh(ff@exprs[,markers]/5),
                ff@exprs[,non_metal])
  data[,c("Center", "Offset", "Width", "Residual")] <- asinh(data[,c("Center", "Offset", "Width", "Residual")]/5)

  # use protein rather than isotope for channel names
  nice_names <- pdata$desc[markers] %>%
    str_split("_") %>%
    sapply(function(x) x[2])
  bead_ch <- which(nice_names=="Bead")
  nice_names[bead_ch] <- pdata$desc[markers[bead_ch]]
  colnames(data)[seq_along(markers)] <- nice_names

  return(data)
}


plot_cleanup_gate <- function(df, cutoffs, marker) {
  gg_flow(df, params=c("Time", marker)) +
    geom_hline(yintercept=cutoffs[[marker]][1], color="red") +
    geom_hline(yintercept=cutoffs[[marker]][2], color="red") +
    guides(fill="none")
}

apply_cleanup_gate <- function(df, event_type, cutoffs, marker, label) {
  filt <- which(event_type=="")
  dat <- df[[marker]][filt]

  idx_out <- which(dat < cutoffs[[marker]][1] | dat > cutoffs[[marker]][2])
  event_type[filt[idx_out]] <- label
  return(event_type)
}


pregate_data <- function(df, fn, dir_out, plot=FALSE) {

  cutoffs <- list("140Ce_Bead" = c(0,asinh(500/5)),
                  "Residual" = c(asinh(20/5), asinh(400/5)),
                  "Center" = c(asinh(300/5), asinh(4000/5)),
                  "Offset" = c(asinh(30/5), asinh(400/5)),
                  "Width" = c(asinh(30/5), asinh(300/5)),
                  "Live" = c(0, asinh(80/5)))
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
    p1 <- plot_cleanup_gate(df, cutoffs, markers[1])
    p2 <- plot_cleanup_gate(df1, cutoffs, markers[2])
    p3 <- plot_cleanup_gate(df2, cutoffs, markers[3])
    p4 <- plot_cleanup_gate(df3, cutoffs, markers[4])
    p5 <- plot_cleanup_gate(df4, cutoffs, markers[5])
    p6 <- plot_cleanup_gate(df5, cutoffs, markers[6])
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
  write_csv(df_stats, file = paste0(dir_out, "cleanup_stats/cleanup_stats_", fn, ".csv"),
            progress=FALSE)

  return(event_type)
}



just_pregate <- function(file, dir_in, dir_out) {
  set.seed(0)
  path <- paste0(dir_in, file)
  fn <- file %>%
    str_remove("_Normalized.fcs") %>%
    str_remove("_normalized.fcs") %>%
    str_remove("_Processed.fcs")
  # message(fn)

  df <- path %>%
    read_data() %>%
    as_tibble() %>%
    pregate_data(fn, dir_out, plot=TRUE)

  return(NULL)
}




scale_data <- function(mat) {
  mat_sc <- apply(mat, 2, function(col) {
    col / quantile(col,0.999)
  })
  return(mat_sc)
}


gg_flow <- function(df, params, nbin=100, xlim=NULL, ylim=NULL, method="bin") {
  mat <- as.matrix(df %>% select(all_of(params)))

  if (is.null(xlim)) {
    m1 <- min(mat[,1])
    M1 <- max(mat[,1])
  } else {
    m1 <- xlim[1]
    M1 <- xlim[2]
  }
  if (is.null(ylim)) {
    m2 <- min(mat[,2])
    M2 <- max(mat[,2])
  } else {
    m2 <- ylim[1]
    M2 <- ylim[2]
  }

  bins <- bin2(mat, nbin = c(nbin,nbin), ab=c(m1,m2,M1,M2))
  df_bin <- melt(bins$nc)
  df_bin$Var1 <- seq(m1,M1,length.out=nbin)[df_bin$Var1]
  df_bin$Var2 <- seq(m2,M2,length.out=nbin)[df_bin$Var2]
  names(df_bin) <- c(params, "count")

  ggplot(df_bin %>% filter(count > 0),
         aes(x=.data[[params[[1]]]], y=.data[[params[[2]]]])) +
    geom_tile(aes(fill=count)) +
    geom_hline(yintercept = 0, linetype = "dotdash") +
    geom_vline(xintercept = 0, linetype = "dotdash") +
    xlim(c(m1-0.1,M1+0.1)) +
    ylim(c(m2-0.1,M2+0.1)) +
    scale_fill_gradientn(colours = hsv(h = seq(0.6667,0, length.out = 11))) +
    theme_bw(base_size=13)
}


gg_flow_thresh <- function(df, params, thresh, thresh_neut, thresh_mono,
                           nbin=100, method="bin",
                           xlim=c(-0.1,7), ylim=c(-0.1,7)) {
  gg_flow(df, params, nbin=nbin, method=method, xlim=xlim, ylim=ylim) +
    geom_hline(yintercept=thresh[params[2]], color="red") +
    geom_vline(xintercept=thresh[params[1]], color="red")+
    geom_hline(yintercept=thresh_neut[params[2]], color="blue", linetype="dashed") +
    geom_vline(xintercept=thresh_neut[params[1]], color="blue", linetype="dashed")+
    geom_hline(yintercept=thresh_mono[params[2]], color="green", linetype="dashed") +
    geom_vline(xintercept=thresh_mono[params[1]], color="green", linetype="dashed")
}

gg_flow_gate <- function(df, params, gate, nbin=100, method="bin",
                         xlim=c(-0.1,7), ylim=c(-0.1,7)) {
  gg_flow(df, params, nbin=nbin, method=method, xlim=xlim, ylim=ylim) +
    geom_path(data=gate, color="red")
}


get_umap <- function(df, x, sel_umap, dir1="DNA1", dir2="CD3") {
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
    scale_color_gradient(low="black", high="red") +
    theme_bw(base_size=base_size) +
    theme(axis.title = element_blank())
}


plot_umap_clust <- function(df, channel, fn, pal=NULL, base_size=8) {
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


plot_umap_major <- function(df, fn) {
  pal_paired <- brewer.pal(12, "Paired")
  names(pal_paired) <- c("neutrophil", "eosinophil",
                         "bcell", "plasmablast",
                         "basophil", "pdc",
                         "tcell", "nkcell",
                         "bcell_lymphoma", "myeloid",
                         "debris", "doublet")

  p1 <- plot_umap_clust(df, "ct", fn, pal_paired)
  p2 <- plot_umap_channel(df, "CD3")
  p3 <- plot_umap_channel(df, "CD19")
  p4 <- plot_umap_channel(df, "CD20")
  p5 <- plot_umap_channel(df, "CD66b")
  p6 <- plot_umap_channel(df, "CD16")
  p7 <- plot_umap_channel(df, "CD123")
  p8 <- plot_umap_channel(df, "CD294")
  p9 <- plot_umap_channel(df, "CD56")
  p10 <- plot_umap_channel(df, "CD11c")
  p11 <- plot_umap_channel(df, "DNA1")
  p12 <- plot_umap_channel(df, "CD45")

  p <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol=4)
  return(p)
}


plot_dna_cd45 <- function(df, fn) {
  dark2 <- brewer.pal(8, "Dark2")
  pal_dark2 <- c("debris"=dark2[6], "single cell"=dark2[1], "doublet"=dark2[2])

  ggplot(df, aes(x=DNA1, y=CD45, color=event)) +
    geom_point(size=0.5, shape=1, alpha=0.5) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape=19))) +
    scale_color_manual(values=pal_dark2) +
    ggtitle(fn) +
    theme_bw()
}


plot_umap_tcell <- function(df, fn) {

  pal_tcell <- brewer.pal(5, "Set1")
  names(pal_tcell) <- c("tcell_cd4", "tcell_cd8", "tcell_gd", "tcell_dn", "tcell_dp")

  p1 <- plot_umap_clust(df, "clust", fn, pal_tcell)
  p2 <- plot_umap_channel(df, "CD4")
  p3 <- plot_umap_channel(df, "CD8a")
  p4 <- plot_umap_channel(df, "TCRgd")
  p <- wrap_plots(p1, p2, p3, p4, ncol=2)
  return(p)
}


plot_umap_mono <- function(df, fn) {

  pal <- brewer.pal(5, "Set1")
  names(pal) <- c("monocyte_classical", "monocyte_nonclassical", "mdc", "pdc", "basophil")

  p1 <- plot_umap_clust(df, "clust", fn, pal)
  p2 <- plot_umap_channel(df, "CD14")
  p3 <- plot_umap_channel(df, "CD38")
  p4 <- plot_umap_channel(df, "CD16")
  p5 <- plot_umap_channel(df, "HLA-DR")
  p6 <- plot_umap_channel(df, "CD11c")
  p7 <- plot_umap_channel(df, "CD66b")
  p8 <- plot_umap_channel(df, "CD123")
  p9 <- plot_umap_channel(df, "CD45RA")
  p <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)
  return(p)
}


run_fastpg <- function(data, resolution=1, n_threads=1) {
  all_knn <- hnsw_knn(data, k=15, distance= 'l2',
                      n_threads=n_threads, M=48)
  el <- get_edgelist(all_knn$idx)

  gr <- graph_from_edgelist(el$edgelist, directed=FALSE)
  E(gr)$weight <- el$jac

  leid <- cluster_leiden(gr, objective_function = "modularity",
                         resolution_parameter = resolution)
  clustering <- as.factor(membership(leid))
  return(clustering)
}


get_medians <- function(data, clustering) {
  medians <- sapply(levels(clustering), function(lev) {
    data[which(clustering==lev),,drop=FALSE] %>%
      apply(2, median)
  }) %>% t()

  return(medians)
}


label_clusters_score_opt <- function(centroids, defs, mdipa_main=TRUE, return_ypred=FALSE) {
  coeff <- defs %>% select(-Phenotype)

  coeff <- as.matrix(coeff)
  rownames(coeff) <- defs$Phenotype

  ypred <- exp(t(coeff[,-1,drop=FALSE] %*% t(centroids[,colnames(coeff[,-1,drop=FALSE]),drop=FALSE]) + coeff[,1])) %>%
    apply(1, function(row) row/sum(row)) %>% t()

  if (return_ypred)
    return(ypred)

  lab <- colnames(ypred)[unname(apply(ypred, 1, which.max))]
  return(lab)
}


get_label_confidence <- function(centroids, defs) {
  ypred <- label_clusters_score_opt(centroids, defs, return_ypred = TRUE)

  conf <- apply(ypred, 1, function(row) {
    row_s <- sort(row, decreasing=TRUE)
    return(row_s[1]-row_s[2])
  })
  return(conf)
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


plot_kdes_all <- function(kdes) {
  ggplot(kdes, aes(x=protein_expression, y=density_estimate)) +
    geom_line() +
    geom_vline(xintercept = 0, linetype="dotted") +
    facet_wrap(~channel) +
    theme_bw()
}

plot_kdes <- function(df_kdes, sel_vals, sel_col="cluster") {
  df_plot <- df_kdes %>%
    filter(.data[[sel_col]] %in% sel_vals)

  ggplot(df_plot, aes(x=intensity, y=density, color=.data[[sel_col]])) +
    geom_path() +
    facet_wrap(~channel, scales="free") +
    theme_bw()
}


get_kdes <- function(df, cols, clustering) {
  levs <- levels(clustering)

  df_kdes <- lapply(cols, function(ch) {
    m <- min(df[[ch]])
    M <- max(df[[ch]])
    lapply(levs, function(lev) {
      dat <- df[[ch]][which(clustering==lev)]

      if (length(dat) < 10)
        dat <- jitter(rep(dat,10))
      else {
        if (sd(dat)==0)
          dat <- jitter(dat)
      }


      kde <- bkde(dat, range.x = c(m,M), gridsize=101L)
      return(tibble(cluster=lev,
                    channel=ch,
                    intensity=kde$x,
                    density=kde$y / sum(kde$y)))
    }) %>% do.call(what=rbind)
  }) %>% do.call(what=rbind) %>%
    mutate(density = pmax(density,0))

  return(df_kdes)
}


compute_bhatt <- function(x,y) {
  return(sum(sqrt(x * y)))
}


conv_sum <- function(f1, f2, intensity, scal="asinh") {

  res <- sapply(intensity, function(z) {
    if (scal=="linear") {
      ## TO DO: understand Event_length better, justify why +15 works
      dif_transf <- z-intensity+15   # intensity = y
      # chain <- 1
    }
    else {
      tmp <- sinh(z)-sinh(intensity)
      dif_transf <- asinh(tmp) # x
      # chain <- cosh(dif_transf)*1/sqrt(1+z^2)
    }

    s <- sum(f1(intensity) * f2(dif_transf))
    # s <- sum(chain * f1(intensity) * f2(dif_transf))
    return(s)
  })
  # return(res)
  return(res/sum(res))
}


merge_clusters_c <- function(df, cols, clustering, min_bhatt=0.1) {
  df_kdes <- get_kdes(df, cols, clustering) %>%
    pivot_wider(names_from="cluster", values_from="density")
  mat <- as.matrix(df_kdes %>% select(-channel, -intensity))

  bhatt_all <- bhatt(mat, length(cols), 101L)
  bhatt_summ <- bhatt_all %>%
    melt() %>%
    filter(Var1 <= Var2) %>%
    mutate(value = replace_na(value, 0)) %>%
    filter(value >= min_bhatt)
  # filter(value >= 0.88^length(cols))

  gr <- graph_from_edgelist(bhatt_summ %>%
                                      select(Var1,Var2) %>%
                                      as.matrix(),
                                    directed=FALSE)
  comp <- components(gr)$membership
  # names(comp) <- levels(clustering)

  cl_merged <- comp[clustering] %>% unname() %>% as.factor()
  return(cl_merged)
}

detect_doublets <- function(df, cols, cell_type, thresh=5) {

  x0 <- as.matrix(df)[,setdiff(cols, "Event_length")]

  set.seed(0)
  sel_for_aug <- which(cell_type != "debris")
  n_doub <- floor(length(sel_for_aug)/2)

  sel1 <- sample(sel_for_aug, n_doub)
  sel2 <- sample(sel_for_aug, n_doub)
  x_doub <- asinh(sinh(x0[sel1,]) + sinh(x0[sel2,]))

  x_aug <- rbind(x_doub, x0[sel_for_aug,])

  all_knn <- hnsw_knn(x_aug, k=15, distance= 'l2',
                      n_threads=1, M=48)
  n_nb_doub <- apply(all_knn$idx, 1, function(row) length(which(row<n_doub)))

  sim_lab <- if_else(cell_type[sel1] < cell_type[sel2],
                     paste0("doublet_", cell_type[sel1], "_", cell_type[sel2]),
                     paste0("doublet_", cell_type[sel2], "_", cell_type[sel1]))
  sim_lab_fac <- as.factor(sim_lab)

  knn_orig <- all_knn$idx[-seq(n_doub),]

  doub_c <- doub_lab_c(knn_orig, sim_lab_fac, n_doub, length(levels(sim_lab_fac)))
  doub_lab <- c("",levels(sim_lab_fac))[doub_c+1]

  ct_final <- cell_type
  ct_final[sel_for_aug] <- if_else(n_nb_doub[-seq(n_doub)] >= thresh,
                                   doub_lab,
                                   cell_type[sel_for_aug])
  return(ct_final)
}


get_bhatt_target_fast <- function(df_kdes, tab, med_len_dna, cols, target) {
  levs_mo <- intersect(names(which( tab > tab[target] )),
                       names(which( med_len_dna[,1] < med_len_dna[target,1] &
                                      med_len_dna[,2] < med_len_dna[target,2]))) %>%
    sort()
  n <- length(levs_mo)

  if(n<1)
    return(NULL)

  bhatt_all <- lapply(cols, function(ch) {
    df_kdes_ch <- df_kdes %>% filter(channel==ch)
    intensity <- df_kdes_ch[["intensity"]]
    dist_target <- df_kdes_ch[[target]]
    scal <- if_else(ch=="length", "linear", "asinh")

    lapply(seq(1,n), function(i) {
      cl1 <- levs_mo[i]
      dist_cl1 <- df_kdes_ch[[cl1]]
      f1 <- approxfun(intensity, dist_cl1/sum(dist_cl1), yleft=0, yright=0)

      bhatt <- lapply(seq(i,n), function(j) {
        cl2 <- levs_mo[j]
        dist_cl2 <- df_kdes_ch[[cl2]]
        f2 <- approxfun(intensity, dist_cl2/sum(dist_cl2), yleft=0, yright=0)

        dist_ref <- conv_sum(f1,f2,intensity, scal=scal)
        compute_bhatt(dist_target, dist_ref)
      }) %>% do.call(what=c)

      return(tibble(target=target, cl1=cl1, cl2=levs_mo[seq(i,n)],
                    channel=ch, bhatt=bhatt))
    }) %>% do.call(what=rbind)
  }) %>% do.call(what=rbind)

  bhatt_summ <- bhatt_all %>%
    group_by(target, cl1, cl2) %>%
    summarise(bhatt = prod(bhatt), .groups = "drop")

  return(bhatt_summ)
}


one_thresh_gaussian <- function(x, sym_zero=TRUE) {
  q <- quantile(x, c(0.25, 0.5, 0.75, 0.95, 0.99, 0.999))
  med <- unname(q["50%"])

  if (sym_zero) {
    xpos <- x[which(x>med & x<q["99%"])]
    xfin <- c(xpos, 2*med - xpos)
  } else {
    xfin <- x[which(x<q["99%"])]
  }
  s <- sd(xfin)

  return(med+2.5*s)
}

get_thresh_gauss <- function(df) {
  thresh <- sapply(names(df), function(col) one_thresh_gaussian(df[[col]]))
  return(thresh)
}

update_pred <- function(cell_type, pred, conf) {
  pred <- as.character(pred)
  ct_update <- c("basophil", "eosinophil", "plasmablast", "pdc")

  for (ct in ct_update) {
    if ((! ct %in% cell_type) & (ct %in% pred)) {
      cells <- which(pred==ct & conf > 0.99)

      if (length(cells) > 0)
        cell_type[cells] <- ct
    }
  }
  return(cell_type)
}


update_clustering <- function(cl, cl_new, to_replace) {
  cl <- as.character(cl)
  cl[to_replace] <- as.character(cl_new)
  return(as.factor(cl))
}


gating_stuff <- function(df, cell_type, dir_out, fn) {
  df_kdes <- get_kdes(df, names(df), as.factor(cell_type))

  thresh_neut <- get_thresh_gauss(df %>% filter(cell_type=="neutrophil"))
  thresh_mono <- get_thresh_gauss(df %>% filter(cell_type=="myeloid"))
  thresh <- get_thresholds_density(df, cell_type)

  ct_plot <- c("neutrophil", "tcell_cd4", "tcell_cd8", "bcell", "myeloid", "nkcell")
  ch_plot <- c("CD45RA", "CD27", "CD197", "CD14", "CD4",
               "CD25", "CD185", "CD183", "CD196", "CD56",
               "CD38", "HLA-DR", "IgD", "CD57", "CD161")
  ch_neut <- setdiff(ch_care, c("CD196"))
  ch_mono <- setdiff(ch_care, c("CD14", "CD196", "CD38", "CD45RA", "HLA-DR"))
  thresh[ch_neut] <- pmax(thresh[ch_neut], thresh_neut[ch_neut])
  thresh[ch_mono] <- pmax(thresh[ch_mono], thresh_mono[ch_mono])

  p <- plot_kdes(df_kdes %>% filter(channel %in% ch_plot), ct_plot) +
    geom_vline(data=tibble(channel=names(thresh), intensity=unname(thresh)) %>% filter(channel %in% ch_plot),
               mapping=aes(xintercept=intensity), linetype="dashed")
  ggsave(p, filename=paste0(dir_out, "gating/kdes_thresh/", fn, ".png"), width=12, height=10)

  n_tot <- length(which(!grepl("debris|doublet", cell_type)))
  n_tcell <- length(which(grepl("tcell", cell_type) & !grepl("doublet", cell_type)))
  df_feat_maj <- tibble(file = fn,
                        `T cell` = n_tcell / n_tot,
                        `T cell CD4` = length(which(cell_type=="tcell_cd4")) / n_tcell,
                        `T cell CD8` = length(which(cell_type=="tcell_cd8")) / n_tcell,
                        `T cell gd` = length(which(cell_type=="tcell_gd")) / n_tcell,
                        `B cell` = length(which(cell_type=="bcell")) / n_tot,
                        `NK cell` = length(which(cell_type=="nkcell")) / n_tot,
                        `Myeloid` = length(which(cell_type=="myeloid")) / n_tot,
                        `pDC` = length(which(cell_type=="pdc")) / n_tot,
                        `Plasmablast` = length(which(cell_type=="plasmablast")) / n_tot,
                        `Neutrophil` = length(which(cell_type=="neutrophil")) / n_tot,
                        `Eosinophil` = length(which(cell_type=="eosinophil")) / n_tot,
                        `Basophil` = length(which(cell_type=="basophil")) / n_tot)
  write_csv(df_feat_maj, file=paste0(dir_out, "/feat_major/feat_major_", fn, ".csv"))

  df_feat_adaptive <- gate_and_plot(df, cell_type, thresh, thresh_neut, thresh_mono,
                                    dir_out, fn, save_plots=TRUE)
  write_csv(df_feat_adaptive, file=paste0(dir_out, "/feat_adaptive/feat_adaptive_", fn, ".csv"))

}


apply_gate_hierarchy <- function(df, cell_type, ct_in="tcell_cd4", ct="T cell CD4") {

  gate_hierarchy <- read_csv("../gate_hierarchy.csv") %>%
    filter(grepl(ct, Parent))
  gr <- graph_from_edgelist(gate_hierarchy %>%
                              select(Parent, GateName) %>%
                              as.matrix())
  ggraph(gr, layout="tree") +
    geom_edge_diagonal() +
    geom_node_label(aes(label=str_wrap(name, width=10))) +
    theme_bw()

  bf <- bfs(gr, root=ct, unreachable = FALSE)
  ord <- V(gr)$name[bf$order[!is.na(bf$order)]]

  # thresholds <- get_thresholds(cell_type)

  ##### learn gates
  df_ct <- df %>% filter(cell_type == ct_in)

  ct_list <- list(ct = seq(nrow(df_ct)))
  names(ct_list) <- ct
  gate_list <- list()

  for (gate_name in ord[-1]) {
    print(gate_name)
    context <- get_context(gate_name, gate_hierarchy)

    if(context$non_naive) {
      parent_idx <- setdiff(ct_list[[context$parent_name]],
                            ct_list[[paste(ct, "Naive")]])
    } else {
      parent_idx <- ct_list[[context$parent_name]]
    }
    df_parent <- df_ct[parent_idx,]

    gate <- get_fd_gate(df_parent, context)
    # flowDensity::deGate(ff, channel="CD57")

    # gate_list[[gate_name]] <- gate
    filt <- apply_gate2(df_parent, gate)
    ct_list[[gate_name]] <- parent_idx[which(filt)]

    p <- plot_gate(df_parent, context, gate)
    ggsave(p, filename = paste0("../figures/gating_train/", ct, "/",
                                gate_name, ".png"),
           width=7, height=6)
  }


}


gate_and_plot <- function(df, cell_type, thresh, thresh_neut, thresh_mono,
                          dir_out, fn, save_plots=TRUE) {
  n_tot <- length(which(!grepl("debris|doublet", cell_type)))

  df_ab <- df %>% filter(grepl("cd4|cd8|dn|dp", cell_type))

  df_neut <- df %>% filter(cell_type=="neutrophil")
  df_myel <- df %>% filter(cell_type=="myeloid")
  # df_mono_cl <- df %>% filter(cell_type=="monocyte_classical")
  df_cd4 <- df %>%
    filter(cell_type=="tcell_cd4") %>%
    mutate(cell_id = seq(length(which(cell_type=="tcell_cd4")))) %>%
    relocate(cell_id)
  df_cd8 <- df %>% filter(cell_type=="tcell_cd8") %>%
    mutate(cell_id = seq(length(which(cell_type=="tcell_cd8")))) %>%
    relocate(cell_id)
  df_gd <- df %>% filter(cell_type=="tcell_gd")
  df_bcell <- df %>% filter(cell_type=="bcell")
  df_nkcell <- df %>% filter(cell_type=="nkcell")

  # df_mait <- df_ab %>% filter(CD4<thresh["CD4"] & CD161>thresh["CD161"])

  df_cd4_naive <- df_cd4 %>% filter(CD45RA>thresh["CD45RA"] & CD27>thresh["CD27"])
  df_cd4_emra <- df_cd4 %>% filter(CD45RA>thresh["CD45RA"] & CD27<thresh["CD27"])
  df_cd4_cd45ralo <- df_cd4 %>% filter(CD45RA<thresh["CD45RA"])

  df_cd4_em1 <- df_cd4_cd45ralo %>% filter(CD27>thresh["CD27"] & CD197<thresh["CD197"])
  df_cd4_em2 <- df_cd4_cd45ralo %>% filter(CD27<thresh["CD27"] & CD197>thresh["CD197"])
  df_cd4_em3 <- df_cd4_cd45ralo %>% filter(CD27<thresh["CD27"] & CD197<thresh["CD197"])
  df_cd4_cm  <- df_cd4_cd45ralo %>% filter(CD27>thresh["CD27"] & CD197>thresh["CD197"])

  df_cd8_naive <- df_cd8 %>% filter(CD45RA>thresh["CD45RA"] & CD27>thresh["CD27"])
  df_cd8_emra <- df_cd8 %>% filter(CD45RA>thresh["CD45RA"] & CD27<thresh["CD27"])
  df_cd8_cd45ralo <- df_cd8 %>% filter(CD45RA<thresh["CD45RA"])

  df_cd8_em1 <- df_cd8_cd45ralo %>% filter(CD27>thresh["CD27"] & CD197<thresh["CD197"])
  df_cd8_em2 <- df_cd8_cd45ralo %>% filter(CD27<thresh["CD27"] & CD197>thresh["CD197"])
  df_cd8_em3 <- df_cd8_cd45ralo %>% filter(CD27<thresh["CD27"] & CD197<thresh["CD197"])
  df_cd8_cm  <- df_cd8_cd45ralo %>% filter(CD27>thresh["CD27"] & CD197>thresh["CD197"])

  params <- c("CD38", "HLA-DR")
  gate_act <- get_act_gate(df_ab, params, thresh)
  df_cd4_act <- apply_gate(df_cd4, params, gate_act)
  df_cd8_act <- apply_gate(df_cd8, params, gate_act)

  df_cd4_nn <- df_cd4 %>% filter(CD45RA < thresh["CD45RA"] | CD27 < thresh["CD27"])

  params <- c("CD25", "CD127")
  gate_treg <- get_treg_gate(df_cd4_nn, params, thresh)
  df_cd4_treg <- apply_gate(df_cd4_nn, params, gate_treg)

  df_cd4_cd185p <- df_cd4_nn %>% filter(CD185 > thresh["CD185"])
  df_cd4_cd185n <- df_cd4_nn %>% filter(CD185 < thresh["CD185"])

  df_cd4_th1  <- df_cd4_cd185n %>% filter(CD183 > thresh["CD183"] & CD196 < thresh["CD196"])
  df_cd4_th17 <- df_cd4_cd185n %>% filter(CD183 < thresh["CD183"] & CD196 > thresh["CD196"])

  df_bcell_naive <- df_bcell %>% filter(IgD > thresh["IgD"] & CD27 < thresh["CD27"])
  df_bcell_mem <- df_bcell %>% filter( CD27 > thresh["CD27"])
  df_bcell_mem_sw <- df_bcell_mem %>% filter(IgD < thresh["IgD"])
  df_bcell_mem_nsw <- df_bcell_mem %>% filter(IgD > thresh["IgD"])

  df_nkcell_late <- df_nkcell %>% filter(CD57 > thresh["CD57"])

  df_mono_cl <- df_myel %>% filter(CD14 > thresh["CD14"] & CD38 > thresh["CD38"])
  df_mono_ncl <- df_myel %>% filter(CD14 < thresh["CD14"] & CD38 < thresh["CD38"])
  df_mdc <- df_myel %>% filter(CD14 < thresh["CD14"] & CD38 > thresh["CD38"])

  df_cd4_all <- df_cd4 %>%
    mutate(naive = cell_id %in% df_cd4_naive$cell_id,
           emra = cell_id %in% df_cd4_emra$cell_id,
           em1 = cell_id %in% df_cd4_em1$cell_id,
           em2 = cell_id %in% df_cd4_em2$cell_id,
           em3 = cell_id %in% df_cd4_em3$cell_id,
           cm = cell_id %in% df_cd4_cm$cell_id,
           treg = cell_id %in% df_cd4_treg$cell_id,
           tfh = cell_id %in% df_cd4_cd185p$cell_id,
           th1 = cell_id %in% df_cd4_th1$cell_id,
           th17 = cell_id %in% df_cd4_th17$cell_id,
           act = cell_id %in% df_cd4_act$cell_id,
           cd57hi = CD57 > thresh["CD57"]) %>%
    select(where(is.logical))

  df_cd8_all <- df_cd8 %>%
    mutate(naive = cell_id %in% df_cd8_naive$cell_id,
           emra = cell_id %in% df_cd8_emra$cell_id,
           em1 = cell_id %in% df_cd8_em1$cell_id,
           em2 = cell_id %in% df_cd8_em2$cell_id,
           em3 = cell_id %in% df_cd8_em3$cell_id,
           cm = cell_id %in% df_cd8_cm$cell_id,
           act = cell_id %in% df_cd8_act$cell_id,
           cd57hi = CD57 > thresh["CD57"]) %>%
    select(where(is.logical))

  df_feat <- tibble(file = fn,
                    # `MAIT/NKT` = nrow(df_mait) / nrow(df_ab),
                    `CD4 Naive` = nrow(df_cd4_naive) / nrow(df_cd4),
                    `CD4 EMRA` = nrow(df_cd4_emra) / nrow(df_cd4),
                    `CD4 EM1` = nrow(df_cd4_em1) / nrow(df_cd4),
                    `CD4 EM2` = nrow(df_cd4_em2) / nrow(df_cd4),
                    `CD4 EM3` = nrow(df_cd4_em3) / nrow(df_cd4),
                    `CD4 CM` = nrow(df_cd4_cm) / nrow(df_cd4),
                    `CD4 act` = nrow(df_cd4_act) / nrow(df_cd4),
                    `CD8 Naive` = nrow(df_cd8_naive) / nrow(df_cd8),
                    `CD8 EMRA` = nrow(df_cd8_emra) / nrow(df_cd8),
                    `CD8 EM1` = nrow(df_cd8_em1) / nrow(df_cd8),
                    `CD8 EM2` = nrow(df_cd8_em2) / nrow(df_cd8),
                    `CD8 EM3` = nrow(df_cd8_em3) / nrow(df_cd8),
                    `CD8 CM` = nrow(df_cd8_cm) / nrow(df_cd8),
                    `CD8 act` = nrow(df_cd8_act) / nrow(df_cd8),
                    `Treg` = nrow(df_cd4_treg) / nrow(df_cd4_nn),
                    `Tfh` = nrow(df_cd4_cd185p) / nrow(df_cd4_nn),
                    `Th1` = nrow(df_cd4_th1) / nrow(df_cd4_nn),
                    `Th17` = nrow(df_cd4_th17) / nrow(df_cd4_nn),
                    `B cell Naive` = nrow(df_bcell_naive) / nrow(df_bcell),
                    `B cell Mem Sw` = nrow(df_bcell_mem_sw) / nrow(df_bcell),
                    `B cell Mem NotSw` = nrow(df_bcell_mem_nsw) / nrow(df_bcell),
                    `B cell Mem` = nrow(df_bcell_mem) / nrow(df_bcell),
                    `NK cell late` = nrow(df_nkcell_late) / nrow(df_nkcell),
                    `Monocyte Classical` = nrow(df_mono_cl) / (nrow(df_mono_cl) + nrow(df_mono_ncl)),
                    `Monocyte Nonclassical` = nrow(df_mono_ncl) / (nrow(df_mono_cl) + nrow(df_mono_ncl)),
                    `mDC` = nrow(df_mdc) / n_tot,
                    ########
                    `CD4 CD57hi` = length(which(df_cd4_all$cd57hi)) / nrow(df_cd4),
                    `CD4 CD57hi out of EMRA` = get_ratio(df_cd4_all, "emra", "cd57hi"),
                    `CD4 CD57hi out of EM1` = get_ratio(df_cd4_all, "em1", "cd57hi"),
                    `CD4 CD57hi out of EM3` = get_ratio(df_cd4_all, "em3", "cd57hi"),
                    `CD4 CD57hi out of Th1` = get_ratio(df_cd4_all, "th1", "cd57hi"),
                    `CD4 CD57hi out of Tfh` = get_ratio(df_cd4_all, "tfh", "cd57hi"),
                    `CD4 CD57hi out of Treg` = get_ratio(df_cd4_all, "treg", "cd57hi"),
                    `CD4 CD57hi out of Act` = get_ratio(df_cd4_all, "act", "cd57hi"),

                    `CD4 Act out of EMRA` = get_ratio(df_cd4_all, "emra", "act"),
                    `CD4 Act out of EM1` = get_ratio(df_cd4_all, "em1", "act"),
                    `CD4 Act out of EM3` = get_ratio(df_cd4_all, "em3", "act"),
                    `CD4 Act out of CM` = get_ratio(df_cd4_all, "cm", "act"),
                    `CD4 Act out of Th1` = get_ratio(df_cd4_all, "th1", "act"),
                    `CD4 Act out of Tfh` = get_ratio(df_cd4_all, "tfh", "act"),
                    `CD4 Act out of Treg` = get_ratio(df_cd4_all, "treg", "act"),

                    `CD4 EMRA out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "emra"),
                    `CD4 EM1 out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "em1"),
                    `CD4 EM3 out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "em3"),
                    `CD4 Th1 out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "th1"),
                    `CD4 Tfh out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "tfh"),
                    `CD4 Treg out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "treg"),
                    `CD4 Act out of CD57hi` = get_ratio(df_cd4_all, "cd57hi", "act"),

                    `CD4 EMRA out of Act` = get_ratio(df_cd4_all, "act", "emra"),
                    `CD4 EM1 out of Act` = get_ratio(df_cd4_all, "act", "em1"),
                    `CD4 EM3 out of Act` = get_ratio(df_cd4_all, "act", "em3"),
                    `CD4 CM out of Act` = get_ratio(df_cd4_all, "act", "cm"),
                    `CD4 Th1 out of Act` = get_ratio(df_cd4_all, "act", "th1"),
                    `CD4 Tfh out of Act` = get_ratio(df_cd4_all, "act", "tfh"),
                    `CD4 Treg out of Act` = get_ratio(df_cd4_all, "act", "treg"),

                    `CD8 CD57hi` = length(which(df_cd8_all$cd57hi)) / nrow(df_cd8),
                    `CD8 CD57hi out of EMRA` = get_ratio(df_cd8_all, "emra", "cd57hi"),
                    `CD8 CD57hi out of EM1` = get_ratio(df_cd8_all, "em1", "cd57hi"),
                    `CD8 CD57hi out of EM3` = get_ratio(df_cd8_all, "em3", "cd57hi"),
                    `CD8 CD57hi out of Act` = get_ratio(df_cd8_all, "act", "cd57hi"),

                    `CD8 Act out of EMRA` = get_ratio(df_cd8_all, "emra", "act"),
                    `CD8 Act out of EM1` = get_ratio(df_cd8_all, "em1", "act"),
                    `CD8 Act out of EM3` = get_ratio(df_cd8_all, "em3", "act"),
                    `CD8 Act out of CM` = get_ratio(df_cd8_all, "cm", "act"),

                    `CD8 EMRA out of CD57hi` = get_ratio(df_cd8_all, "cd57hi", "emra"),
                    `CD8 EM1 out of CD57hi` = get_ratio(df_cd8_all, "cd57hi", "em1"),
                    `CD8 EM3 out of CD57hi` = get_ratio(df_cd8_all, "cd57hi", "em3"),
                    `CD8 Act out of CD57hi` = get_ratio(df_cd8_all, "cd57hi", "act"),

                    `CD8 EMRA out of Act` = get_ratio(df_cd8_all, "act", "emra"),
                    `CD8 EM1 out of Act` = get_ratio(df_cd8_all, "act", "em1"),
                    `CD8 EM3 out of Act` = get_ratio(df_cd8_all, "act", "em3"),
                    `CD8 CM out of Act` = get_ratio(df_cd8_all, "act", "cm"))

  if (save_plots) {

    # params <- c("CD4", "CD161")
    # p1 <- gg_flow_thresh(df_ab, params, thresh) + ggtitle("ab T cells")
    # p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    # p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    # p <- p1 + p2 + p3 + plot_layout(ncol=2)
    # ggsave(p, filename=paste0(dir_out, "gating/mait_nkt/", fn, ".png"), width=12, height=10)

    params <- c("CD45RA", "CD27")
    p3 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p1 <- gg_flow_thresh(df_cd4, params, thresh, thresh_neut, thresh_mono) + ggtitle("CD4 T cells")
    p2 <- gg_flow_thresh(df_cd8, params, thresh, thresh_neut, thresh_mono) + ggtitle("CD8 T cells")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tcell_mem/", fn, ".png"), width=12, height=10)

    params <- c("CD27", "CD197")
    p1 <- gg_flow_thresh(df_cd4_cd45ralo, params, thresh, thresh_neut, thresh_mono) + ggtitle("CD4 T cells")
    p2 <- gg_flow_thresh(df_cd8_cd45ralo, params, thresh, thresh_neut, thresh_mono) + ggtitle("CD8 T cells")
    p3 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tcell_em_cm/", fn, ".png"), width=12, height=10)

    params <- c("CD38", "HLA-DR")
    p1 <- gg_flow_gate(df_cd4, params, gate_act) + ggtitle("CD4 T cells")
    p2 <- gg_flow_gate(df_cd8, params, gate_act) + ggtitle("CD8 T cells")
    p3 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_myel, params, thresh, thresh_neut, thresh_mono) + ggtitle("Myeloid")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tcell_act/", fn, ".png"), width=12, height=10)

    params <- c("CD27", "IgD")
    p1 <- gg_flow_thresh(df_bcell, params, thresh, thresh_neut, thresh_mono) + ggtitle("B cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_myel, params, thresh, thresh_neut, thresh_mono) + ggtitle("Myeloid")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/bcell_mem/", fn, ".png"), width=12, height=10)

    params <- c("CD25", "CD127")
    p1 <- gg_flow_gate(df_cd4_nn, params, gate_treg) + ggtitle("Non-naive CD4 T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_myel, params, thresh, thresh_neut, thresh_mono) + ggtitle("Myeloid")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/treg/", fn, ".png"), width=12, height=10)

    params <- c("CD185", "CD4")
    p1 <- gg_flow_thresh(df_cd4_nn, params, thresh, thresh_neut, thresh_mono) + ggtitle("Non-naive CD4 T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_myel, params, thresh, thresh_neut, thresh_mono) + ggtitle("Myeloid")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tfh/", fn, ".png"), width=12, height=10)

    params <- c("CD196", "CD183")
    p1 <- gg_flow_thresh(df_cd4_cd185n, params, thresh, thresh_neut, thresh_mono) + ggtitle("Non-naive CD185lo CD4 T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_myel, params, thresh, thresh_neut, thresh_mono) + ggtitle("Myeloid")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/cd4_func/", fn, ".png"), width=12, height=10)

    params <- c("CD57", "CD56")
    p1 <- gg_flow_thresh(df_nkcell, params, thresh, thresh_neut, thresh_mono) + ggtitle("NK cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_cd4, params, thresh, thresh_neut, thresh_mono) + ggtitle("CD4 T cells")
    p4 <- gg_flow_thresh(df_cd8, params, thresh, thresh_neut, thresh_mono) + ggtitle("CD8 T cells")
    p <- p1 + p2 + p3 + p4 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/nk_late_early/", fn, ".png"), width=12, height=10)

    params <- c("CD14", "CD38")
    p1 <- gg_flow_thresh(df_myel, params, thresh, thresh_neut, thresh_mono) + ggtitle("Myeloid")
    p2 <- gg_flow_thresh(df_neut, params, thresh, thresh_neut, thresh_mono) + ggtitle("Neutrophils")
    p <- p1 + p2 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/mono/", fn, ".png"), width=12, height=5)
  }

  return(df_feat)
}


get_ratio <- function(df, col1, col2) {
  num <- length(which(df[[col1]] & df[[col2]]))
  denom <- length(which(df[[col1]]))
  return(num/denom)
}


apply_gate <- function(df, params, gate) {
  in_gate <- point.in.polygon(df[[params[1]]],
                                  df[[params[2]]],
                                  gate[[params[1]]],
                                  gate[[params[2]]])
  return(df[which(in_gate>0),])
}


get_treg_gate <- function(df, params, thresh) {

  t1 <- thresh[params[1]]
  t2 <- thresh[params[2]]
  med <- median(df[[params[2]]])
  M <- max(df[[params[1]]])


  x1 <- df %>% pull(params[1])
  x2 <- df %>% pull(params[2])
  x <- x2-x1
  inter <- deGate(x, NA)

  gate <- tibble(var1 = c(t1, med-inter, M, M, t1, t1),
                 var2 = c(t1+inter, med, med, 0, 0, t1+inter))
  names(gate) <- params
  return(gate)
}

get_act_gate <- function(df, params, thresh) {
  t1 <- thresh[params[1]]
  t2 <- thresh[params[2]]
  M1 <- max(df[[params[1]]])
  M2 <- max(df[[params[2]]])

  gate <- tibble(var1 = c(M1, 0.75*t1, 0.75*t1, t1, M1, M1),
                 var2 = c(M2, M2, t2, 0.6*t2, 0.6*t2, M2))
  names(gate) <- params
  return(gate)
}


estimate_distributions <- function(cell_type, df, fn,
                                   channels, cell_types) {

  tab <- table(cell_type)
  to_use <- intersect(cell_types, names(which(tab > 100)))
  to_use <- union(to_use, "all")

  df_kdes <- lapply(to_use, function(ct) {

    if (ct == "all") {
      cells <- which(!grepl("debris|doublet", cell_type))
    } else {
      cells <- which(cell_type == ct)
    }

    lapply(channels, function(ch) {
      x <- df[[ch]][cells]
      kde <- bkde(x, range.x=c(-1,8), gridsize = 101L)

      return(tibble(file=fn, cell_type=ct, channel=ch,
                    expression = kde$x,
                    density = kde$y / max(kde$y)))
    }) %>% do.call(what=rbind)
  }) %>% do.call(what=rbind)

  return(df_kdes)
}


##### from Mito

get_one_th_fd <- function(df, cell_type, ct, channel) {
  x <- df %>% filter(grepl(ct, cell_type)) %>% pull(channel)
  deGate(x, channel)
}

get_two_th_fd <- function(df, cell_type, ct, channels, position) {
  df_ct <- df %>% filter(grepl(ct, cell_type)) %>% select(all_of(channels))
  x <- df_ct %>% as.matrix() %>% flowFrame()
  fd <- flowDensity(x, channels, position)

  if (position[1]) {
    t1 <- quantile(df_ct[[channels[1]]][fd@index], 0.01)
  } else {
    t1 <- quantile(df_ct[[channels[1]]][fd@index], 0.99)
  }

  if (position[2]) {
    t2 <- quantile(df_ct[[channels[2]]][fd@index], 0.01)
  } else {
    t2 <- quantile(df_ct[[channels[2]]][fd@index], 0.99)
  }

  return(c(t1,t2))
}

get_thresholds_density <- function(df, cell_type) {
  thresholds <- setNames(rep(0, length(df)), names(df))
  thresholds["CD57"] <- get_one_th_fd(df, cell_type, "nkcell", "CD57")

  thresholds["CD45RA"] <- get_one_th_fd(df, cell_type, "tcell_cd4", "CD45RA")
  thresholds["CD27"] <- get_one_th_fd(df, cell_type, "tcell_cd4", "CD27")
  thresholds["CD197"] <- get_one_th_fd(df, cell_type, "tcell_cd4", "CD197")

  thresholds["CD185"] <- get_one_th_fd(df, cell_type, "tcell_cd4", "CD185")
  thresholds["CD183"] <- get_one_th_fd(df, cell_type, "tcell_cd4", "CD183")
  thresholds["CD196"] <- get_one_th_fd(df, cell_type, "tcell_cd4", "CD196")

  thresholds["CD14"] <- get_one_th_fd(df, cell_type, "myeloid", "CD14")
  thresholds["CD38"] <- get_one_th_fd(df, cell_type, "myeloid", "CD38")

  thresholds["IgD"] <- get_one_th_fd(df, cell_type, "bcell", "IgD")

  thresholds[c("CD14", "CD38")] <- get_two_th_fd(df, cell_type, "myeloid",
                                                 c("CD14", "CD38"), c(TRUE,TRUE))
  x <- df %>% filter(cell_type == "tcell_cd8" & CD38 > thresholds["CD38"]) %>%
    pull("HLA-DR")
  thresholds["HLA-DR"] <- deGate(x,"HLA-DR")

  x <- df %>% filter(cell_type == "tcell_cd4") %>%
    pull("CD25")
  thresholds["CD25"] <- deGate(x,"CD25")

  return(thresholds)
}


get_hull <- function(data, context, gate_col) {
  df <- data %>%
    filter(.data[[gate_col]] %in% context$this_descendants) %>%
    select(all_of(c(context$x_col, context$y_col)))
  idx <- chull(df)
  hull <- df[c(idx, idx[1]),]
  return(hull)
}


plot_gate <- function(data, context, hull) {

  x_col <- context$x_col
  y_col <- context$y_col

  hull[[x_col]] <- pmax(pmin(hull[[x_col]], max(data[[x_col]])),
                        min(data[[x_col]]))
  hull[[y_col]] <- pmax(pmin(hull[[y_col]], max(data[[y_col]])),
                        min(data[[y_col]]))

  p <- gg_flow(data, params = c(x_col, y_col)) +
    geom_path(data=hull, color="red") +
    ggtitle(paste0("Gate: ", context$this_name,
                   "; parent: ", context$parent_name)) +
    theme_bw(base_size=16)
  return(p)
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

  non_naive <- if_else(!is.na(gate_hierarchy %>% filter(GateName==gate_name) %>% pull(NonNaive)),
                       TRUE, FALSE)

  context <- list(this_name = gate_name,
                  this_descendants = this_descendants,
                  parent_name = parent_name,
                  parent_descendants = parent_descendants,
                  x_col = x_col,
                  y_col = y_col,
                  x_loc = x_loc,
                  y_loc = y_loc,
                  non_naive = non_naive)
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

get_fd_gate <- function(df, context) {
  ff <- flowFrame(as.matrix(df))
  fd <- flowDensity(ff, channels = c(context$x_col, context$y_col),
              position = c(context$x_loc, context$y_loc)=="hi",
              ellip.gate=FALSE)
  return(as_tibble(fd@filter))
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

apply_gate2 <- function(df, gate) {
  x_col <- names(gate)[1]
  y_col <- names(gate)[2]

  pip <- point.in.polygon(df %>% pull(x_col),
                          df %>% pull(y_col),
                          gate %>% pull(x_col),
                          gate %>% pull(y_col))
  return(pip > 0)
}

