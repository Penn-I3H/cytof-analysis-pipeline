read_data <- function(path) {
  ff <- read.FCS(path, truncate_max_range = FALSE)
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

apply_cleanup_gate <- function(df, cutoffs, marker) {
  df %>%
    filter(.data[[marker]] >= cutoffs[[marker]][1] &
             .data[[marker]] <= cutoffs[[marker]][2])
}


pregate_data <- function(df, fn, dir_out, plot=FALSE) {

  cutoffs <- list("140Ce_Bead" = c(0,asinh(500/5)),
                  "Residual" = c(asinh(20/5), asinh(400/5)),
                  "Center" = c(asinh(300/5), asinh(4000/5)),
                  "Offset" = c(asinh(30/5), asinh(400/5)),
                  "Width" = c(asinh(30/5), asinh(300/5)),
                  "Live" = c(0, asinh(80/5)))
  markers <- names(cutoffs)

  df1 <- apply_cleanup_gate(df, cutoffs, markers[1])
  df2 <- apply_cleanup_gate(df1, cutoffs, markers[2])
  df3 <- apply_cleanup_gate(df2, cutoffs, markers[3])
  df4 <- apply_cleanup_gate(df3, cutoffs, markers[4])
  df5 <- apply_cleanup_gate(df4, cutoffs, markers[5])
  df6 <- apply_cleanup_gate(df5, cutoffs, markers[6])

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

  markers_keep <- names(df)[!grepl("Time|Bead|Live|Center|Offset|Residual|Width", names(df))]

  df_clean <- df6 %>%
    select(all_of(markers_keep))

  return(df_clean)
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


gg_flow_thresh <- function(df, params, thresh, nbin=100, method="bin",
                           xlim=c(-0.1,7), ylim=c(-0.1,7)) {
  gg_flow(df, params, nbin=nbin, method=method, xlim=xlim, ylim=ylim) +
    geom_hline(yintercept=thresh[params[2]], color="red") +
    geom_vline(xintercept=thresh[params[1]], color="red")
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
                         "debris", "agg")

  p1 <- plot_umap_clust(df, "ct", fn, pal_paired)
  p2 <- plot_umap_channel(df, "CD3")
  p3 <- plot_umap_channel(df, "CD19")
  p4 <- plot_umap_channel(df, "CD66b")
  p5 <- plot_umap_channel(df, "CD123")
  p6 <- plot_umap_channel(df, "CD294")
  p7 <- plot_umap_channel(df, "CD56")
  p8 <- plot_umap_channel(df, "CD11c")
  p9 <- plot_umap_channel(df, "DNA1")

  p <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)
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


# separate_doublets_feb <- function(clustering, df, cols) {
#   clustering_new <- as.character(clustering)
#   tab <- table(clustering)
#   candidates <- setdiff(names(which(tab>200)), grep("agg|debris", names(tab), value=TRUE))
#
#   to_check <- c("tcell", "neutrophil", "bcell")
#   to_check <- to_check[which(!paste0("agg_", to_check) %in% levels(clustering))]
#
#   if (length(to_check)==0)
#     return(clustering)
#
#   candidates <- candidates[grep(paste(to_check, collapse="|"), candidates)]
#
#   for (target in candidates) {
#     cells <- which(clustering==target)
#
#     m10 <- median(df[["DNA1"]][cells])
#     m20 <- median(df[["DNA2"]][cells])
#     iqr10 <- IQR(df[["DNA1"]][cells])
#     iqr20 <- IQR(df[["DNA2"]][cells])
#
#     dna <- 0.5*(df[["DNA1"]][cells] + df[["DNA2"]][cells])
#     m <- median(dna)
#     iqr <- IQR(dna)
#     doub <- which((dna > m+2*iqr) & (sinh(dna) > 1.9*sinh(m)))
#
#     if (length(doub) > 50) {
#       bh <- test_pure_doub(df, cells_sing=setdiff(cells, cells[doub]),
#                            cells_doub=cells[doub], cols=cols)
#       if (bh > 0.2)
#         clustering_new[cells[doub]] <- paste0("agg_", target)
#     }
#   }
#
#   return(as.factor(clustering_new))
# }
#
#
# test_pure_doub <- function(df, cells_sing, cells_doub, cols) {
#
#   bhatt <- lapply(cols, function(ch) {
#     m <- min(df[[ch]])
#     M <- max(df[[ch]])
#
#     kde_sing <- bkde(df[[ch]][cells_sing], range.x = c(m,M), gridsize=101L)
#     kde_doub <- bkde(df[[ch]][cells_doub], range.x = c(m,M), gridsize=101L)
#
#     intensity <- kde_sing$x
#     dist_sing <- pmax( kde_sing$y / sum(kde_sing$y) , 0)
#     dist_doub <- pmax( kde_doub$y / sum(kde_doub$y) , 0)
#
#     f1 <- approxfun(intensity, dist_sing, yleft=0, yright=0)
#     scal <- if_else(ch=="length", "linear", "asinh")
#     dist_ref <- conv_sum(f1,f1,intensity, scal=scal)
#
#     bh <- compute_bhatt(dist_ref, dist_doub)
#
#     return(bh)
#   }) %>% Reduce(f=prod)
#
#   return(bhatt)
# }


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
    reshape2::melt() %>%
    filter(Var1 <= Var2) %>%
    mutate(value = replace_na(value, 0)) %>%
    filter(value >= min_bhatt)
  # filter(value >= 0.88^length(cols))

  gr <- igraph::graph_from_edgelist(bhatt_summ %>%
                                      select(Var1,Var2) %>%
                                      as.matrix(),
                                    directed=FALSE)
  comp <- components(gr)$membership
  # names(comp) <- levels(clustering)

  cl_merged <- comp[clustering] %>% unname() %>% as.factor()
  return(cl_merged)
}

detect_doublets_apr <- function(df, cols, cell_type, thresh=5) {

  x0 <- as.matrix(df)[,setdiff(cols, "Event_length")]

  set.seed(0)
  n_doub <- floor(nrow(x0)/10)
  sel_for_aug <- which(cell_type != "debris")
  sel1 <- sample(sel_for_aug, n_doub)
  sel2 <- sample(sel_for_aug, n_doub)
  x_doub <- asinh(sinh(x0[sel1,]) + sinh(x0[sel2,]))

  x_aug <- rbind(x_doub, x0[sel_for_aug,])

  all_knn <- RcppHNSW::hnsw_knn(x_aug, k=15, distance= 'l2',
                                n_threads=1, M=48)
  n_nb_doub <- apply(all_knn$idx, 1, function(row) length(which(row<n_doub)))


  sim_lab <- if_else(cell_type[sel1] < cell_type[sel2],
                     paste0("agg_", cell_type[sel1], "_", cell_type[sel2]),
                     paste0("agg_", cell_type[sel2], "_", cell_type[sel1]))
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


# get_labels_with_doublets_mar <- function(clustering, df, gr, defs, centroids_sc) {
#
#   labels <- levels(clustering)
#
#   el <- as_edgelist(gr)
#   df_el <- tibble(source = as.integer(el[,1]),
#                   target = as.integer(el[,2]),
#                   weight = E(gr)$weight)
#
#   is_doub <- labels %in% labels[unique(df_el$target)]
#
#   lab_not_doub <- labels[which(!is_doub)]
#   ev_not_doub <- which(clustering %in% lab_not_doub)
#
#   mat <- data.matrix(df[ev_not_doub,]) %>% scale()
#
#   tmp <- as.factor(as.character(clustering[ev_not_doub]))
#   centroids <- get_medians(mat, tmp)
#   labels <- rep("agg", length(labels))
#   # labels[which(!is_doub)] <- label_clusters_score_opt(centroids[lab_not_doub,], defs, mdipa_main = FALSE)
#   labels[which(!is_doub)] <- label_clusters_score_opt(centroids_sc[which(!is_doub),], defs, mdipa_main = FALSE)
#
#
#   df_el <- df_el %>%
#     mutate(source_name = labels[source]) %>%
#     arrange(target, -weight)
#
#   ts <- as.integer(V(gr)$name[as.integer(topological.sort(gr))])
#   for (i in ts) {
#     if (labels[i]=="agg") {
#       components <- df_el %>%
#         filter(target==i & !(grepl("debris|agg", source_name))) %>%
#         pull(source_name) %>%
#         unique() %>%
#         sort()
#
#       labels[i] <- paste(c("agg", components), collapse="_")
#       if (labels[i] == "agg")
#         labels[i] <- "debris"
#     }
#   }
#
#   return(make.unique(labels))
# }



# plot_graph_suga <- function(gr) {
#   layout <- create_layout(gr, layout = "sugiyama")
#   ggraph(layout) +
#     geom_edge_fan(aes(alpha = weight),
#                   arrow = arrow(length = unit(4, 'mm')),
#                   end_cap = circle(3, 'mm')) +
#     geom_node_label(aes(label=name), size=1) +
#     theme_graph(base_family = "sans") +
#     theme(text = element_text(size = 18))
# }


# detect_doublets <- function(df, cols, clustering) {
#
#   df_kdes <- get_kdes(df, cols, clustering) %>%
#     pivot_wider(names_from="cluster", values_from="density")
#
#   tab <- table(clustering)
#   med_len_dna <- get_medians(df[,c("length","DNA1")], clustering)
#
#   system.time(bhatt_summ <- lapply(names(tab), get_bhatt_target_fast,
#                                    df_kdes=df_kdes, tab=tab,
#                                    med_len_dna=med_len_dna, cols=cols) %>%
#                 do.call(what=rbind))
#
#   bhatt_graph <- bhatt_summ %>%
#     filter(bhatt > 0.2) %>%
#     pivot_longer(all_of(c("cl1","cl2")))
#
#   gr <- graph_from_edgelist(bhatt_graph %>% select(value, target) %>% as.matrix())
#   E(gr)$weight <- bhatt_graph$bhatt
#
#   return(gr)
# }


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

get_thresh_gauss <- function(df_neut, df_myel, df_cd4, df_bcell) {
  thresh <- sapply(names(df_neut), function(col) one_thresh_gaussian(df_neut[[col]]))
  thresh_myel <- sapply(names(df_myel), function(col) one_thresh_gaussian(df_myel[[col]]))
  thresh_cd4 <- -sapply(names(df_cd4), function(col) one_thresh_gaussian(-df_cd4[[col]], sym_zero = FALSE))
  thresh_bcell <- -sapply(names(df_bcell), function(col) one_thresh_gaussian(-df_bcell[[col]], sym_zero = FALSE))

  mark_m <- c("IgD", "TCRgd", "CD183")
  thresh[mark_m] <- thresh_myel[mark_m]

  mark_cd4 <- c("CD4")
  thresh[mark_cd4] <- thresh_cd4[mark_cd4]

  mark_bcell <- c("HLA-DR")
  thresh[mark_bcell] <- thresh_bcell[mark_bcell]

  thresh["CD196"] <- 2
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
  df_neut <- df %>% filter(cell_type=="neutrophil")
  df_mono <- df %>% filter(grepl("monocyte", cell_type))
  df_mono_cl <- df %>% filter(cell_type=="monocyte_classical")
  df_cd4 <- df %>% filter(cell_type=="tcell_cd4")
  df_cd8 <- df %>% filter(cell_type=="tcell_cd8")
  df_gd <- df %>% filter(cell_type=="tcell_gd")
  df_bcell <- df %>% filter(cell_type=="bcell")
  df_nkcell <- df %>% filter(cell_type=="nkcell")

  df_kdes <- get_kdes(df, names(df), as.factor(cell_type))

  thresh <- get_thresh_gauss(df_neut, df_mono_cl, df_cd4, df_bcell)

  ch_care <- c("CD4", "CD161", "CD45RA", "CD27", "CD197",
               "CD25", "CD127", "CD185", "CD183", "CD196",
               "CD38", "HLA-DR", "IgD", "CD57")
  p <- plot_kdes(df_kdes %>% filter(channel %in% ch_care),
                 c("neutrophil", "tcell_cd4", "tcell_cd8", "bcell", "monocyte_classical", "nkcell")) +
    geom_vline(data=tibble(channel=names(thresh), intensity=unname(thresh)) %>% filter(channel %in% ch_care),
               mapping=aes(xintercept=intensity), linetype="dashed")
  ggsave(p, filename=paste0(dir_out, "gating/kdes_thresh/", fn, ".png"), width=12, height=10)

  n_tot <- length(which(!grepl("debris|agg", cell_type)))
  n_tcell <- length(which(grepl("tcell", cell_type) & !grepl("agg", cell_type)))
  df_feat_maj <- tibble(file = fn,
                        `T cell` = n_tcell / n_tot,
                        `T cell CD4` = nrow(df_cd4) / n_tcell,
                        `T cell CD8` = nrow(df_cd8) / n_tcell,
                        `T cell gd` = nrow(df_gd) / n_tcell,
                        `B cell` = nrow(df_bcell) / n_tot,
                        `NK cell` = nrow(df_nkcell) / n_tot,
                        `Monocyte` = nrow(df_mono) / n_tot,
                        `Monocyte Classical` = length(which(cell_type=="monocyte_classical")) / nrow(df_mono),
                        `Monocyte Nonclassical` = length(which(cell_type=="monocyte_nonclassical")) / nrow(df_mono),
                        `mdc` = length(which(cell_type=="mdc")) / n_tot,
                        `pdc` = length(which(cell_type=="pdc")) / n_tot,
                        `Plasmablast` = length(which(cell_type=="plasmablast")) / n_tot,
                        `Neutrophil` = nrow(df_neut) / n_tot,
                        `Eosinophil` = length(which(cell_type=="eosinophil")) / n_tot,
                        `Basophil` = length(which(cell_type=="basophil")) / n_tot)
  write_csv(df_feat_maj, file=paste0(dir_out, "/feat_major/feat_major_", fn, ".csv"))

  df_feat_adaptive <- gate_and_plot(df, cell_type, thresh, dir_out, fn, save_plots=TRUE)
  write_csv(df_feat_adaptive, file=paste0(dir_out, "/feat_adaptive/feat_adaptive_", fn, ".csv"))

}


gate_and_plot <- function(df, cell_type, thresh, dir_out, fn, save_plots=TRUE) {

  df_ab <- df %>% filter(grepl("cd4|cd8|dn|dp", cell_type))

  df_neut <- df %>% filter(cell_type=="neutrophil")
  df_mono <- df %>% filter(grepl("monocyte", cell_type))
  df_mono_cl <- df %>% filter(cell_type=="monocyte_classical")
  df_cd4 <- df %>% filter(cell_type=="tcell_cd4")
  df_cd8 <- df %>% filter(cell_type=="tcell_cd8")
  df_gd <- df %>% filter(cell_type=="tcell_gd")
  df_bcell <- df %>% filter(cell_type=="bcell")
  df_nkcell <- df %>% filter(cell_type=="nkcell")

  df_mait <- df_ab %>% filter(CD4<thresh["CD4"] & CD161>thresh["CD161"])

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

  df_feat <- tibble(file = fn,
                    `MAIT/NKT` = nrow(df_mait) / nrow(df_ab),
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
                    `NK cell late` = nrow(df_nkcell_late) / nrow(df_nkcell))

  if (save_plots) {

    params <- c("CD4", "CD161")
    p1 <- gg_flow_thresh(df_ab, params, thresh) + ggtitle("ab T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/mait_nkt/", fn, ".png"), width=12, height=10)

    params <- c("CD45RA", "CD27")
    p3 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p1 <- gg_flow_thresh(df_cd4, params, thresh) + ggtitle("CD4 T cells")
    p2 <- gg_flow_thresh(df_cd8, params, thresh) + ggtitle("CD8 T cells")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tcell_mem/", fn, ".png"), width=12, height=10)

    params <- c("CD27", "CD197")
    p1 <- gg_flow_thresh(df_cd4_cd45ralo, params, thresh) + ggtitle("CD4 T cells")
    p2 <- gg_flow_thresh(df_cd8_cd45ralo, params, thresh) + ggtitle("CD8 T cells")
    p3 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tcell_em_cm/", fn, ".png"), width=12, height=10)

    params <- c("CD38", "HLA-DR")
    p1 <- gg_flow_gate(df_cd4, params, gate_act) + ggtitle("CD4 T cells")
    p2 <- gg_flow_gate(df_cd8, params, gate_act) + ggtitle("CD8 T cells")
    p3 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_bcell, params, thresh) + ggtitle("B cells")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tcell_act/", fn, ".png"), width=12, height=10)

    params <- c("CD27", "IgD")
    p1 <- gg_flow_thresh(df_bcell, params, thresh) + ggtitle("B cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/bcell_mem/", fn, ".png"), width=12, height=10)

    params <- c("CD25", "CD127")
    p1 <- gg_flow_gate(df_cd4_nn, params, gate_treg) + ggtitle("Non-naive CD4 T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/treg/", fn, ".png"), width=12, height=10)

    params <- c("CD185", "CD4")
    p1 <- gg_flow_thresh(df_cd4_nn, params, thresh) + ggtitle("Non-naive CD4 T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/tfh/", fn, ".png"), width=12, height=10)

    params <- c("CD196", "CD183")
    p1 <- gg_flow_thresh(df_cd4_cd185n, params, thresh) + ggtitle("Non-naive CD185lo CD4 T cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/cd4_func/", fn, ".png"), width=12, height=10)

    params <- c("CD57", "CD56")
    p1 <- gg_flow_thresh(df_nkcell, params, thresh) + ggtitle("NK cells")
    p2 <- gg_flow_thresh(df_neut, params, thresh) + ggtitle("Neutrophils")
    p3 <- gg_flow_thresh(df_mono, params, thresh) + ggtitle("Monocytes")
    p <- p1 + p2 + p3 + plot_layout(ncol=2)
    ggsave(p, filename=paste0(dir_out, "gating/nk_late_early/", fn, ".png"), width=12, height=10)

  }

  return(df_feat)
}


apply_gate <- function(df, params, gate) {
  in_gate <- sp::point.in.polygon(df[[params[1]]],
                                  df[[params[2]]],
                                  gate[[params[1]]],
                                  gate[[params[2]]])
  return(df[which(in_gate>0),])
}


get_treg_gate <- function(df, params, thresh) {
  t1 <- thresh[params[1]]
  t2 <- thresh[params[2]]
  M <- max(df[[params[1]]])

  gate <- tibble(var1 = c(t1, t1/2, t1/4, t1/4, M, M, t1),
                 var2 = c(t2, t2, 0.75*t2, 0, 0, t2, t2))
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
      cells <- which(!grepl("debris|agg", cell_type))
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




