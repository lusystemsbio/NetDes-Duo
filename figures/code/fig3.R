library(Seurat)
library(ggplot2)
library(dplyr)
library(visNetwork)

DATA_DIR <- "../data"
OUT_DIR  <- "../plots"


# figure 3A
nnb <- readRDS(file.path(DATA_DIR, "nnb_pca.rds"))
ncb <- readRDS(file.path(DATA_DIR, "ncb_pca.rds"))

#PCA plot by tissue
sdev <- nnb[["pca"]]@stdev
pct.var <- (sdev^2 / sum(sdev^2)) * 100
pc1_pct = formatC(signif(pct.var[[1]], 3), format = "fg", digits = 3, flag = "#")
pc2_pct = formatC(signif(pct.var[[2]], 3), format = "fg", digits = 3, flag = "#")

p <- PCAPlot(nnb, group.by = "MULTI_ID",pt.size = 0.25) +
  labs(title = NULL, x = paste0("PC1 (",pc1_pct, "%)")
       , y = paste0("PC2 (", pc2_pct, "%)"), color = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.title = element_text(size = 22.5),
        axis.text  = element_text(size = 22.5, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  scale_x_reverse() + scale_y_reverse(limits = c(5.5, -7))
ggsave(
  filename = file.path(OUT_DIR, "fig3A_naive.pdf"),
  plot = p,
  width = 300,
  height = 150,
  units = "mm"
)

sdev <- ncb[["pca"]]@stdev
pct.var <- (sdev^2 / sum(sdev^2)) * 100
pc1_pct = formatC(signif(pct.var[[1]], 3), format = "fg", digits = 3, flag = "#")
pc2_pct = formatC(signif(pct.var[[2]], 3), format = "fg", digits = 3, flag = "#")

p <- PCAPlot(ncb, group.by = "MULTI_ID",pt.size = 0.25) +
  labs(title = NULL, x = paste0("PC1 (",pc1_pct, "%)")
       , y = paste0("PC2 (", pc2_pct, "%)"), color = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.title = element_text(size = 22.5),
        axis.text  = element_text(size = 22.5, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  scale_x_reverse() + scale_y_reverse(limits = c(5.5, -7))
ggsave(
  filename = file.path(OUT_DIR, "fig3A_cancer.pdf"),
  plot = p,
  width = 300,
  height = 150,
  units = "mm"
)


# figure 3B
plot_clusters_mean_scaled <- function(coord_x, coord_y, gene_clusters, xp, cutoff = 20,
                                      title = "Mean-Scaled Trajectory Clusters",
                                      ylim = c(0, 5), n_col = 5) {

  cluster_counts <- sort(table(gene_clusters[gene_clusters != ""]), decreasing = TRUE)
  keep_clusters <- names(cluster_counts)[cluster_counts > cutoff]

  n_row <- ceiling(length(keep_clusters) / n_col)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_row, n_col), oma = c(0, 0, 3, 0), mar = c(2.2, 2.2, 1.5, 0.8),
      mgp = c(1.2, 0.4, 0), tcl = -0.25, cex.axis = 1.5, cex.main = 1.5)

  grid_x <- seq(0, 1, length.out = 80)

  for (cl in keep_clusters) {
    idx <- which(gene_clusters == cl)
    plot(NA, xlim = c(0, 1), ylim = ylim, xlab = "", ylab = "",
         main = paste0(cl, " (n = ", as.character(length(idx)), ")"), type = "n")

    y_mat <- matrix(NA, nrow = length(idx), ncol = length(grid_x))
    print("here")
    for (j in seq_along(idx)) {
      i <- idx[j]
      x <- as.numeric(coord_x[[i]])
      y <- as.numeric(coord_y[[i]])

      m <- mean(y, na.rm = TRUE)
      y_scaled <- y / m

      lines(x, y_scaled, type = "o",
            pch = 1, cex = 0.5, lwd = 0.5, col = rgb(0, 0, 0, 0.1))

      y_mat[j, ] <- approx(x, y_scaled, xout = grid_x, rule = 2)$y
    }

    avg_y <- colMeans(y_mat, na.rm = TRUE)
    lines(grid_x, avg_y, col = "red", lwd = 2)
  }

  mtext(title, side = 3, outer = TRUE, line = 1, cex = 1.4)
}

# --- naive ---
traj_clusters <- readRDS(file.path(DATA_DIR, "naive_traj_clusters.rds"))

#take the automatic clusters and manually clean them (seperate di and dfi into 2 clusters)
tc <- traj_clusters
cluster_labels <- as.character(tc$cluster)

cluster1  <- c("did")
cluster2  <- c("d", "df", "difd", "dif")
cluster3  <- c("id")
cluster4  <- c("fi", "i", "fdi")
cluster67 <- c("di", "dfi")

# initialize manual cluster column
tc$manual_cluster <- NA

# combine into 6 broad clusters first
tc$manual_cluster[tc$cluster %in% cluster1]  <- "Cluster 2"
tc$manual_cluster[tc$cluster %in% cluster2]  <- "Cluster 1"
tc$manual_cluster[tc$cluster %in% cluster3]  <- "Cluster 3"
tc$manual_cluster[tc$cluster %in% cluster4]  <- "Cluster 5"
tc$manual_cluster[tc$cluster %in% cluster67] <- "cluster23_di_dfi"

idx67 <- which(tc$manual_cluster == "cluster23_di_dfi")

for (i in idx67) {
  x <- as.numeric(tc$x[[i]])
  y <- as.numeric(tc$y[[i]])

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  min_pos <- x[which.min(y)]

  if (min_pos < 0.5) {
    tc$manual_cluster[i] <- "Cluster 4"
  } else {
    tc$manual_cluster[i] <- "Cluster 6"
  }
}
tc <- tc[!is.na(tc$manual_cluster), ]

png(file.path(OUT_DIR, "fig3B_naive.png"),
    width = 6, height = 4, units = "in", res = 600, bg = "white")

plot_clusters_mean_scaled(
  coord_x = tc$x,
  coord_y = tc$y,
  gene_clusters = tc$manual_cluster,
  xp = tc$gene,
  cutoff = 20,
  title = "Final Naive Clusters",
  ylim = c(0, 3), n_col = 3)
dev.off()

# --- tumor-bearing ---
traj_clusters <- readRDS(file.path(DATA_DIR, "cancer_traj_clusters.rds"))

#take the automatic clusters and manually clean them (seperate di and dfi into 2 clusters)
tc <- traj_clusters
cluster_labels <- as.character(tc$cluster)

cluster1  <- c("fi", "i", "fdi", "ifi", "ifdi")
cluster5  <- c("d", "df")
cluster4  <- c("id")
cluster6  <- c("idi")
cluster7  <- c("did")
cluster23 <- c("di", "dfi")

tc$manual_cluster <- NA

# combine into broad clusters first
tc$manual_cluster[tc$cluster %in% cluster1]  <- "Cluster 1"
tc$manual_cluster[tc$cluster %in% cluster5]  <- "Cluster 5"
tc$manual_cluster[tc$cluster %in% cluster4]  <- "Cluster 4"
tc$manual_cluster[tc$cluster %in% cluster6]  <- "Cluster 6"
tc$manual_cluster[tc$cluster %in% cluster7]  <- "Cluster 7"
tc$manual_cluster[tc$cluster %in% cluster23] <- "cluster23_di_dfi"

idx23 <- which(tc$manual_cluster == "cluster23_di_dfi")

for (i in idx23) {
  x <- as.numeric(tc$x[[i]])
  y <- as.numeric(tc$y[[i]])

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  min_pos <- x[which.min(y)]

  if (min_pos < 0.5) {
    tc$manual_cluster[i] <- "Cluster 2"
  } else {
    tc$manual_cluster[i] <- "Cluster 3"
  }
}
tc <- tc[!is.na(tc$manual_cluster), ]

png(file.path(OUT_DIR, "fig3B_cancer.png"),
    width = 6, height = 4, units = "in", res = 600, bg = "white")

plot_clusters_mean_scaled(
  coord_x = tc$x,
  coord_y = tc$y,
  gene_clusters = tc$manual_cluster,
  xp = tc$gene,
  cutoff = 20,
  title = "Final Cancer Clusters",
  ylim = c(0, 3), n_col = 4)
dev.off()


# figure 3C
tfs_nnb <- read.csv(file.path(DATA_DIR, "tfs_nnb.csv"), row.names = 1)
tfs_ncb <- read.csv(file.path(DATA_DIR, "tfs_ncb.csv"), row.names = 1)

# -----------------------------------
# 1. Clean each TF list
# -----------------------------------

prep_tfs = function(x, list_name) {

  x = as.data.frame(x, stringsAsFactors = FALSE)

  df = data.frame(
    tf = as.character(x[[1]]),
    q_value = as.numeric(x[[4]]),
    stringsAsFactors = FALSE
  )

  df = df %>%
    dplyr::filter(!is.na(tf), !is.na(q_value)) %>%
    dplyr::group_by(tf) %>%
    dplyr::slice_min(order_by = q_value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(q_value) %>%
    dplyr::mutate(
      rank = dplyr::row_number(),
      neglogq = -log(pmax(q_value, .Machine$double.xmin)),
      list = list_name
    )

  return(df)
}

nnb = prep_tfs(tfs_nnb, "NNB")
ncb = prep_tfs(tfs_ncb, "NCB")

# -----------------------------------
# 2. Overlap across full lists
# -----------------------------------

overlap_tfs = intersect(nnb$tf, ncb$tf)

nnb = nnb %>%
  dplyr::mutate(overlap = tf %in% overlap_tfs)

ncb = ncb %>%
  dplyr::mutate(overlap = tf %in% overlap_tfs)

# -----------------------------------
# 3. Keep top 20 and arrange horizontally
# -----------------------------------

top_n = 20

# more horizontal spacing between circles
x_spacing = 2.35
x_vals = seq(1, by = x_spacing, length.out = top_n)

# two rows
y_top = 1.22
y_bottom = 0.92

nnb_top = nnb %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::mutate(
    x = x_vals,
    y = y_top
  )

ncb_top = ncb %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::mutate(
    x = x_vals,
    y = y_bottom
  )

plot_df = dplyr::bind_rows(nnb_top, ncb_top)

# -----------------------------------
# 4. Colors
# -----------------------------------

col_overlap = "#D55E00"
col_nonoverlap = "#A8A8A8"
col_guides = "#D9D9D9"

# -----------------------------------
# 5. Plot
# -----------------------------------

x_max = max(x_vals)

p = ggplot(plot_df, aes(x = x, y = y)) +

  # horizontal guide lines
  geom_segment(
    aes(x = min(x_vals), xend = max(x_vals), y = y_top, yend = y_top),
    color = col_guides,
    linewidth = 0.9
  ) +
  geom_segment(
    aes(x = min(x_vals), xend = max(x_vals), y = y_bottom, yend = y_bottom),
    color = col_guides,
    linewidth = 0.9
  ) +

  # circles
  geom_point(
    aes(size = neglogq, color = overlap),
    alpha = 0.9
  ) +

  # top-row labels: close to circles
  geom_text(
    data = nnb_top %>% dplyr::filter(!overlap),
    aes(x = x, y = y_top + 0.1, label = tf),
    angle = 50,
    vjust = 0.9,
    hjust = 0,
    size = 3.7,
    color = "black"
  ) +
  geom_text(
    data = nnb_top %>% dplyr::filter(overlap),
    aes(x = x, y = y_top + 0.1, label = tf),
    angle = 50,
    vjust = 0.9,
    hjust = 0,
    size = 3.7,
    color = col_overlap
  ) +

  # bottom-row labels: also close to circles and above them
  geom_text(
    data = ncb_top %>% dplyr::filter(!overlap),
    aes(x = x, y = y_bottom + 0.1, label = tf),
    angle = 50,
    vjust = 0.9,
    hjust = 0,
    size = 3.7,
    color = "black"
  ) +
  geom_text(
    data = ncb_top %>% dplyr::filter(overlap),
    aes(x = x, y = y_bottom + 0.1, label = tf),
    angle = 50,
    vjust = 0.9,
    hjust = 0,
    size = 3.7,
    color = col_overlap
  ) +

  scale_x_continuous(
    breaks = NULL,
    labels = NULL,
    limits = c(min(x_vals) - 1.0, x_max + 16)
  ) +

  scale_y_continuous(
    breaks = NULL,
    labels = NULL,
    limits = c(0.78, 1.48)
  ) +

  scale_color_manual(
    values = c(`TRUE` = col_overlap, `FALSE` = col_nonoverlap),
    labels = c(`TRUE` = "Significant in Both", `FALSE` = "Significant in One"),
    name = NULL
  ) +

  scale_size_continuous(
    range = c(4.2, 8.6),
    breaks = c(5, 10, 15),
    name = expression(-log(q-value))
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  coord_cartesian(clip = "off") +

  theme_classic(base_size = 16) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),

    # smaller legend
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = c(0.78, 0.58),
    legend.justification = c(0, 0.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.7, "cm"),

    plot.margin = margin(15, 120, 15, 15)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    size = guide_legend(order = 2, override.aes = list(color = "black"))
  )


ggsave(
  filename = file.path(OUT_DIR, "fig3C.pdf"),
  plot = p,
  width = 12,
  height = 2.2,
  units = "in"
)


# figure 3D
plot_network_both = function(tf_links1, tf_links2){
  topology1 = data.frame(as.matrix(tf_links1), stringsAsFactors = FALSE)
  topology2 = data.frame(as.matrix(tf_links2), stringsAsFactors = FALSE)
  node_list <- unique(c(topology1[,1], topology1[,2], topology2[,1], topology2[,2]))
  nodes <- data.frame(id = node_list,
                      label = node_list,
                      font.size = min(150/nchar(node_list), 30),
                      shape='circle',
                      shapeProperties = list(useBorderWithImage = TRUE),
                      widthConstraint = list(minimum = 75, maximum = 75),
                      heightConstraint = list(minimum = 75, maximum = 75))

  only1 = "rgba(100, 0, 0, 0.4)"
  only2 = "rgba(0, 100,0, 0.4)"
  both = "rgba(0, 0, 255, 1)"
  edges <- data.frame(from = c(topology1[,1]),
                      to = c(topology1[,2]),
                      arrows.to.type = "arrow",
                      width = 2, color = only1)
  for (i in 1:length(topology2[,1])){
    found = FALSE
    for (j in 1:length(edges[,1])){
      if (topology2[,1][i] == edges$from[j] && topology2[,2][i] == edges$to[j]){
        edges$color[j] = both
        edges$width = 4
        found = TRUE
      }
    }
    if (!found){
      new_row = data.frame (from = topology2[,1][i], to = topology2[,2][i], arrows.to.type = "arrow", width = 3, color = only2)
      edges = rbind(edges, new_row)
    }
  }
  #edges$color <- replicate(nrow(edges), list(color = "green", opacity = 0.1), simplify = FALSE)
  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = TRUE)
}

naive_three_networks_8   <- read.csv(file.path(DATA_DIR, "naive_initial_network.csv"),  row.names = 1)
cancer_three_networks_10 <- read.csv(file.path(DATA_DIR, "cancer_initial_network.csv"), row.names = 1)

plot_network_both(cancer_three_networks_10, naive_three_networks_8)
