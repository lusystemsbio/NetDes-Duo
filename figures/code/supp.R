library(Seurat)
library(princurve)
library(ggplot2)
library(readr)
library(visNetwork)

DATA_DIR <- "../data"
OUT_DIR  <- "../plots"


# figure S1
neutrophil <- readRDS(file.path(DATA_DIR, "combined_neutrophil_pca.rds"))
# level order pins the palette
neutrophil$condition <- factor(neutrophil$condition, levels = c("Tumor-Bearing", "Naive"))

# PCA plot by MULTI_ID (6 colors)
sdev <- neutrophil[["pca"]]@stdev
pct.var <- (sdev^2 / sum(sdev^2)) * 100
pc1_pct = signif(pct.var[[1]], digits = 3)
pc2_pct = signif(pct.var[[2]], digits = 3)

p <- PCAPlot(neutrophil, group.by = "condition", pt.size = 0.25) +
  labs(title = NULL, x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"), color = NULL) +
  theme_classic(base_size = 7) +
  theme(axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  coord_cartesian(xlim = c(-15, 20), ylim = c(-10, 15)) +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))

ggsave(
  filename = file.path(OUT_DIR, "figS1.pdf"),
  plot = p, width = 150, height = 100, units = "mm")


# figures S2 and S3
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

# figure S2
traj_clusters <- readRDS(file.path(DATA_DIR, "cancer_traj_clusters.rds"))

png(
  file.path(OUT_DIR, "figS2.png"),
  width = 14,
  height = 8,
  units = "in",
  res = 600,
  bg = "white"
)

plot_clusters_mean_scaled(
  coord_x = traj_clusters$x,
  coord_y = traj_clusters$y,
  gene_clusters = traj_clusters$cluster,
  xp = traj_clusters$gene,
  cutoff = 20,
  title = "Initial Tumor-Bearing Clusters",
  ylim = c(0, 5)
)

dev.off()

# figure S3
traj_clusters <- readRDS(file.path(DATA_DIR, "naive_traj_clusters.rds"))

png(
  file.path(OUT_DIR, "figS3.png"),
  width = 14,
  height = 8,
  units = "in",
  res = 600,
  bg = "white"
)

plot_clusters_mean_scaled(
  coord_x = traj_clusters$x,
  coord_y = traj_clusters$y,
  gene_clusters = traj_clusters$cluster,
  xp = traj_clusters$gene,
  cutoff = 20,
  title = "Initial Naive Clusters",
  ylim = c(0, 5)
)

dev.off()


# figure S10
mats <- list(); labs <- character()
for (sd in c(125, 177, 629, 753, 1419, 1434, 1477, 1594, 1612, 2120,
             2569, 3219, 3646, 4119, 4250, 4510, 5042, 5236, 6072, 6534)) for (i in 1:2) {
  df <- read_csv(file.path(DATA_DIR, sprintf("figS10_synth_%d_%d.csv", sd, i)), show_col_types = FALSE)
  m <- as.matrix(df[, -1]); storage.mode(m) <- "numeric"; rownames(m) <- as.character(df[[1]])
  mats[[length(mats) + 1]] <- m; labs <- c(labs, sprintf("Seed %d - Net %d", sd, i))
}
ylim <- range(vapply(mats, range, numeric(2)))

pdf(file.path(OUT_DIR, "figS10.pdf"), width = 20, height = 12)
par(mfrow = c(4, 5), mar = c(2.4, 2.8, 2.0, 0.4), oma = c(0.2, 0.2, 0.2, 0.2),
    cex.axis = 1, cex.main = 1.45)
for (k in seq_along(mats))
  matplot(seq_len(ncol(mats[[k]])), t(mats[[k]]), xlab = "", ylab = "", type = "l",
          lty = 1, lwd = 2, col = seq_len(nrow(mats[[k]])), ylim = ylim, main = labs[k])
dev.off()


# figure S5
plot_network = function(tf_links = tf_links){
  topology = data.frame(as.matrix(tf_links), stringsAsFactors = FALSE)

  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list, label = node_list,
                      font.size = min(150/nchar(node_list), 30),
                      shape='circle',
                      shapeProperties = list(useBorderWithImage = TRUE),
                      widthConstraint = list(minimum = 75, maximum = 75),
                      heightConstraint = list(minimum = 75, maximum = 75))
  edge_col <- data.frame(c(1, 2), c("blue", "darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1, 2), c("arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")

  edges <- data.frame(from = c(topology[,1]), to = c(topology[,2]),
                      arrows.to.type = arrow_type$color[c(as.numeric(topology[,3]))],
                      width = 3,
                      color = edge_col$color[c(as.numeric(topology[,3]))])

  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
}

naive_signed <- read_csv(file.path(DATA_DIR, "naive_final_network_signed.csv"), show_col_types = FALSE)[, c("Source","Target", "Interaction")]
cancer_signed <- read_csv(file.path(DATA_DIR, "cancer_final_network_signed.csv"), show_col_types = FALSE)[, c("Source","Target", "Interaction")]

plot_network(naive_signed)
plot_network(cancer_signed)


# figure S7
# setup: baseline clouds -> principal curve -> pseudotime per state
naive_base <- read.csv(
  file.path(DATA_DIR, "naive_baseline_cloud.csv"),
  row.names = 1, check.names = FALSE)

naive_base <- as.matrix(naive_base)
storage.mode(naive_base) <- "numeric"

naive_pca <- prcomp(t(naive_base), center = TRUE, scale. = FALSE)

naive_base_pcs <- data.frame(
  PC1 = naive_pca$x[, 1],
  PC2 = naive_pca$x[, 2]
)

naive_curve_fit <- principal_curve(
  as.matrix(naive_base_pcs[, c("PC1", "PC2")]),
  smoother = "smooth_spline",
  maxit = 50,
  stretch = 0
)

naive_curve <- naive_curve_fit$s[naive_curve_fit$ord, , drop = FALSE]

naive_base_proj <- project_to_curve(
  x = as.matrix(naive_base_pcs[, c("PC1", "PC2")]),
  s = naive_curve,
  stretch = 0
)

naive_base_pt <- naive_base_proj$lambda
naive_base_pt <- (naive_base_pt - min(naive_base_pt)) / (max(naive_base_pt) - min(naive_base_pt))
naive_base_pt <- 1 - naive_base_pt
naive_base_pcs$pseudotime <- naive_base_pt

cancer_base <- read.csv(file.path(DATA_DIR, "cancer_baseline_cloud.csv"),
                        row.names = 1, check.names = FALSE)

cancer_base <- as.matrix(cancer_base)
storage.mode(cancer_base) <- "numeric"

cancer_pca <- prcomp(t(cancer_base), center = TRUE, scale. = FALSE)

cancer_base_pcs <- data.frame(
  PC1 = cancer_pca$x[, 1],
  PC2 = cancer_pca$x[, 2]
)

cancer_curve_fit <- principal_curve(
  as.matrix(cancer_base_pcs[, c("PC1", "PC2")]),
  smoother = "smooth_spline",
  maxit = 50,
  stretch = 0
)

cancer_curve <- cancer_curve_fit$s[cancer_curve_fit$ord, , drop = FALSE]

cancer_base_proj <- project_to_curve(
  x = as.matrix(cancer_base_pcs[, c("PC1", "PC2")]),
  s = cancer_curve,
  stretch = 0
)

cancer_base_pt <- cancer_base_proj$lambda
cancer_base_pt <- (cancer_base_pt - min(cancer_base_pt)) /
  (max(cancer_base_pt) - min(cancer_base_pt))
cancer_base_pt <- 1 - cancer_base_pt
cancer_base_pcs$pseudotime <- cancer_base_pt

nnb <- readRDS(file.path(DATA_DIR, "nnb_expr.rds"))
ncb <- readRDS(file.path(DATA_DIR, "ncb_expr.rds"))

organ_levels <- c("Lung", "Blood", "BM")
organ_cols   <- c(Blood = "#F8766D", BM = "#00BA38", Lung = "#619CFF")

big <- theme(axis.title = element_text(size = 22),
             axis.text  = element_text(size = 18),
             strip.text = element_blank(),
             legend.title = element_text(size = 20),
             legend.text  = element_text(size = 18),
             strip.background = element_blank())


# Project real cells into the simulated PCA by z-scoring each gene, then predict().
project_real <- function(cloud, pca, obj) {
  genes  <- rownames(cloud)
  expr   <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
  shared <- intersect(genes, rownames(expr))

  aligned <- matrix(0, nrow = length(genes), ncol = ncol(expr),
                    dimnames = list(genes, colnames(expr)))
  for (g in shared) {
    x <- expr[g, ]
    if (sd(x) > 0) aligned[g, ] <- (x - mean(x)) / sd(x)
  }
  list(scores = predict(pca, t(aligned))[, 1:2],
       organ = factor(as.character(obj$MULTI_ID), levels = organ_levels))
}

naive_proj  <- project_real(naive_base,  naive_pca,  nnb)
cancer_proj <- project_real(cancer_base, cancer_pca, ncb)


# Pseudotime along the simulated curve; stretch = 2 lets the ends extrapolate.
real_pseudotime <- function(scores, curve, base_proj) {
  pr <- project_to_curve(as.matrix(scores), s = curve, stretch = 2)
  1 - (pr$lambda - min(base_proj$lambda)) / (max(base_proj$lambda) - min(base_proj$lambda))
}

naive_real_pt  <- real_pseudotime(naive_proj$scores,  naive_curve,  naive_base_proj)
cancer_real_pt <- real_pseudotime(cancer_proj$scores, cancer_curve, cancer_base_proj)


# Plot 2 -- bar: tissue makeup across ten pseudotime bins.
bin_table <- function(pt, organ, n = 10) {
  bins <- cut(pt, breaks = seq(min(pt), max(pt), length.out = n + 1),
              labels = FALSE, include.lowest = TRUE)
  out <- as.data.frame(prop.table(table(bin = factor(bins, 1:n), organ = organ), margin = 1))
  out$bin <- as.integer(as.character(out$bin))
  out
}
bar_df <- rbind(
  cbind(bin_table(naive_real_pt,  naive_proj$organ),  condition = "Naive"),
  cbind(bin_table(cancer_real_pt, cancer_proj$organ), condition = "Tumor-Bearing")
)
bar_df$condition <- factor(bar_df$condition, levels = c("Naive", "Tumor-Bearing"))

p_bar <- ggplot(bar_df, aes(bin, Freq, fill = organ)) +
  geom_col(width = 0.8) +
  facet_wrap(~ condition) +
  scale_fill_manual(values = organ_cols) +
  scale_x_continuous(breaks = c(0.5, 3, 5.5, 8, 10.5),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  labs(x = "Principal Curve Coordinate", y = "Tissue Fraction", fill = "Organ") +
  theme_classic() +
  big +
  theme(panel.spacing.x = unit(3, "lines"))


# Plot 3 -- PCA overlay (one per condition); PC1 flipped, axes clipped to 1-99%.
pca_plot <- function(proj, pcs, curve) {
  real <- data.frame(PC1 = -proj$scores[, 1], PC2 = proj$scores[, 2], organ = proj$organ)
  sim  <- data.frame(PC1 = -pcs$PC1, PC2 = pcs$PC2)
  line <- data.frame(PC1 = -curve[, "PC1"], PC2 = curve[, "PC2"])
  ggplot() +
    geom_point(data = real, aes(PC1, PC2, color = organ), size = 0.5, alpha = 0.45) +
    geom_point(data = sim, aes(PC1, PC2), color = "grey25", size = 0.4, alpha = 0.6) +
    geom_path(data = line, aes(PC1, PC2), color = "black", linewidth = 1) +
    coord_cartesian(xlim = quantile(real$PC1, c(0.01, 0.99)),
                    ylim = quantile(real$PC2, c(0.01, 0.99))) +
    scale_color_manual(values = organ_cols) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "PC1", y = "PC2", color = "Organ") +
    theme_classic() + big
}
p_pca_naive  <- pca_plot(naive_proj,  naive_base_pcs,  naive_curve)
p_pca_cancer <- pca_plot(cancer_proj, cancer_base_pcs, cancer_curve)


ggsave(file.path(OUT_DIR, "figS7A.pdf"),        p_bar,        width = 8, height = 5)
ggsave(file.path(OUT_DIR, "figS7B_naive.pdf"),  p_pca_naive,  width = 5, height = 4)
ggsave(file.path(OUT_DIR, "figS7B_cancer.pdf"), p_pca_cancer, width = 5, height = 4)

# figure S11
x          <- c(27, 30, 33, 36, 39, 42)
x_chr      <- as.character(x)

cols <- c(
  `NetDes-Duo` = "#B07AA1",
  NetDes   = "#F28E2B",
  ppcor         = "#E15759",
  SCODE         = "#76B7B2",
  GENIE3 = "#59A14F",
  DeepSEM = "#FF9DA7"
)

dfs <- list(
  `NetDes-Duo` = read.csv(file.path(DATA_DIR, "figS11_jaccard_NetDes_Duo.csv")),
  NetDes       = read.csv(file.path(DATA_DIR, "figS11_jaccard_NetDes.csv")),
  GENIE3       = read.csv(file.path(DATA_DIR, "figS11_jaccard_GENIE3.csv")),
  ppcor        = read.csv(file.path(DATA_DIR, "figS11_jaccard_ppcor.csv")),
  SCODE        = read.csv(file.path(DATA_DIR, "figS11_jaccard_SCODE.csv")),
  DeepSEM      = read.csv(file.path(DATA_DIR, "figS11_jaccard_DeepSEM.csv"))
)

methods <- names(dfs)
mu   <- matrix(NA_real_, nrow = length(methods), ncol = length(x), dimnames = list(methods, x_chr))
sdv  <- mu
nmat <- mu

for (m in methods) {
  d <- dfs[[m]]
  for (j in seq_along(x)) {
    v <- d$jaccard[d$total_edges == x[j]]
    mu[m, j]   <- mean(v, na.rm = TRUE)
    sdv[m, j]  <- sd(v, na.rm = TRUE)
    nmat[m, j] <- sum(!is.na(v))
  }
}
sem <- sdv / sqrt(nmat)

ymin <- 0
ymax <- 1

pdf(file.path(OUT_DIR, "figS11.pdf"), width = 8, height = 6)

par(mfrow = c(1, 1),
    mar = c(5.2, 5.4, 3.2, 1.0),
    mgp = c(3.2, 0.9, 0),
    cex.lab = 1.45,
    cex.axis = 1.25,
    cex.main = 1.55)

bp <- barplot(mu,
              beside = TRUE,
              col = cols[rownames(mu)],
              border = "black",
              ylim = c(ymin, ymax),
              names.arg = x_chr,
              xlab = "Initial Network Size",
              ylab = "Average Jaccard Index",
              main = "Network Agreement",
              xpd = FALSE)

legend("topright",
       legend = rownames(mu),
       fill   = cols[rownames(mu)],
       bty    = "n",
       cex    = 1.05)

arrows(x0 = bp, y0 = mu - sem,
       x1 = bp, y1 = mu + sem,
       angle = 90, code = 3, length = 0.04, lwd = 1.2)

dev.off()


# figure S14
cols <- c(
  `NetDes-Duo` = "#B07AA1",
  NetDes  = "#F28E2B",
  GENIE3  = "#59A14F",
  ppcor   = "#E15759",
  SCODE   = "#76B7B2",
  DeepSEM = "#FF9DA7"
)

# naive vs tumor-bearing overlap per method, cut to the 12 core TFs
d <- read.csv(file.path(DATA_DIR, "figS14_overlap_jaccard.csv"))
d <- d[match(names(cols), d$method), ]

pdf(file.path(OUT_DIR, "figS14.pdf"), width = 8, height = 6)

par(mfrow = c(1, 1),
    mar = c(5.2, 5.4, 3.2, 1.0),
    mgp = c(3.2, 0.9, 0),
    cex.lab = 1.45,
    cex.axis = 1.25,
    cex.main = 1.55)

barplot(d$jaccard,
        names.arg = d$method,
        col = cols[d$method],
        border = "black",
        ylim = c(0, 1),
        ylab = "Jaccard Index",
        main = "Naive vs Tumor-Bearing Network Overlap",
        cex.names = 1.05,
        xpd = FALSE)

dev.off()
