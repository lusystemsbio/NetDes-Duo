library(princurve)
library(SeuratObject)
library(ggplot2)
library(ggridges)

DATA_DIR <- "../data"
OUT_DIR  <- "../plots"


# shared setup: baseline clouds -> principal curve -> pseudotime per state
big <- theme(plot.title  = element_text(size = 16, face = "bold"),
             axis.title  = element_text(size = 15),
             axis.text   = element_text(size = 12),
             strip.text  = element_text(size = 13, face = "bold"),
             legend.text = element_text(size = 12),
             strip.background = element_blank())

naive_base <- read.csv(
  file.path(DATA_DIR, "naive_baseline_cloud.csv"),
  row.names = 1, check.names = FALSE)

naive_base <- as.matrix(naive_base)
storage.mode(naive_base) <- "numeric"

# PCA on baseline only
naive_pca <- prcomp(t(naive_base), center = TRUE, scale. = FALSE)

naive_base_pcs <- data.frame(
  PC1 = naive_pca$x[, 1],
  PC2 = naive_pca$x[, 2]
)

# principal curve on baseline PCA only
naive_curve_fit <- principal_curve(
  as.matrix(naive_base_pcs[, c("PC1", "PC2")]),
  smoother = "smooth_spline",
  maxit = 50,
  stretch = 0
)

naive_curve <- naive_curve_fit$s[naive_curve_fit$ord, , drop = FALSE]

# baseline pseudotime
naive_base_proj <- project_to_curve(
  x = as.matrix(naive_base_pcs[, c("PC1", "PC2")]),
  s = naive_curve,
  stretch = 0
)

naive_base_pt <- naive_base_proj$lambda
naive_base_pt <- (naive_base_pt - min(naive_base_pt)) / (max(naive_base_pt) - min(naive_base_pt))
#fix sign
naive_base_pt <- 1 - naive_base_pt
naive_base_pcs$pseudotime <- naive_base_pt
naive_base_mean <- mean(naive_base_pt)

#cancer case
cancer_base <- read.csv(file.path(DATA_DIR, "cancer_baseline_cloud.csv"),
                        row.names = 1, check.names = FALSE)

cancer_base <- as.matrix(cancer_base)
storage.mode(cancer_base) <- "numeric"

# PCA on baseline only
cancer_pca <- prcomp(t(cancer_base), center = TRUE, scale. = FALSE)

cancer_base_pcs <- data.frame(
  PC1 = cancer_pca$x[, 1],
  PC2 = cancer_pca$x[, 2]
)

# principal curve on baseline PCA only
cancer_curve_fit <- principal_curve(
  as.matrix(cancer_base_pcs[, c("PC1", "PC2")]),
  smoother = "smooth_spline",
  maxit = 50,
  stretch = 0
)

cancer_curve <- cancer_curve_fit$s[cancer_curve_fit$ord, , drop = FALSE]

# baseline pseudotime
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
cancer_base_mean <- mean(cancer_base_pt)


# figure 5B
# --- helper: compute quasi-potential curve from baseline pseudotime values ---
quasi_potential <- function(pt, n_grid = 512, bw_adjust = 1) {
  # pt: vector of baseline states' pseudotime (lambda), already in [0,1]
  d <- density(pt, adjust = bw_adjust, n = n_grid, from = 0, to = 1)
  P <- d$y
  P[P <= 0] <- NA                     # guard against log(0); KDE rarely hits 0 but be safe
  U <- -log(P)
  U <- U - min(U, na.rm = TRUE)       # zero the global minimum for cross-condition comparability
  data.frame(pseudotime = d$x, U = U)
}

# --- build curves for both conditions ---

naive_U  <- quasi_potential(naive_base_pt,  bw_adjust = 1)
cancer_U <- quasi_potential(cancer_base_pt, bw_adjust = 1)

naive_U$condition  <- "Naive"
cancer_U$condition <- "Tumor-Bearing"
both <- rbind(naive_U, cancer_U)

p_potential <- ggplot(both, aes(pseudotime, U)) +
  geom_area(fill = "#CFE0EE") +
  geom_line(linewidth = 1.2, color = "#1F3B57") +
  facet_wrap(~ condition) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Quasi-Potential Landscape Along Maturation Pseudotime",
       x = "Principal Curve Coordinate",
       y = "Quasi-Potential") +
  theme_classic() + big +
  theme(panel.spacing.x = unit(2, "lines"))

ggsave(file.path(OUT_DIR, "fig5B.pdf"), p_potential, width = 5.0, height = 3.0)


# figure 5C
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


# Plot 1 -- ridge: each tissue's pseudotime distribution, with the simulated one for reference.
ridge_df <- rbind(
  data.frame(pseudotime = naive_base_pt,  group = "Simulated",                     condition = "Naive"),
  data.frame(pseudotime = naive_real_pt,  group = as.character(naive_proj$organ),  condition = "Naive"),
  data.frame(pseudotime = cancer_base_pt, group = "Simulated",                     condition = "Tumor-Bearing"),
  data.frame(pseudotime = cancer_real_pt, group = as.character(cancer_proj$organ), condition = "Tumor-Bearing")
)
ridge_df$group     <- factor(ridge_df$group, levels = c("Simulated", organ_levels))
ridge_df$condition <- factor(ridge_df$condition, levels = c("Naive", "Tumor-Bearing"))

p_ridge <- ggplot(ridge_df, aes(pseudotime, group, fill = group)) +
  geom_density_ridges(alpha = 0.8, scale = 1.35, rel_min_height = 0.012,
                      quantile_lines = TRUE, quantiles = 2, color = "grey20", linewidth = 0.5) +
  facet_wrap(~ condition) +
  scale_fill_manual(values = c(Baseline = "grey75", organ_cols)) +
  labs(x = NULL, y = NULL) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(legend.position = "none", panel.spacing.x = unit(3, "lines")) +
  big

ggsave(file.path(OUT_DIR, "fig5C.pdf"),  p_ridge,      width = 8.0, height = 4.0)
