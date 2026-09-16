library(princurve)
library(ggplot2)
library(patchwork)
library(dplyr)

DATA_DIR <- "../data"
OUT_DIR  <- "../plots"


# shared setup: baseline clouds -> principal curve -> pseudotime per state
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


# figure 6A
top_naive  <- c("Egr1_Cebpb", "Ets1")
top_cancer <- c("Cebpb_Jun",  "Ets1")

# pseudotime for a perturbed cloud, using the matching baseline curve
pt_of <- function(pert_name, base, pca, curve, base_proj, dir_, prefix) {
  f    <- file.path(dir_, paste0(prefix, pert_name, ".csv"))
  pert <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE)); storage.mode(pert) <- "numeric"
  pert <- pert[rownames(base), , drop = FALSE]
  pc   <- predict(pca, newdata = t(pert))[, 1:2]
  lam  <- project_to_curve(pc, s = curve, stretch = 0)$lambda
  pt   <- 1 - (lam - min(base_proj$lambda)) / (max(base_proj$lambda) - min(base_proj$lambda))
  pmin(pmax(pt, 0), 1)
}

# one panel: baseline vs one perturbation, NO legend, label in the title only
panel <- function(base_pt, pert_pt, pert_name, condition) {
  lab <- gsub("_", " + ", pert_name)
  df  <- rbind(data.frame(pseudotime = base_pt, group = "Baseline"),
               data.frame(pseudotime = pert_pt, group = "Perturbed"))
  df$group <- factor(df$group, levels = c("Baseline", "Perturbed"))
  ggplot(df, aes(pseudotime, fill = group, color = group)) +
    geom_density(alpha = 0.45, linewidth = 1, adjust = 1) +
    scale_fill_manual(values = c(Baseline = "grey60", Perturbed = "orange"),
                      name = NULL) +
    scale_color_manual(values = c(Baseline = "grey60", Perturbed = "orange"),
                       name = NULL) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = paste0(condition, ": ", lab), x = "Principal Curve Coordinate",
         y = "State Density") +
    theme_classic(base_size = 13) +
    theme(plot.title  = element_text(size = 13, face = "bold"),
          plot.margin = margin(3, 3, 3, 3))
}
np1 <- panel(naive_base_pt,  pt_of(top_naive[1],  naive_base,  naive_pca,  naive_curve,  naive_base_proj,  DATA_DIR,  "naive_perturbed_"),  top_naive[1],  "Naive")
np2 <- panel(naive_base_pt,  pt_of(top_naive[2],  naive_base,  naive_pca,  naive_curve,  naive_base_proj,  DATA_DIR,  "naive_perturbed_"),  top_naive[2],  "Naive")
cp1 <- panel(cancer_base_pt, pt_of(top_cancer[1], cancer_base, cancer_pca, cancer_curve, cancer_base_proj, DATA_DIR, "cancer_perturbed_"), top_cancer[1], "Tumor-Bearing")
cp2 <- panel(cancer_base_pt, pt_of(top_cancer[2], cancer_base, cancer_pca, cancer_curve, cancer_base_proj, DATA_DIR, "cancer_perturbed_"), top_cancer[2], "Tumor-Bearing")

fig <- (np1 | cp1) / (np2 | cp2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.margin = margin(2, 2, 2, 2))

ggsave(file.path(OUT_DIR, "fig6A.pdf"), fig, width = 6.2, height = 4.4)


# figure 6B
pretty_pair <- function(p) {
  g <- strsplit(p, "_")[[1]]
  if (length(g) == 2 && g[1] == g[2]) g <- g[1]   # _gene1_gene1 -> gene1
  paste(g, collapse = " + ")
}

naive_sweep  <- read.csv(file.path(DATA_DIR, "fig6B_naive_sweep.csv"))
cancer_sweep <- read.csv(file.path(DATA_DIR, "fig6B_cancer_sweep.csv"))

# anchor each pair to its own per = 0 cloud
naive_sweep <- naive_sweep %>%
  group_by(pair) %>%
  mutate(mean_diff = perturbed_mean - perturbed_mean[per == 0]) %>%
  ungroup() %>%
  arrange(pair, per)

print(naive_sweep)

naive_sweep$label <- vapply(naive_sweep$pair, pretty_pair, character(1))
p_naive <- ggplot(naive_sweep, aes(x = per, y = mean_diff, color = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  labs(
    title = "Naive",
    x = "Perturbation Strength",
    y = "Principal Curve Shift"
  ) +
  theme_classic()


# anchor each pair to its own per = 0 cloud
cancer_sweep <- cancer_sweep %>%
  group_by(pair) %>%
  mutate(mean_diff = perturbed_mean - perturbed_mean[per == 0]) %>%
  ungroup() %>%
  arrange(pair, per)

print(cancer_sweep)

cancer_sweep$label <- vapply(cancer_sweep$pair, pretty_pair, character(1))
p_cancer <- ggplot(cancer_sweep, aes(x = per, y = mean_diff, color = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  labs(
    title = "Tumor-Bearing",
    x = "Perturbation Strength",
    y = "Principal Curve Shift"
  ) +
  theme_classic()



ggsave(file.path(OUT_DIR, "fig6B_cancer.pdf"), p_cancer, width = 4, height = 3)
ggsave(file.path(OUT_DIR, "fig6B_naive.pdf"), p_naive, width = 4, height = 3)
