library(readr)
library(dplyr)
library(visNetwork)

DATA_DIR <- "../data"
OUT_DIR  <- "../plots"


# figure 2A
noise_levels = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)

plot_noise_comparison <- function() {

  base <- DATA_DIR

  data_list <- vector("list", length(noise_levels))
  names(data_list) <- as.character(noise_levels)

  # read all data first
  for (i in c(1,2)){
    f <- paste0(base, "/fig2A_data_", i, ".csv")

    df <- read_csv(f, show_col_types = FALSE)
    genes <- as.character(df[[1]])
    mat <- as.matrix(df[, -1])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- genes

    data_list[[i]] <- mat
  }
  # same y-axis for all panels
  y_min <- min(vapply(data_list, min, numeric(1), na.rm = TRUE))
  y_max <- max(vapply(data_list, max, numeric(1), na.rm = TRUE))

  # colors for genes
  gene_names <- rownames(data_list[[1]])
  gene_cols <- seq_along(gene_names)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(
    mfrow = c(1,2),
    mar = c(2.4, 2.8, 2.0, 0.4),
    oma = c(0.2, 0.2, 0.2, 0.2),
    cex.axis = 1,
    cex.main = 1.45)
  for (i in c(1,2)) {
    mat <- data_list[[i]]

    # x-axis is time index
    x <- seq_len(ncol(mat))
    matplot(x, t(mat),
            xlab = "",
            ylab = "",
            type = "l", lty = 1, lwd = 2,
            col = gene_cols,
            ylim = c(y_min, y_max),
            main = paste0("Network ", i))

  }
}

pdf(file.path(OUT_DIR, "fig2A.pdf"), width = 8, height = 3)
plot_noise_comparison()
dev.off()


# figure 2B
plot_network_both_signed = function(tf_links1, tf_links2) {
  topology1 = data.frame(as.matrix(tf_links1), stringsAsFactors = FALSE)
  topology2 = data.frame(as.matrix(tf_links2), stringsAsFactors = FALSE)
  topology1[,3] = as.integer(topology1[,3])
  topology2[,3] = as.integer(topology2[,3])

  node_list = unique(c(topology1[,1], topology1[,2], topology2[,1], topology2[,2]))
  nodes <- data.frame(id = node_list,
                      label = node_list,
                      font.size = min(150/nchar(node_list), 30),
                      shape='circle',
                      shapeProperties = list(useBorderWithImage = TRUE),
                      widthConstraint = list(minimum = 75, maximum = 75),
                      heightConstraint = list(minimum = 75, maximum = 75))

  only1_col = "rgba(100, 0, 0, 0.5)"
  only2_col = "rgba(0, 100, 0, 0.5)"
  both_col  = "rgba(0, 0, 255, 1)"

  arrow_type_from = function(type_vec) ifelse(type_vec == 2, "circle", "arrow")
  edges = data.frame(
    from = topology1[,1],
    to   = topology1[,2],
    reg_type = topology1[,3],  # 1=activation, 2=repression
    arrows.to.type = arrow_type_from(topology1[,3]),
    width  = 4,
    color  = only1_col,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(topology2))) {
    f = topology2[i, 1]
    t = topology2[i, 2]
    r = topology2[i, 3]

    j = which(edges$from == f & edges$to == t & edges$reg_type == r)
    if (length(j) > 0) {
      edges$color[j] = both_col
      edges$width[j] = 5
    } else {
      edges = rbind(edges, data.frame(
        from = f,
        to   = t,
        reg_type = r,
        arrows.to.type = arrow_type_from(r),
        width  = 4,
        color  = only2_col,
        stringsAsFactors = FALSE
      ))
    }
  }

  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = TRUE)
}

gt1 = read.csv(file.path(DATA_DIR, "fig2B_ground_truth_1.csv"))
gt2 = read.csv(file.path(DATA_DIR, "fig2B_ground_truth_2.csv"))

plot_network_both_signed(gt1,gt2)


# figure 2C
metric <- "AUPRC"
x      <- c(27, 30, 33, 36, 39, 42)
true_edges <- 18

cols <- c(
  `NetDes-Duo` = "#B07AA1",
  NetDes   = "#F28E2B",
  ppcor         = "#E15759",
  SCODE         = "#76B7B2",
  Random        = "grey35",
  GENIE3 = "#59A14F",
  DeepSEM = "#FF9DA7",
  Random = "grey35"
)

dfs <- list(
  `NetDes-Duo`   = read.csv(file.path(DATA_DIR, "fig2C_auprc_NetDes_Overlap.csv")),
  NetDes         = read.csv(file.path(DATA_DIR, "fig2C_auprc_NetDes.csv")),
  GENIE3         = read.csv(file.path(DATA_DIR, "fig2C_auprc_GENIE3.csv")),
  ppcor          = read.csv(file.path(DATA_DIR, "fig2C_auprc_ppcor.csv")),
  SCODE          = read.csv(file.path(DATA_DIR, "fig2C_auprc_SCODE.csv")),
  DeepSEM        = read.csv(file.path(DATA_DIR, "fig2C_auprc_DeepSEM.csv"))
)


collapse_networks <- function(df, method, metric = "AUPRC") {
  out <- df %>%
    group_by(seed, total_edges) %>%
    summarise(
      val = mean(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )

  out$Method <- method
  out <- out[, c("Method", "seed", "total_edges", "val")]

  return(out)
}

all <- do.call(rbind, Map(collapse_networks, dfs, names(dfs), MoreArgs = list(metric = metric)))
all <- all[all$total_edges %in% x, ]
all$Method <- factor(all$Method, levels = names(dfs))

methods <- levels(all$Method)
x_chr   <- as.character(x)

mu   <- matrix(NA_real_, nrow = length(methods), ncol = length(x), dimnames = list(methods, x_chr))
sdv  <- mu
nmat <- mu

for (m in methods) {
  for (j in seq_along(x)) {
    te <- x[j]
    v  <- all$val[all$Method == m & all$total_edges == te]
    mu[m, j]   <- mean(v, na.rm = TRUE)
    sdv[m, j]  <- sd(v, na.rm = TRUE)
    nmat[m, j] <- sum(!is.na(v))
  }
}
sem <- sdv / sqrt(40)

mu   <- rbind(mu,  Random = true_edges / x)
sem  <- rbind(sem, Random = rep(0, length(x)))

ymin <- 0
ymax <- 1

pdf(file.path(OUT_DIR, "fig2C.pdf"), width = 8, height = 6)

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
              ylab = paste0(metric),
              main = "All Interactions",
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


# figure 2D
x          <- c(27, 30, 33, 36, 39, 42)
x_chr      <- as.character(x)
true_edges <- 18

cols <- c(
  `NetDes-Duo` = "#B07AA1",
  NetDes   = "#F28E2B",
  ppcor         = "#E15759",
  SCODE         = "#76B7B2",
  Random        = "grey35"
)

dfs <- list(
  `NetDes-Duo`   = read.csv(file.path(DATA_DIR, "fig2D_auprc_typed_NetDes_Overlap.csv")),
  NetDes         = read.csv(file.path(DATA_DIR, "fig2D_auprc_typed_NetDes.csv")),
  ppcor          = read.csv(file.path(DATA_DIR, "fig2D_auprc_typed_ppcor.csv")),
  SCODE          = read.csv(file.path(DATA_DIR, "fig2D_auprc_typed_SCODE.csv"))
)

collapse_networks <- function(df, method, metric) {
  out <- df %>%
    group_by(seed, total_edges) %>%
    summarise(val = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")
  out$Method <- method
  out[, c("Method", "seed", "total_edges", "val")]
}
build_mu_sem <- function(metric) {
  all <- do.call(rbind, Map(collapse_networks, dfs, names(dfs), MoreArgs = list(metric = metric)))
  all <- all[all$total_edges %in% x, ]
  all$Method <- factor(all$Method, levels = names(dfs))
  methods <- levels(all$Method)

  mu   <- matrix(NA_real_, nrow=length(methods), ncol=length(x), dimnames=list(methods, x_chr))
  sdv  <- mu
  nmat <- mu

  for (m in methods) {
    for (j in seq_along(x)) {
      v <- all$val[all$Method == m & all$total_edges == x[j]]
      mu[m, j]   <- mean(v, na.rm=TRUE)
      sdv[m, j]  <- sd(v,   na.rm=TRUE)
      nmat[m, j] <- sum(!is.na(v))
    }
  }
  sem <- sdv / sqrt(nmat)

  # random baseline: n_act / (n_act + n_false)
  gt_s <- read_csv(file.path(DATA_DIR, "fig2B_ground_truth_1.csv"), show_col_types=FALSE)
  gt_s <- gt_s[gt_s$Target != "gene_4", ]
  n_act <- sum(gt_s$Interaction == 1)
  n_inh <- sum(gt_s$Interaction == 2)

  if (metric == "AUPRC_act") {
    # pool = n_act true + (total_edges - n_act - n_inh) false
    random_base <- n_act / (n_act + (x - n_act - n_inh))
  } else {
    random_base <- n_inh / (n_inh + (x - n_act - n_inh))
  }

  mu  <- rbind(mu,  Random = random_base)
  sem <- rbind(sem, Random = rep(0, length(x)))

  list(mu = mu, sem = sem)
}

ms_act <- build_mu_sem("AUPRC_act")
ms_inh <- build_mu_sem("AUPRC_inh")

plot_bar_panel <- function(ms, title) {
  bp <- barplot(ms$mu,
                beside     = TRUE,
                col        = cols[rownames(ms$mu)],
                border     = "black",
                ylim       = c(0, 1),
                names.arg  = x_chr,
                xlab       = "Initial Network Size",
                ylab       = "Average AUPRC",
                main       = title,
                xpd        = FALSE)
  legend("topright", legend = rownames(ms$mu),
         fill = cols[rownames(ms$mu)], bty = "n", cex = 1.05)
  arrows(x0 = bp, y0 = ms$mu - ms$sem,
         x1 = bp, y1 = ms$mu + ms$sem,
         angle = 90, code = 3, length = 0.04, lwd = 1.2)
}

pdf(file.path(OUT_DIR, "fig2D_activating.pdf"), width = 8, height = 6)
plot_bar_panel(ms_act, "Activating Interactions")
dev.off()

pdf(file.path(OUT_DIR, "fig2D_inhibiting.pdf"), width = 8, height = 6)
plot_bar_panel(ms_inh, "Inhibiting Interactions")
dev.off()
