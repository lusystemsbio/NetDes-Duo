# Step 4 - AUPRC evaluation of all six methods. Per-edge scores in, per-method
# AUPRC tables out. Run from inside this directory: Rscript 4_auprc_evaluation.R

library(readr)
library(GENIE3)
library(PRROC)
library(dplyr)
library(ppcor)

DATA_DIR <- "../synthetic_data/"
OUT_DIR  <- "../synthetic_data"

edge_id <- function(df){
  paste(df$Source, "to", df$Target)
}
precision <- function(gt, chosen){
  gt_id = edge_id(gt)
  chosen_id = edge_id(chosen)

  tp = sum(chosen_id %in% gt_id)
  return (tp/length(chosen_id))
}
recall <- function(gt, chosen){
  gt_id = edge_id(gt)
  chosen_id = edge_id(chosen)

  tp = sum(chosen_id %in% gt_id)
  return (tp/length(gt_id))
}
fpr <- function(init, gt, chosen){
  init_id = edge_id(init)
  gt_id = edge_id(gt)
  chosen_id = edge_id(chosen)

  fp = length(chosen_id) - sum(chosen_id %in% gt_id)
  negs = length(init_id) - length(gt_id)
  return (fp/negs)
}
#calculate partial AUPRC given a set of initial edges, ground truth, and score for each edge (0 is deleted first)
partial_AUPRC <- function(initial_edges, ground_truth, score){
  K = seq(0, 1, 0.005)
  pr_df <- data.frame(K = K, precision = NA_real_, recall = NA_real_)

  for (i in seq_along(K)){
    chosen = score[score$Interaction >= K[i], ]
    if (nrow(chosen) > 0){
      pr_df$precision[i] = precision(ground_truth, chosen)
      pr_df$recall[i] = recall(ground_truth, chosen)
    }
  }
  pr_df <- pr_df[!is.na(pr_df$precision) & !is.na(pr_df$recall), ]

  if (nrow(pr_df) == 1) return(pr_df$precision[1])

  area_under = 0
  total_area = 0
  for (i in 1:(nrow(pr_df) - 1)){
    dx = pr_df$recall[i] - pr_df$recall[i + 1]
    area_under = area_under + dx * (pr_df$precision[i+1] + pr_df$precision[i]) / 2
    total_area = total_area + dx
  }
  ratio = area_under / total_area
  return (ratio)
}

#this function does not work
partial_AUROC <- function(initial_edges, ground_truth, score){
  K = seq(0, 1, 0.005)
  roc_df <- data.frame(K = K, tpr = NA_real_, fpr = NA_real_)

  for (i in seq_along(K)){
    chosen = score[score$Interaction >= K[i], ]
    if (nrow(chosen) >= nrow(ground_truth)){
      roc_df$tpr[i] = recall(ground_truth, chosen)
      roc_df$fpr[i] = fpr(initial_edges, ground_truth, chosen)
    }
  }
  roc_df <- roc_df[!is.na(roc_df$tpr) & !is.na(roc_df$fpr), ]

  if (nrow(roc_df) == 1) return(roc_df$tpr[1])
  if (nrow(roc_df) == 0) return(0)

  area_under = 0
  total_area = 0
  for (i in 1:(nrow(roc_df) - 1)){
    dx = roc_df$fpr[i] - roc_df$fpr[i + 1]
    area_under = area_under + dx * (roc_df$tpr[i+1] + roc_df$tpr[i]) / 2
    total_area = total_area + dx
  }

  if (total_area == 0) return(roc_df$tpr[nrow(roc_df)])

  return (area_under / total_area)
}


# UNSIGNED

netdes_signscore_eval <- function(seeds, total_edges) {
  forced_target = "gene_4"
  out = data.frame(seed = integer(0), total_edges = integer(0), network = integer(0), AUPRC = numeric(0), AUROC = numeric(0))
  for (seed in seeds) {
    base = paste0(DATA_DIR, "ODE_Case_", seed)
    for (init_edge in total_edges){
      out_dir = paste0(base, "/total_", init_edge, "_edges")
      for (num in c(1, 2)) {
        init_edges = read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types = FALSE)
        gt   = read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types = FALSE)
        tab = read_csv(paste0(out_dir, "/auprc_table_", num, "_signscore.csv"), show_col_types = FALSE)
        init_edges = init_edges[init_edges$Target != forced_target, , drop = FALSE]
        edge_cols = setdiff(colnames(tab), c("n_val", "k_val"))
        scores = vapply(edge_cols, function(col) mean(as.numeric(tab[[col]]), na.rm = TRUE), numeric(1))
        parts = strsplit(edge_cols, " to ", fixed = TRUE)
        src = vapply(parts, `[`, character(1), 1)
        tgt = vapply(parts, `[`, character(1), 2)
        score_df = data.frame(Source = src, Target = tgt, Interaction = scores, stringsAsFactors = FALSE)
        score_df = score_df[score_df$Target != forced_target, , drop = FALSE]
        aup = partial_AUPRC(init_edges, gt, score_df)
        aur = partial_AUROC(init_edges, gt, score_df)
        out = rbind(out, data.frame(seed = seed, total_edges = init_edge, network = num, AUPRC = aup, AUROC = aur))
      }
    }
  }
  return (out)
}
netdes_base_eval <- function(seeds, total_edges) {
  forced_target = "gene_4"
  out = data.frame(seed = integer(0),
                   total_edges = integer(0),
                   network = integer(0),
                   AUPRC = numeric(0))

  for (seed in seeds) {
    base = paste0(DATA_DIR, "ODE_Case_", seed)
    for (init_edge in total_edges){
      out_dir = paste0(base, "/total_", init_edge, "_edges")
      for (num in c(1, 2)) {
        init_edges = read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types = FALSE)
        gt   = read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types = FALSE)
        tab = read_csv(paste0(out_dir, "/auprc_table_base_", num, ".csv"), show_col_types = FALSE)

        init_edges = init_edges[init_edges$Target != forced_target, , drop = FALSE]

        #turn table into df format
        edge_cols = setdiff(colnames(tab), c("CV"))
        #score is fraction of networks where edge is kept
        scores = vapply(edge_cols, function(col) {
          mean(as.numeric(tab[[col]]), na.rm = TRUE)
        }, numeric(1))

        parts = strsplit(edge_cols, " to ", fixed = TRUE)
        src = vapply(parts, `[`, character(1), 1)
        tgt = vapply(parts, `[`, character(1), 2)

        score_df = data.frame(
          Source = src,
          Target = tgt,
          Interaction = scores,
          stringsAsFactors = FALSE
        )
        score_df = score_df[score_df$Target != forced_target, , drop = FALSE]
        aup = partial_AUPRC(init_edges, gt, score_df)
        aur = partial_AUROC(init_edges, gt, score_df)
        out = rbind(out, data.frame(seed = seed, total_edges = init_edge, network = num, AUPRC = aup, AUROC = aur))
      }
    }
  }
  return (out)
}
genie3_eval <- function(seeds, total_edges) {
  forced_target <- "gene_4"
  out = data.frame(seed = integer(0), total_edges = integer(0), network = integer(0), AUPRC = numeric(0), AUROC = numeric(0))

  for (seed in seeds) {
    base <- paste0(DATA_DIR, "ODE_Case_", seed)
    for (num in c(1, 2)) {
      # --- expression is shared across all total_edges ---
      expr_df = read_csv(paste0(base, "/data_", num, ".csv"), show_col_types = FALSE)

      gene_names = as.character(expr_df[[1]])
      expr = as.matrix(expr_df[, -1])
      storage.mode(expr) = "numeric"
      rownames(expr) = gene_names

      gt = read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types = FALSE)

      w = GENIE3(expr)

      # evaluate each initial network size
      for (te in total_edges) {
        out_dir = paste0(base, "/total_", te, "_edges")

        init_edges = read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types = FALSE)
        init_edges =  init_edges[init_edges$Target != forced_target, , drop = FALSE]

        score_df = init_edges[, c("Source", "Target"), drop = FALSE]
        score_df$Interaction = vapply(seq_len(nrow(score_df)), function(i) {
          s <- score_df$Source[i]
          t <- score_df$Target[i]
          w[t, s]
        }, numeric(1))

        mn = min(score_df$Interaction, na.rm = TRUE)
        mx = max(score_df$Interaction, na.rm = TRUE)
        score_df$Interaction <- (score_df$Interaction - mn) / (mx - mn)
        aup = partial_AUPRC(init_edges, gt, score_df)
        aur = partial_AUROC(init_edges, gt, score_df)

        out = rbind(out, data.frame(seed = seed, total_edges = te, network = num, AUPRC = aup, AUROC = aur))
      }
    }
  }

  return (out)
}
ppcor_eval <- function(seeds, total_edges) {
  forced_target <- "gene_4"
  out <- data.frame(seed = integer(0), total_edges = integer(0), network = integer(0), AUPRC = numeric(0), AUROC = numeric(0))

  for (seed in seeds) {
    base <- paste0(DATA_DIR, "ODE_Case_", seed)
    for (num in c(1, 2)) {
      expr_df <- read_csv(paste0(base, "/data_", num, ".csv"), show_col_types = FALSE)
      genes <- as.character(expr_df[[1]])
      expr_mat <- as.matrix(expr_df[, -1])
      storage.mode(expr_mat) <- "numeric"
      rownames(expr_mat) <- genes

      X <- as.data.frame(t(expr_mat))
      colnames(X) <- genes

      pc <- spcor(X)$estimate
      dimnames(pc) <- list(genes, genes)

      gt <- read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types = FALSE)
      # evaluate each initial network size
      for (te in total_edges) {
        out_dir <- paste0(base, "/total_", te, "_edges")

        init_edges <- read_csv(paste0(out_dir, "/initial_network_", num, ".csv"),
                               show_col_types = FALSE)
        init_edges <- init_edges[init_edges$Target != forced_target, , drop = FALSE]

        score_df <- init_edges[, c("Source", "Target"), drop = FALSE]
        score_df$Interaction <- vapply(seq_len(nrow(score_df)), function(i) {
          s <- score_df$Source[i]
          t <- score_df$Target[i]
          abs(pc[t, s])
        }, numeric(1))

        mn = min(score_df$Interaction, na.rm = TRUE)
        mx = max(score_df$Interaction, na.rm = TRUE)
        score_df$Interaction <- (score_df$Interaction - mn) / (mx - mn)
        aup <- partial_AUPRC(init_edges, gt, score_df)
        aur <- partial_AUROC(init_edges, gt, score_df)

        out <- rbind(out, data.frame(seed = seed, total_edges = te, network = num,
                                     AUPRC = aup, AUROC = aur))
      }
    }
  }
  return(out)
}
scode_eval <- function(seeds, total_edges) {
  forced_target <- "gene_4"
  out <- data.frame(seed = integer(0), total_edges = integer(0),
                    network = integer(0), AUPRC = numeric(0), AUROC = numeric(0))

  for (seed in seeds) {
    base <- paste0(DATA_DIR, "ODE_Case_", seed)
    for (te in total_edges) {
      out_dir <- paste0(base, "/total_", te, "_edges")
      for (num in c(1, 2)) {
        init_edges <- read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types = FALSE)
        gt <- read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types = FALSE)
        init_edges <- init_edges[init_edges$Target != forced_target, , drop = FALSE]

        weights_file <- paste0(out_dir, "/SCODE_weights_", num, ".csv")
        if (!file.exists(weights_file)) next

        w <- read_csv(weights_file, show_col_types = FALSE)

        # Build lookup: key = "source|target" → abs(weight)
        w_lookup <- setNames(abs(w$weight), paste0(w$source, "|", w$target))

        score_df <- init_edges[, c("Source", "Target"), drop = FALSE]
        score_df$Interaction <- vapply(seq_len(nrow(score_df)), function(i) {
          key <- paste0(score_df$Source[i], "|", score_df$Target[i])
          val <- w_lookup[key]
          if (is.na(val)) 0 else val
        }, numeric(1))

        mn <- min(score_df$Interaction, na.rm = TRUE)
        mx <- max(score_df$Interaction, na.rm = TRUE)
        if (mx > mn) {
          score_df$Interaction <- (score_df$Interaction - mn) / (mx - mn)
        }

        aup <- partial_AUPRC(init_edges, gt, score_df)
        aur <- partial_AUROC(init_edges, gt, score_df)
        out <- rbind(out, data.frame(seed = seed, total_edges = te,
                                     network = num, AUPRC = aup, AUROC = aur))
      }
    }
  }
  return(out)
}
deepsem_eval <- function(seeds, total_edges) {
  forced_target <- "gene_4"
  out <- data.frame(seed = integer(0), total_edges = integer(0),
                    network = integer(0), AUPRC = numeric(0), AUROC = numeric(0))

  for (seed in seeds) {
    base <- paste0(DATA_DIR, "ODE_Case_", seed)
    for (te in total_edges) {
      out_dir <- paste0(base, "/total_", te, "_edges")
      for (num in c(1, 2)) {
        init_edges <- read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types = FALSE)
        gt <- read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types = FALSE)
        init_edges <- init_edges[init_edges$Target != forced_target, , drop = FALSE]

        weights_file <- paste0(out_dir, "/DeepSEM_weights_", num, ".csv")
        if (!file.exists(weights_file)) next

        w <- read_csv(weights_file, show_col_types = FALSE)

        # Build lookup: key = "source|target" → abs(weight)
        w_lookup <- setNames(abs(w$weight), paste0(w$source, "|", w$target))

        score_df <- init_edges[, c("Source", "Target"), drop = FALSE]
        score_df$Interaction <- vapply(seq_len(nrow(score_df)), function(i) {
          key <- paste0(score_df$Source[i], "|", score_df$Target[i])
          val <- w_lookup[key]
          if (is.na(val)) 0 else val
        }, numeric(1))

        mn <- min(score_df$Interaction, na.rm = TRUE)
        mx <- max(score_df$Interaction, na.rm = TRUE)
        if (mx > mn) {
          score_df$Interaction <- (score_df$Interaction - mn) / (mx - mn)
        }

        aup <- partial_AUPRC(init_edges, gt, score_df)
        aur <- partial_AUROC(init_edges, gt, score_df)
        out <- rbind(out, data.frame(seed = seed, total_edges = te,
                                     network = num, AUPRC = aup, AUROC = aur))
      }
    }
  }
  return(out)
}

# run all evals
seeds = c(125,177,629,753,1419, 1434, 1477,1594,1612,
          2120,2569,3219,3646,4119,4250,4510,5042,5236,6072,6534)
total_edges = c(27, 30, 33, 36, 39, 42)
netdes_base_results = netdes_base_eval(seeds, total_edges)
genie3_results = genie3_eval(seeds, total_edges)
ppcor_results = ppcor_eval(seeds, total_edges)
netdes_signscore_results = netdes_signscore_eval(seeds, total_edges)
scode_results  <- scode_eval(seeds, total_edges)
deepsem_results <- deepsem_eval(seeds, total_edges)

# figure 2C reads these six tables (as figures/data/fig2C_auprc_<method>.csv)
write.csv(netdes_signscore_results, file.path(OUT_DIR, "auprc_NetDes_Overlap.csv"), row.names = FALSE)
write.csv(netdes_base_results,      file.path(OUT_DIR, "auprc_NetDes.csv"),         row.names = FALSE)
write.csv(genie3_results,           file.path(OUT_DIR, "auprc_GENIE3.csv"),         row.names = FALSE)
write.csv(ppcor_results,            file.path(OUT_DIR, "auprc_ppcor.csv"),          row.names = FALSE)
write.csv(scode_results,            file.path(OUT_DIR, "auprc_SCODE.csv"),          row.names = FALSE)
write.csv(deepsem_results,          file.path(OUT_DIR, "auprc_DeepSEM.csv"),        row.names = FALSE)


# SIGNED

partial_AUPRC_typed <- function(initial_edges, ground_truth, score_df, sign_val) {
  gt_signed <- ground_truth[ground_truth$Interaction == sign_val, , drop=FALSE]
  if (nrow(gt_signed) == 0) return(NA_real_)

  # remove the other sign class from candidate pool and scores
  other_val   <- if (sign_val == 1) 2 else 1
  other_ids   <- paste(ground_truth[ground_truth$Interaction == other_val, ]$Source,
                       ground_truth[ground_truth$Interaction == other_val, ]$Target)
  init_filt   <- initial_edges[!paste(initial_edges$Source, initial_edges$Target) %in% other_ids, , drop=FALSE]
  score_filt  <- score_df[!paste(score_df$Source, score_df$Target) %in% other_ids, , drop=FALSE]

  if (nrow(init_filt) == 0 || nrow(score_filt) == 0) return(NA_real_)

  K <- seq(0, 1, 0.005)
  pr_df <- data.frame(K = K, precision = NA_real_, recall = NA_real_)
  for (i in seq_along(K)) {
    chosen <- score_filt[score_filt$Interaction >= K[i], ]
    if (nrow(chosen) > 0) {
      pr_df$precision[i] <- precision(gt_signed, chosen)
      pr_df$recall[i]    <- recall(gt_signed, chosen)
    }
  }
  pr_df <- pr_df[!is.na(pr_df$precision) & !is.na(pr_df$recall), ]
  if (nrow(pr_df) == 1) return(pr_df$precision[1])

  area_under <- 0; total_area <- 0
  for (i in 1:(nrow(pr_df) - 1)) {
    dx <- pr_df$recall[i] - pr_df$recall[i + 1]
    area_under <- area_under + dx * (pr_df$precision[i+1] + pr_df$precision[i]) / 2
    total_area  <- total_area + dx
  }
  return(area_under / total_area)
}

netdes_signscore_typed_eval <- function(seeds, total_edges) {
  forced_target = "gene_4"
  out = data.frame(seed=integer(0), total_edges=integer(0), network=integer(0),
                   AUPRC_act=numeric(0), AUPRC_inh=numeric(0))
  for (seed in seeds) {
    base = paste0(DATA_DIR, "ODE_Case_", seed)
    for (init_edge in total_edges) {
      out_dir = paste0(base, "/total_", init_edge, "_edges")
      for (num in c(1, 2)) {
        init_edges = read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types=FALSE)
        gt         = read_csv(paste0(base,    "/ground_truth_",    num, ".csv"), show_col_types=FALSE)
        tab        = read_csv(paste0(out_dir, "/auprc_table_",     num, "_signscore_typed.csv"), show_col_types=FALSE)
        init_edges = init_edges[init_edges$Target != forced_target, , drop=FALSE]
        edge_cols  = setdiff(colnames(tab), c("n_val", "k_val"))
        parts = strsplit(edge_cols, " to ", fixed=TRUE)
        src = vapply(parts, `[`, character(1), 1)
        tgt = vapply(parts, `[`, character(1), 2)
        scores_act = vapply(edge_cols, function(col) mean(as.numeric(tab[[col]]) == 1, na.rm=TRUE), numeric(1))
        scores_inh = vapply(edge_cols, function(col) mean(as.numeric(tab[[col]]) == 2, na.rm=TRUE), numeric(1))
        score_act_df = data.frame(Source=src, Target=tgt, Interaction=scores_act, stringsAsFactors=FALSE)
        score_inh_df = data.frame(Source=src, Target=tgt, Interaction=scores_inh, stringsAsFactors=FALSE)
        score_act_df = score_act_df[score_act_df$Target != forced_target, , drop=FALSE]
        score_inh_df = score_inh_df[score_inh_df$Target != forced_target, , drop=FALSE]

        aup_act <- partial_AUPRC_typed(init_edges, gt, score_act_df, 1)
        aup_inh <- partial_AUPRC_typed(init_edges, gt, score_inh_df, 2)
        out = rbind(out, data.frame(seed=seed, total_edges=init_edge, network=num,
                                    AUPRC_act=aup_act, AUPRC_inh=aup_inh))
      }
    }
  }
  return(out)
}
netdes_base_typed_eval <- function(seeds, total_edges) {
  forced_target = "gene_4"
  out = data.frame(seed=integer(0), total_edges=integer(0), network=integer(0),
                   AUPRC_act=numeric(0), AUPRC_inh=numeric(0))
  for (seed in seeds) {
    base = paste0(DATA_DIR, "ODE_Case_", seed)
    for (init_edge in total_edges) {
      out_dir = paste0(base, "/total_", init_edge, "_edges")
      for (num in c(1, 2)) {
        init_edges = read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types=FALSE)
        gt         = read_csv(paste0(base,    "/ground_truth_",    num, ".csv"), show_col_types=FALSE)
        tab        = read_csv(paste0(out_dir, "/auprc_table_base_", num, "_typed.csv"), show_col_types=FALSE)
        init_edges = init_edges[init_edges$Target != forced_target, , drop=FALSE]
        edge_cols  = setdiff(colnames(tab), c("CV"))
        parts = strsplit(edge_cols, " to ", fixed=TRUE)
        src = vapply(parts, `[`, character(1), 1)
        tgt = vapply(parts, `[`, character(1), 2)
        scores_act = vapply(edge_cols, function(col) mean(as.numeric(tab[[col]]) == 1, na.rm=TRUE), numeric(1))
        scores_inh = vapply(edge_cols, function(col) mean(as.numeric(tab[[col]]) == 2, na.rm=TRUE), numeric(1))
        score_act_df = data.frame(Source=src, Target=tgt, Interaction=scores_act, stringsAsFactors=FALSE)
        score_inh_df = data.frame(Source=src, Target=tgt, Interaction=scores_inh, stringsAsFactors=FALSE)
        score_act_df = score_act_df[score_act_df$Target != forced_target, , drop=FALSE]
        score_inh_df = score_inh_df[score_inh_df$Target != forced_target, , drop=FALSE]

        aup_act <- partial_AUPRC_typed(init_edges, gt, score_act_df, 1)
        aup_inh <- partial_AUPRC_typed(init_edges, gt, score_inh_df, 2)
        out = rbind(out, data.frame(seed=seed, total_edges=init_edge, network=num,
                                    AUPRC_act=aup_act, AUPRC_inh=aup_inh))
      }
    }
  }
  return(out)
}
ppcor_typed_eval <- function(seeds, total_edges) {
  forced_target <- "gene_4"
  out <- data.frame(seed=integer(0), total_edges=integer(0), network=integer(0),
                    AUPRC_act=numeric(0), AUPRC_inh=numeric(0))
  for (seed in seeds) {
    base <- paste0(DATA_DIR, "ODE_Case_", seed)
    for (num in c(1, 2)) {
      expr_df  <- read_csv(paste0(base, "/data_", num, ".csv"), show_col_types=FALSE)
      genes    <- as.character(expr_df[[1]])
      expr_mat <- as.matrix(expr_df[, -1])
      storage.mode(expr_mat) <- "numeric"
      rownames(expr_mat) <- genes
      X <- as.data.frame(t(expr_mat))
      colnames(X) <- genes
      pc <- spcor(X)$estimate
      dimnames(pc) <- list(genes, genes)
      gt <- read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types=FALSE)
      for (te in total_edges) {
        out_dir    <- paste0(base, "/total_", te, "_edges")
        init_edges <- read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types=FALSE)
        init_edges <- init_edges[init_edges$Target != forced_target, , drop=FALSE]

        # activating score: positive partial correlation magnitude
        score_act_df <- init_edges[, c("Source", "Target"), drop=FALSE]
        score_act_df$Interaction <- vapply(seq_len(nrow(score_act_df)), function(i) {
          v <- pc[score_act_df$Target[i], score_act_df$Source[i]]
          if (v > 0) v else 0
        }, numeric(1))
        mn <- min(score_act_df$Interaction); mx <- max(score_act_df$Interaction)
        if (mx > mn) score_act_df$Interaction <- (score_act_df$Interaction - mn) / (mx - mn)

        # inhibiting score: negative partial correlation magnitude
        score_inh_df <- init_edges[, c("Source", "Target"), drop=FALSE]
        score_inh_df$Interaction <- vapply(seq_len(nrow(score_inh_df)), function(i) {
          v <- pc[score_inh_df$Target[i], score_inh_df$Source[i]]
          if (v < 0) abs(v) else 0
        }, numeric(1))
        mn <- min(score_inh_df$Interaction); mx <- max(score_inh_df$Interaction)
        if (mx > mn) score_inh_df$Interaction <- (score_inh_df$Interaction - mn) / (mx - mn)

        aup_act <- partial_AUPRC_typed(init_edges, gt, score_act_df, 1)
        aup_inh <- partial_AUPRC_typed(init_edges, gt, score_inh_df, 2)
        out <- rbind(out, data.frame(seed=seed, total_edges=te, network=num,
                                     AUPRC_act=aup_act, AUPRC_inh=aup_inh))
      }
    }
  }
  return(out)
}
scode_typed_eval <- function(seeds, total_edges) {
  forced_target <- "gene_4"
  out <- data.frame(seed=integer(0), total_edges=integer(0), network=integer(0),
                    AUPRC_act=numeric(0), AUPRC_inh=numeric(0))
  for (seed in seeds) {
    base <- paste0(DATA_DIR, "ODE_Case_", seed)
    for (te in total_edges) {
      out_dir <- paste0(base, "/total_", te, "_edges")
      for (num in c(1, 2)) {
        init_edges <- read_csv(paste0(out_dir, "/initial_network_", num, ".csv"), show_col_types=FALSE)
        gt         <- read_csv(paste0(base, "/ground_truth_", num, ".csv"), show_col_types=FALSE)
        init_edges <- init_edges[init_edges$Target != forced_target, , drop=FALSE]

        weights_file <- paste0(out_dir, "/SCODE_weights_", num, ".csv")
        if (!file.exists(weights_file)) next

        w <- read_csv(weights_file, show_col_types=FALSE)
        w_lookup <- setNames(w$weight, paste0(w$source, "|", w$target))

        # activating score: positive weight magnitude
        score_act_df <- init_edges[, c("Source", "Target"), drop=FALSE]
        score_act_df$Interaction <- vapply(seq_len(nrow(score_act_df)), function(i) {
          val <- w_lookup[paste0(score_act_df$Source[i], "|", score_act_df$Target[i])]
          if (is.na(val) || val < 0) 0 else val
        }, numeric(1))
        mn <- min(score_act_df$Interaction); mx <- max(score_act_df$Interaction)
        if (mx > mn) score_act_df$Interaction <- (score_act_df$Interaction - mn) / (mx - mn)

        # inhibiting score: negative weight magnitude
        score_inh_df <- init_edges[, c("Source", "Target"), drop=FALSE]
        score_inh_df$Interaction <- vapply(seq_len(nrow(score_inh_df)), function(i) {
          val <- w_lookup[paste0(score_inh_df$Source[i], "|", score_inh_df$Target[i])]
          if (is.na(val) || val > 0) 0 else abs(val)
        }, numeric(1))
        mn <- min(score_inh_df$Interaction); mx <- max(score_inh_df$Interaction)
        if (mx > mn) score_inh_df$Interaction <- (score_inh_df$Interaction - mn) / (mx - mn)

        aup_act <- partial_AUPRC_typed(init_edges, gt, score_act_df, 1)
        aup_inh <- partial_AUPRC_typed(init_edges, gt, score_inh_df, 2)
        out <- rbind(out, data.frame(seed=seed, total_edges=te, network=num,
                                     AUPRC_act=aup_act, AUPRC_inh=aup_inh))
      }
    }
  }
  return(out)
}

# -------------------- run --------------------
seeds = c(125,177,629,753,1419,1434,1477,1594,1612,
          2120,2569,3219,3646,4119,4250,4510,5042,5236,6072,6534)
total_edges = c(27, 30, 33, 36, 39, 42)

netdes_signscore_typed_results <- netdes_signscore_typed_eval(seeds, total_edges)
netdes_base_typed_results <- netdes_base_typed_eval(seeds, total_edges)
ppcor_typed_results      <- ppcor_typed_eval(seeds, total_edges)
scode_typed_results <- scode_typed_eval(seeds, total_edges)

# figure 2D reads these four tables (as figures/data/fig2D_auprc_typed_<method>.csv)
write.csv(netdes_signscore_typed_results, file.path(OUT_DIR, "auprc_typed_NetDes_Overlap.csv"), row.names = FALSE)
write.csv(netdes_base_typed_results,      file.path(OUT_DIR, "auprc_typed_NetDes.csv"),         row.names = FALSE)
write.csv(ppcor_typed_results,            file.path(OUT_DIR, "auprc_typed_ppcor.csv"),          row.names = FALSE)
write.csv(scode_typed_results,            file.path(OUT_DIR, "auprc_typed_SCODE.csv"),          row.names = FALSE)
