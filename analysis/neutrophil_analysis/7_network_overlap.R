# Step 7 - naive vs tumour-bearing network overlap per method, each scoring its
# condition's initial-network candidates at NetDes-Duo's edge count (figure S14).

library(GENIE3)
library(ppcor)

DATA_DIR <- "../neutrophil_data"
# step-3 run tree - not committed, edit for your setup
WORK_DIR <- "/projects/lulab/alex/neutrophil_case/network_optimization"
CONDS    <- c("naive", "cancer")
SEED     <- 1

# binned expression, transposed
load_expr <- function(cond) {
  p <- file.path(WORK_DIR, paste0(cond, "_results"))
  d <- read.csv(file.path(p, paste0(cond, "_expression_data.csv")))
  g <- read.csv(file.path(p, paste0(cond, "_genes_expressed.csv")))[["x"]]
  d <- as.data.frame(t(d))
  d <- d[2:nrow(d), , drop = FALSE]
  rownames(d) <- g
  d
}

init  <- lapply(setNames(CONDS, CONDS), function(c) read.csv(file.path(DATA_DIR, paste0(c, "_initial_network.csv"))))
final <- lapply(setNames(CONDS, CONDS), function(c) read.csv(file.path(DATA_DIR, paste0(c, "_final_network_signed.csv"))))
netdes <- lapply(setNames(CONDS, CONDS), function(c) read.csv(file.path(DATA_DIR, paste0(c, "_netdes_network.csv"))))
core  <- intersect(unique(init$naive$Target), unique(init$cancer$Target))
cat(sprintf("core TFs (%d): %s\n", length(core), paste(sort(core), collapse = ", ")))

in_core <- function(df) df[df$Source %in% core & df$Target %in% core & df$Source != df$Target,
                           c("Source", "Target"), drop = FALSE]
key <- function(df) paste(df$Source, df$Target, sep = "->")

# per condition: candidate pool, target size, and one score column per method
sel <- list()
for (cond in CONDS) {
  expr  <- load_expr(cond)
  genes <- unique(init[[cond]]$Target)
  X     <- as.matrix(expr[genes, , drop = FALSE])
  storage.mode(X) <- "numeric"

  set.seed(SEED)
  g3 <- GENIE3(X)                      # [regulator, target]
  pc <- spcor(t(X))$estimate; dimnames(pc) <- list(genes, genes)
  ext <- lapply(setNames(c("scode", "deepsem"), c("SCODE", "DeepSEM")), function(m) {
    w <- read.csv(file.path(DATA_DIR, paste0(cond, "_", m, "_weights.csv")))
    setNames(abs(w$weight), paste(w$Source, w$Target, sep = "->"))
  })

  cand <- in_core(init[[cond]])
  n_keep <- nrow(in_core(final[[cond]]))
  score <- list(
    GENIE3  = g3[cbind(cand$Source, cand$Target)],
    ppcor   = abs(pc[cbind(cand$Target, cand$Source)]),
    SCODE   = ext$SCODE[key(cand)],
    DeepSEM = ext$DeepSEM[key(cand)]
  )
  top <- lapply(score, function(s) key(cand)[order(s, decreasing = TRUE)][seq_len(n_keep)])
  top[["NetDes-Duo"]] <- key(in_core(final[[cond]]))
  top[["NetDes"]]     <- key(in_core(netdes[[cond]]))
  sel[[cond]] <- list(top = top, n_cand = nrow(cand), n_keep = n_keep)
  cat(sprintf("%-7s %3d candidates in core, keeping %d\n", cond, nrow(cand), n_keep))
}

methods <- c("NetDes-Duo", "NetDes", "GENIE3", "ppcor", "SCODE", "DeepSEM")
out <- do.call(rbind, lapply(methods, function(m) {
  a <- sel$naive$top[[m]]; b <- sel$cancer$top[[m]]
  data.frame(method = m, naive_edges = length(a), cancer_edges = length(b),
             shared = length(intersect(a, b)), union = length(union(a, b)),
             jaccard = length(intersect(a, b)) / length(union(a, b)))
}))

print(out, row.names = FALSE, digits = 4)
write.csv(out, file.path(DATA_DIR, "neutrophil_overlap_jaccard.csv"), row.names = FALSE)
