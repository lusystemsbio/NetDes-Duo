# Step 1 - demultiplexing, QC, PCA, diffusion pseudotime and trajectory-shape
# gene clustering for the naive and tumour-bearing neutrophil datasets.

library(Seurat)
library(destiny)
library(ggplot2)
library(parallel)
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(destiny)
library(Matrix)

# 10x matrices from GEO GSE216387 - not committed, edit for your setup
RAW_DIR <- "/projects/lulab/alex/neutrophil_case/neutrophil_raw_data"
OUT_DIR <- "../neutrophil_data"

#Helper functions
find_coords = function(obj, expr, gene){
  #use smoothing on pseudotime to get expression
  line = qplot(obj$dpt, as.numeric(expr[gene,]), alpha = I(0)) + geom_smooth(method = "loess")
  vals = ggplot_build(line)$data[[2]]
  coords = list()
  coords[[1]] = vals$x/max(vals$x)
  coords[[2]] = vals$y
  return(coords)
}
find_slopes = function(coords){
  #find/scale slopes from list of coords
  x_vals = coords[[1]]
  y_vals = coords[[2]]
  slopes = numeric(length(x_vals) - 1)
  for (i in 1:length(slopes)){
    slopes[i] = ((y_vals[i+1] - y_vals[i])/(x_vals[i+1] - x_vals[i]))/(max(y_vals)-min(y_vals))
  }
  return(slopes)
}
find_turns = function(obj, expr, gene){
  #find turns given gene name, returning cluster name (such as 'idi')
  coords = find_coords(obj, expr, gene)
  slopes = find_slopes(coords)
  turns = vector(mode = "character", length = length(slopes))
  for (i in 1:length(slopes)){
    if (abs(slopes[i]) <= 0.15){
      turns[i] = "f"
    } else if (slopes[i] > 0.15) {
      turns[i] = "i"
    } else {
      turns[i] = "d"
    }
  }
  ans = ""
  last_attached = ""
  temp_count = 0
  for (i in 1:(length(turns)-1)){
    if (turns[i+1] != turns[i] || i == length(turns) - 1){
      if (temp_count > 6 && last_attached != turns[i]){
        ans = paste(ans, turns[i], sep = "")
        last_attached = turns[i]
      }
      temp_count = 0
    } else{
      temp_count = temp_count + 1
    }
  }
  ret = list()
  ret[[1]] = ans
  ret[[2]] = coords[[1]]
  ret[[3]] = coords[[2]]

  return(ret)
}


# NAIVE

rna_barcodes = file.path(RAW_DIR, "naive/GSM6705648_GR20015_barcodes.tsv.gz")
rna_features = file.path(RAW_DIR, "naive/GSM6705648_GR20015_features.tsv.gz")
rna_mtx      = file.path(RAW_DIR, "naive/GSM6705648_GR20015_matrix.mtx.gz")

hto_barcodes = file.path(RAW_DIR, "naive/GSM6705649_GR20016-HTO_barcodes.tsv.gz")
hto_features = file.path(RAW_DIR, "naive/GSM6705649_GR20016-HTO_features.tsv.gz")
hto_mtx      = file.path(RAW_DIR, "naive/GSM6705649_GR20016-HTO_matrix.mtx.gz")

rna_counts = ReadMtx(mtx = rna_mtx, features = rna_features, cells = rna_barcodes)
hto_counts = ReadMtx(mtx = hto_mtx, features = hto_features, cells = hto_barcodes, feature.column = 1)
colnames(rna_counts) = sub("-1$", "", colnames(rna_counts))
common = intersect(colnames(rna_counts), colnames(hto_counts))
rna_counts = rna_counts[, common, drop = FALSE]
hto_counts = hto_counts[, common, drop = FALSE]
hto_counts = hto_counts[, colnames(rna_counts), drop = FALSE]

nnb = CreateSeuratObject(counts = rna_counts)
nnb[["HTO"]] = CreateAssayObject(counts = hto_counts)

nnb = MULTIseqDemux(nnb, assay = "HTO", quantile = 0.99, autoThresh = TRUE)
nnb = subset(nnb, subset = !(MULTI_ID %in% c("Negative", "Doublet", "unmapped")))

nnb$MULTI_ID = as.character(nnb$MULTI_ID)
nnb$MULTI_ID[nnb$MULTI_ID == "MHTO4-AAAGCATTCTTCACG"] = "BM"
nnb$MULTI_ID[nnb$MULTI_ID == "MHTO3-CTTGCCGCATGTCAT"] = "Blood"
nnb$MULTI_ID[nnb$MULTI_ID == "MHTO5-CTTTGTCTTTGTGAG"] = "Lung"
nnb$MULTI_ID = factor(nnb$MULTI_ID)


nnb = subset(nnb, cells = WhichCells(nnb, expression =
                                       nCount_RNA   >= quantile(nnb$nCount_RNA,   0.01, na.rm = TRUE) &
                                       nCount_RNA   <= quantile(nnb$nCount_RNA,   0.99, na.rm = TRUE) &
                                       nFeature_RNA >= quantile(nnb$nFeature_RNA, 0.01, na.rm = TRUE) &
                                       nFeature_RNA <= quantile(nnb$nFeature_RNA, 0.99, na.rm = TRUE)
))
nnb = subset(nnb, features = rownames(nnb)[Matrix::rowSums(GetAssayData(nnb, layer="counts") > 0) >= 5])

nnb = NormalizeData(nnb)
nnb = FindVariableFeatures(nnb, nfeatures = 500)
all_genes = rownames(nnb)
nnb = ScaleData(nnb, features = all_genes)
nnb = RunPCA(nnb, features = VariableFeatures(nnb))


nnb = FindNeighbors(nnb, dims = 1:20)
nnb = FindClusters(nnb, resolution = 0.5)
nnb = subset(nnb, idents = c(0, 1, 2, 3, 4, 5, 6))

nnb = RunUMAP(nnb, features = VariableFeatures(nnb))
nnb = RunTSNE(nnb, features = VariableFeatures(nnb))


dm_n = DiffusionMap(Embeddings(nnb, "pca")[,1:5])
dpt_n = DPT(dm_n)
nnb$dpt = rank(dpt_n$dpt)
nnb$dpt = 1 - (nnb$dpt / max(nnb$dpt))

# Raw DPT, not rank-normalized
dpt_raw <- dpt_n$dpt
nnb$dpt_raw <- (dpt_raw - min(dpt_raw, na.rm = TRUE)) /
  (max(dpt_raw, na.rm = TRUE) - min(dpt_raw, na.rm = TRUE))


nnb <- FindVariableFeatures(nnb, nfeatures = 5000)
top_hvg_5000 <- VariableFeatures(nnb)
expr_mat <- GetAssayData(nnb, layer = "data")
avg_expr <- Matrix::rowMeans(expr_mat)
genes_to_cluster <- top_hvg_5000[avg_expr[top_hvg_5000] > 0.01]


traj_list <- vector("list", length(genes_to_cluster))
for (i in seq_along(genes_to_cluster)) {
  gene <- genes_to_cluster[i]
  if (i %% 100 == 0) {
    print(i)
  }
  res <- find_turns(nnb, expr_mat, gene)
  traj_list[[i]] <- list(gene = gene, cluster = res[[1]], x = res[[2]], y = res[[3]])
}

traj_clusters <- data.frame(
  gene = sapply(traj_list, `[[`, "gene"),
  cluster = sapply(traj_list, `[[`, "cluster"),
  stringsAsFactors = FALSE
)
traj_clusters$x <- lapply(traj_list, `[[`, "x")
traj_clusters$y <- lapply(traj_list, `[[`, "y")

saveRDS(traj_clusters, file.path(OUT_DIR, "naive_traj_clusters.rds"))


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


# TUMOUR-BEARING

rna_barcodes = file.path(RAW_DIR, "cancer/GSM6705650_GR20017_barcodes.tsv.gz")
rna_features = file.path(RAW_DIR, "cancer/GSM6705650_GR20017_features.tsv.gz")
rna_mtx      = file.path(RAW_DIR, "cancer/GSM6705650_GR20017_matrix.mtx.gz")

hto_barcodes = file.path(RAW_DIR, "cancer/GSM6705649_GR20018-HTO_barcodes.tsv.gz")
hto_features = file.path(RAW_DIR, "cancer/GSM6705649_GR20018-HTO_features.tsv.gz")
hto_mtx      = file.path(RAW_DIR, "cancer/GSM6705649_GR20018-HTO_matrix.mtx.gz")

rna_counts = ReadMtx(mtx = rna_mtx, features = rna_features, cells = rna_barcodes)
hto_counts = ReadMtx(mtx = hto_mtx, features = hto_features, cells = hto_barcodes, feature.column = 1)
colnames(rna_counts) = sub("-1$", "", colnames(rna_counts))
common = intersect(colnames(rna_counts), colnames(hto_counts))
rna_counts = rna_counts[, common, drop = FALSE]
hto_counts = hto_counts[, common, drop = FALSE]
hto_counts = hto_counts[, colnames(rna_counts), drop = FALSE]

ncb = CreateSeuratObject(counts = rna_counts)
ncb[["HTO"]] = CreateAssayObject(counts = hto_counts)

ncb = MULTIseqDemux(ncb, assay = "HTO", quantile = 0.99, autoThresh = TRUE)
ncb = subset(ncb, subset = !(MULTI_ID %in% c("Negative", "Doublet", "unmapped")))

ncb$MULTI_ID = as.character(ncb$MULTI_ID)
ncb$MULTI_ID[ncb$MULTI_ID == "MHTO6-TATGCTGCCACGGTA"] = "Blood"
ncb$MULTI_ID[ncb$MULTI_ID == "MHTO7-GAGTCTGCCAGTATC"] = "BM"
ncb$MULTI_ID[ncb$MULTI_ID == "MHTO8-TATAGAACGCCAGGC"] = "Lung"
ncb$MULTI_ID = factor(ncb$MULTI_ID)


ncb = subset(ncb, cells = WhichCells(ncb, expression =
                                       nCount_RNA   >= quantile(ncb$nCount_RNA,   0.01, na.rm = TRUE) &
                                       nCount_RNA   <= quantile(ncb$nCount_RNA,   0.99, na.rm = TRUE) &
                                       nFeature_RNA >= quantile(ncb$nFeature_RNA, 0.01, na.rm = TRUE) &
                                       nFeature_RNA <= quantile(ncb$nFeature_RNA, 0.99, na.rm = TRUE)
))
ncb = subset(ncb, features = rownames(ncb)[Matrix::rowSums(GetAssayData(ncb, layer="counts") > 0) >= 5])

ncb = NormalizeData(ncb)
ncb = FindVariableFeatures(ncb, nfeatures = 500)
all_genes = rownames(ncb)
ncb = ScaleData(ncb, features = all_genes)
ncb = RunPCA(ncb, features = VariableFeatures(ncb))


ncb = FindNeighbors(ncb, dims = 1:20)
ncb = FindClusters(ncb, resolution = 0.5)
ncb = subset(ncb, idents = c(0, 1, 2, 3, 4, 5, 6))

ncb = RunUMAP(ncb, features = VariableFeatures(ncb))
ncb = RunTSNE(ncb, features = VariableFeatures(ncb))


dm_n = DiffusionMap(Embeddings(ncb, "pca")[,1:5])
dpt_n = DPT(dm_n)
ncb$dpt = rank(dpt_n$dpt)
ncb$dpt = 1 - (ncb$dpt / max(ncb$dpt))

# Raw DPT, not rank-normalized
dpt_raw <- dpt_n$dpt
ncb$dpt_raw <- (dpt_raw - min(dpt_raw, na.rm = TRUE)) /
  (max(dpt_raw, na.rm = TRUE) - min(dpt_raw, na.rm = TRUE))


ncb <- FindVariableFeatures(ncb, nfeatures = 5000)
top_hvg_5000 <- VariableFeatures(ncb)
expr_mat <- GetAssayData(ncb, layer = "data")
avg_expr <- Matrix::rowMeans(expr_mat)
genes_to_cluster <- top_hvg_5000[avg_expr[top_hvg_5000] > 0.01]


traj_list <- vector("list", length(genes_to_cluster))
for (i in seq_along(genes_to_cluster)) {
  gene <- genes_to_cluster[i]
  if (i %% 100 == 0) {
    print(i)
  }
  res <- find_turns(ncb, expr_mat, gene)
  traj_list[[i]] <- list(gene = gene, cluster = res[[1]], x = res[[2]], y = res[[3]])
}

traj_clusters <- data.frame(
  gene = sapply(traj_list, `[[`, "gene"),
  cluster = sapply(traj_list, `[[`, "cluster"),
  stringsAsFactors = FALSE
)
traj_clusters$x <- lapply(traj_list, `[[`, "x")
traj_clusters$y <- lapply(traj_list, `[[`, "y")

saveRDS(traj_clusters, file.path(OUT_DIR, "cancer_traj_clusters.rds"))


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


# BOTH CONDITIONS MERGED

rna_barcodes = file.path(RAW_DIR, "naive/GSM6705648_GR20015_barcodes.tsv.gz")
rna_features = file.path(RAW_DIR, "naive/GSM6705648_GR20015_features.tsv.gz")
rna_mtx      = file.path(RAW_DIR, "naive/GSM6705648_GR20015_matrix.mtx.gz")
hto_barcodes = file.path(RAW_DIR, "naive/GSM6705649_GR20016-HTO_barcodes.tsv.gz")
hto_features = file.path(RAW_DIR, "naive/GSM6705649_GR20016-HTO_features.tsv.gz")
hto_mtx      = file.path(RAW_DIR, "naive/GSM6705649_GR20016-HTO_matrix.mtx.gz")

rna_counts = ReadMtx(mtx = rna_mtx, features = rna_features, cells = rna_barcodes)
hto_counts = ReadMtx(mtx = hto_mtx, features = hto_features, cells = hto_barcodes, feature.column = 1)
colnames(rna_counts) = sub("-1$", "", colnames(rna_counts))
common = intersect(colnames(rna_counts), colnames(hto_counts))
rna_counts = rna_counts[, common, drop = FALSE]
hto_counts = hto_counts[, common, drop = FALSE]
hto_counts = hto_counts[, colnames(rna_counts), drop = FALSE]

nnb_all = CreateSeuratObject(counts = rna_counts)
nnb_all[["HTO"]] = CreateAssayObject(counts = hto_counts)

nnb_all = MULTIseqDemux(nnb_all, assay = "HTO", quantile = 0.99, autoThresh = TRUE)
nnb_all = subset(nnb_all, subset = !(MULTI_ID %in% c("Negative", "Doublet", "unmapped")))

nnb_all$MULTI_ID = as.character(nnb_all$MULTI_ID)
nnb_all$MULTI_ID[nnb_all$MULTI_ID == "MHTO4-AAAGCATTCTTCACG"] = "BM-N"
nnb_all$MULTI_ID[nnb_all$MULTI_ID == "MHTO3-CTTGCCGCATGTCAT"] = "Blood-N"
nnb_all$MULTI_ID[nnb_all$MULTI_ID == "MHTO5-CTTTGTCTTTGTGAG"] = "Lung-N"

rna_barcodes = file.path(RAW_DIR, "cancer/GSM6705650_GR20017_barcodes.tsv.gz")
rna_features = file.path(RAW_DIR, "cancer/GSM6705650_GR20017_features.tsv.gz")
rna_mtx      = file.path(RAW_DIR, "cancer/GSM6705650_GR20017_matrix.mtx.gz")
hto_barcodes = file.path(RAW_DIR, "cancer/GSM6705649_GR20018-HTO_barcodes.tsv.gz")
hto_features = file.path(RAW_DIR, "cancer/GSM6705649_GR20018-HTO_features.tsv.gz")
hto_mtx      = file.path(RAW_DIR, "cancer/GSM6705649_GR20018-HTO_matrix.mtx.gz")

rna_counts = ReadMtx(mtx = rna_mtx, features = rna_features, cells = rna_barcodes)
hto_counts = ReadMtx(mtx = hto_mtx, features = hto_features, cells = hto_barcodes, feature.column = 1)
colnames(rna_counts) = sub("-1$", "", colnames(rna_counts))
common = intersect(colnames(rna_counts), colnames(hto_counts))
rna_counts = rna_counts[, common, drop = FALSE]
hto_counts = hto_counts[, common, drop = FALSE]
hto_counts = hto_counts[, colnames(rna_counts), drop = FALSE]

ncb_all = CreateSeuratObject(counts = rna_counts)
ncb_all[["HTO"]] = CreateAssayObject(counts = hto_counts)

ncb_all = MULTIseqDemux(ncb_all, assay = "HTO", quantile = 0.99, autoThresh = TRUE)
ncb_all = subset(ncb_all, subset = !(MULTI_ID %in% c("Negative", "Doublet", "unmapped")))

ncb_all$MULTI_ID = as.character(ncb_all$MULTI_ID)
ncb_all$MULTI_ID[ncb_all$MULTI_ID == "MHTO6-TATGCTGCCACGGTA"] = "Blood-C"
ncb_all$MULTI_ID[ncb_all$MULTI_ID == "MHTO7-GAGTCTGCCAGTATC"] = "BM-C"
ncb_all$MULTI_ID[ncb_all$MULTI_ID == "MHTO8-TATAGAACGCCAGGC"] = "Lung-C"

#merge + rerun normalize/pca
nnb_all$condition = "Naive"
ncb_all$condition = "Tumor-Bearing"
neutrophil = merge(nnb_all, ncb_all)
neutrophil = NormalizeData(neutrophil, nfeatures = 500)
neutrophil = FindVariableFeatures(neutrophil)
neutrophil = ScaleData(neutrophil)
neutrophil = RunPCA(neutrophil)
