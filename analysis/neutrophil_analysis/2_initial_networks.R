# Step 2 - TF enrichment on the trajectory clusters from step 1, then
# RcisTarget + TRRUST + NetAct construction of the two initial GRNs.

library(Seurat)
library(dplyr)
library(BiocManager)
library(devtools)
library(NetAct)
library(edgeR)
library(devtools)
library(RcisTarget)
library(readr)
library(ggplot2)
library(sctransform)
library(mgcv)
library(scales)
library(tibble)
library(stringr)
library(parallel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(doParallel)
library(foreach)
library(Matrix)

# cisTarget rankings, TRRUST table and the step-1 Seurat objects - not committed
DB_DIR  <- "~/netdes_test/Neutrophil"
OBJ_DIR <- "~/netdes_test/Neutrophil"
OUT_DIR <- "../neutrophil_data"

nnb_clusters = readRDS(file.path(OUT_DIR, "nnb_clusters.rds"))
ncb_clusters = readRDS(file.path(OUT_DIR, "ncb_clusters.rds"))

both_tfs = readRDS(file.path(OUT_DIR, "both_tfs.rds"))
cancer_tfs = readRDS(file.path(OUT_DIR, "cancer_tfs.rds"))
naive_tfs = readRDS(file.path(OUT_DIR, "naive_tfs.rds"))
nnb_not = readRDS(file.path(OUT_DIR, "nnb_not.rds"))
ncb_not = readRDS(file.path(OUT_DIR, "ncb_not.rds"))


# NAIVE TF ENRICHMENT

library(NetAct)
library(stringr)
library(ggplot2)

#ncb = readRDS("~/netdes_test/Neutrophil/with_blood_naive_ncb.rds")
nnb = readRDS(file.path(OBJ_DIR, "with_blood_naive_nnb.rds"))


library(Seurat)
library(Matrix)

nnb = FindVariableFeatures(
  nnb,
  selection.method = "vst",
  nfeatures = 10000
)

hvg_10000 = VariableFeatures(nnb)

# 2. Get normalized expression matrix
# Seurat v5:
expr = GetAssayData(
  nnb,
  assay = DefaultAssay(nnb),
  layer = "data"
)

hvg_10000 = intersect(hvg_10000, rownames(expr))
avg_expr = Matrix::rowMeans(expr[hvg_10000, , drop = FALSE])

# 4. Select top 3,000 HVGs by average expression
back_gens = names(sort(avg_expr, decreasing = TRUE))[1:10000]

# --- databases ---
data(mDB)


# --- clean cluster gene names ---
b2 = naive_tfs
TF = intersect(naive_tfs, names(mDB))
TF.mDB = mDB[TF]
length(TF.mDB)
final.results = c(NULL)

for (k in 1:length(TF)) {
  nameB      = TF.mDB[[k]]
  listA      = nnb_clusters
  names(listA) = paste0("g", 1:length(listA))
  background = back_gens
  number     = vector(mode = "list", length = length(listA))

  for (i in 1:length(listA)) {
    a = length(background[ background %in% listA[[i]] &  background %in% nameB])
    b = length(background[ background %in% listA[[i]] & !(background %in% nameB)])
    c = length(background[ background %in% nameB      & !(background %in% listA[[i]])])
    d = length(background[!(background %in% c(listA[[i]], nameB))])
    fish.mat <- matrix(c(a, c, b, d), nrow = 2, ncol = 2)
    number[[i]] = fisher.test(fish.mat, alternative = "greater")$p.value
  }

  k1 = cbind(names(TF.mDB)[[k]], names(listA), as.matrix(number))
  final.results = rbind(final.results, k1)
}

# BH per group
q_scale_groups = list()
for (i in 1:length(nnb_clusters)) {
  name_group = paste0("g", i)
  gn = final.results[final.results[,2] == name_group, ]
  q_value = p.adjust(as.numeric(gn[,3]), method = "BH")
  q_scale_groups[[i]] = cbind(gn, q_value)
}
final.results = do.call(rbind, q_scale_groups)
unique(unlist(final.results[,1]))
naive_tfs
# --- select core TFs ---
qv  = as.numeric(final.results[,4])
tfs = final.results[qv <= 0.25, , drop = FALSE]
tfs
tfs_nnb = tfs[order(as.numeric(tfs[,4])), , drop = FALSE]

q = as.numeric(tfs_nnb[, 4])
tfs_nnb_ordered = tfs_nnb[order(q), , drop = FALSE]
tfs_nnb_unique = tfs_nnb_ordered[!duplicated(tfs_nnb_ordered[, 1]), , drop = FALSE]
tfs_nnb = tfs_nnb_unique[order(as.numeric(tfs_nnb_unique[, 4])), , drop = FALSE]
tfs_nnb

tfs_2 = unique(unlist(tfs_nnb[,1]))
tfs_2
naive_tfs
length(names(mDB))

length(naive_tfs)

length(intersect(tfs_2,naive_tfs))
length(both_tfs)
length(intersect(tfs_2, both_tfs))
both_tfs

write.csv(tfs_nnb, file.path(OUT_DIR, "tfs_nnb.csv"))


# TUMOUR-BEARING TF ENRICHMENT

library(NetAct)
library(stringr)
library(ggplot2)

#back_gens = read.csv("~/netdes_test/Neutrophil/cancer_genes_expressed.csv")[[2]]
ncb = readRDS(file.path(OBJ_DIR, "with_blood_naive_ncb.rds"))
nnb = readRDS(file.path(OBJ_DIR, "with_blood_naive_nnb.rds"))

# --- databases ---
data(mDB)


ncb = FindVariableFeatures(
  ncb,
  selection.method = "vst",
  nfeatures = 10000
)

hvg_10000 = VariableFeatures(ncb)

# 2. Get normalized expression matrix
# Seurat v5:
expr = GetAssayData(
  ncb,
  assay = DefaultAssay(ncb),
  layer = "data"
)

hvg_10000 = intersect(hvg_10000, rownames(expr))
avg_expr = Matrix::rowMeans(expr[hvg_10000, , drop = FALSE])

# 4. Select top 3,000 HVGs by average expression
back_gens = names(sort(avg_expr, decreasing = TRUE))[1:10000]

#motifRankings = importRankings("~/netdes_test/Neutrophil/mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")


# --- clean cluster gene names ---

b2 = cancer_tfs

TF = intersect(cancer_tfs, names(mDB))

TF.mDB = mDB[TF]
final.results = c(NULL)

for (k in 1:length(TF)) {
  nameB      = TF.mDB[[k]]
  listA      = ncb_clusters
  names(listA) = paste0("g", 1:length(listA))
  background = back_gens
  number     = vector(mode = "list", length = length(listA))

  for (i in 1:length(listA)) {
    a = length(background[ background %in% listA[[i]] &  background %in% nameB])
    b = length(background[ background %in% listA[[i]] & !(background %in% nameB)])
    c = length(background[ background %in% nameB      & !(background %in% listA[[i]])])
    d = length(background[!(background %in% c(listA[[i]], nameB))])
    fish.mat <- matrix(c(a, c, b, d), nrow = 2, ncol = 2)
    number[[i]] = fisher.test(fish.mat, alternative = "greater")$p.value
  }

  k1 = cbind(names(TF.mDB)[[k]], names(listA), as.matrix(number))
  final.results = rbind(final.results, k1)
}

# BH per group
q_scale_groups = list()
for (i in 1:length(ncb_clusters)) {
  name_group = paste0("g", i)
  gn = final.results[final.results[,2] == name_group, ]
  q_value = p.adjust(as.numeric(gn[,3]), method = "BH")
  q_scale_groups[[i]] = cbind(gn, q_value)
}
final.results = do.call(rbind, q_scale_groups)
cancer_tfs
# --- select core TFs ---
qv  = as.numeric(final.results[,4])
tfs = final.results[qv <= 0.25, , drop = FALSE]
tfs_ncb = tfs[order(as.numeric(tfs[,4])), , drop = FALSE]

q = as.numeric(tfs_ncb[, 4])
tfs_ncb_ordered = tfs_ncb[order(q), , drop = FALSE]
tfs_ncb_unique = tfs_ncb_ordered[!duplicated(tfs_ncb_ordered[, 1]), , drop = FALSE]
tfs_ncb = tfs_ncb_unique[order(as.numeric(tfs_ncb_unique[, 4])), , drop = FALSE]

tfs_1 = unique(unlist(tfs_ncb[,1]))
tfs_1
cancer_tfs

length(intersect(tfs_1, cancer_tfs))
length(tfs_1)

length(cancer_tfs)


length(intersect(tfs_1, both_tfs))
#intersect(tfs_2, to_add_n)
both_tfs

write.csv(tfs_ncb, file.path(OUT_DIR, "tfs_ncb.csv"))


# INITIAL NETWORK CONSTRUCTION

# --- shared helpers ---

genenames<-function(listgene){
  genes=list(length=length(listgene))
  for (i in 1:length(listgene)){
    gene <- gsub(" \\(.*\\). ", "; ", listgene[i], fixed=FALSE)
    k<- unique(unlist(strsplit(gene, "; ")))
    if(length(k)==0){
      genes[[i]]=c('NULL')}else{
        genes[[i]]=k
      }

  }
  genes}
targetnames<-function(listgene2){
  genes2=list(length=length(listgene2))
  for (i in 1:length(listgene2)){
    gene <- listgene2[i]
    k<- unique(unlist(strsplit(gene, ";")))
    if(length(k)==0){
      genes2[[i]]=c('NULL')}else{
        genes2[[i]]=k
      }
  }
  genes2
}
updatalist<-function(gene,number,tfs){
  genef=gene
  for (i in number){
    genef[[i]]=intersect(genef[[i]],tfs)
  }
  genef
}
getnrow<-function(genes,tfs){

  number=c(NULL)
  for (c in (1:length(genes))){
    if (!(length(intersect(genes[[c]],tfs))==0)){
      number=c(number,c)
    }
  }
  number
}
tflinks<-function(number,genes,genes2){
  dtotal=matrix((NA),ncol=3)
  for(i in number){
    if(length(genes[[i]])==1){
      d=matrix(NA,ncol = 3,nrow = length(genes2[[i]]))
      d[,1]=genes[[i]]
      d[,2]=genes2[[i]]
      dtotal=rbind(dtotal,d)
    }else if (length(genes[[i]])>1){ for (t in 1:length(genes[[i]])){
      d=matrix(NA,ncol = 3,nrow = length(genes2[[i]]))
      d[,1]=genes[[i]][t]
      d[,2]=genes2[[i]]
      dtotal=rbind(dtotal,d)
    }
    }
  }
  dtotal=dtotal[-1,]
  colnames(dtotal)=c('Source','Target','Interaction')
  dtotal
}
grn_score = function(network){
  val = 0
  nodes = unique(network[,1])
  for (i in 1:length(nodes)){
    num_tar = 0
    for (j in 1:length(network[,2])){
      if (network[,2][j] == nodes[i]){
        num_tar = num_tar + 1
      }
    }
    val = val + (num_tar)**(-2)
  }
  return (val)
}

#tfs = c(both_tfs_2, to_add_n[[1]])
#clusters = nnb_clusters

initial_grn = function(clusters, tfs, motifRankings){
  data(motifAnnotations_mgi_v9)
  for (i in 1:length(clusters)){
    clusters[[i]] <- intersect(clusters[[i]], colnames(motifRankings))
  }
  motiftable <- list()
  motifwgene <- list()
  for(i in 1:length(clusters)) {
    tryCatch({
      motifEnrichmentTable_wGenes = cisTarget(clusters[[i]], motifRankings,motifAnnot=motifAnnotations_mgi_v9)
      motifs_AUC = calcAUC(clusters[[i]], motifRankings)
      motifEnrichmentTable = addMotifAnnotation(motifs_AUC, motifAnnot=motifAnnotations_mgi_v9, nesThreshold = 1)
      motifEnrichmentTable_wGenes = addSignificantGenes(motifEnrichmentTable, geneSets=clusters[[i]], rankings=motifRankings, nCores=1, method="aprox")

      motiftable[[i]] <- motifEnrichmentTable
      motifwgene[[i]] <- motifEnrichmentTable_wGenes
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
    })
  }
  motifwgenetotal = do.call(rbind,motifwgene)
  motiftabletotal = do.call(rbind,motiftable)
  lim = 2
  ae.wgene = motifwgenetotal[NES > lim, ]
  ae = motiftabletotal[motiftabletotal$NES>lim,]
  TFhighConf = split(ae.wgene$TF_highConf,
                     ae$geneSet)
  listgene = TFhighConf$geneSet
  TFhighConf2 = split(ae.wgene$enrichedGenes,
                      ae$geneSet)
  listgene2 = TFhighConf2$geneSet


  #intial network with Rcis target
  genestest = genenames(listgene = listgene)
  genestest2 = targetnames(listgene2 = listgene2)
  numbertest = getnrow(genestest,tfs=tfs)
  genestest = updatalist(gene = genestest,number = numbertest,tfs = tfs)
  genestest2 = updatalist(gene = genestest2,number = numbertest,tfs = tfs)
  tflinks.results = tflinks(numbertest,genes = genestest,genes2=genestest2)
  tflinks.results = as.data.frame(tflinks.results)
  dtotal.unique = distinct(tflinks.results, .keep_all = TRUE)
  dtotal.unique[,3]=1

  #adding more connections using Trrust
  trrust_ <- read.table(file.path(DB_DIR, 'trrust_rawdata.mouse.tsv'), sep='\t', header=F)
  colnames(trrust_) = c('Source','Target','Interaction','n')
  trrust_tf=trrust_[trrust_[,1]%in%tfs,]
  trrust_tf_2=trrust_tf[trrust_tf[,2]%in%tfs,]
  trrust_tf_3=trrust_tf_2[!(trrust_tf_2[,1]==trrust_tf_2[,2]),]
  interactions=unique(trrust_tf_3[,1:2])
  inter=cbind(interactions,1)
  colnames(inter) = colnames(dtotal.unique)
  dtotal.unique  =  rbind(dtotal.unique, inter)
  dtotal.unique = unique(dtotal.unique)

  tfs_mdb2 = mDB[unique(dtotal.unique[,2])]
  tfs_mdb = mDB[unique(dtotal.unique[,2])]
  for (i in 1:length(tfs_mdb)){
    tfs_mdb2[[i]]  =  tfs_mdb[[i]][tfs_mdb[[i]]%in%tfs]
  }
  inter_NATACT = matrix(NA,nrow = 1,ncol = 3)
  inter_NATACT[1,] = c("Source","Target","Interaction")
  for (i in 1:length(tfs_mdb2)){
    if (length(tfs_mdb2[[i]])>0){
      interaction = matrix(NA,nrow = length(tfs_mdb2[[i]]),ncol = 3)
      interaction[,1] = names(tfs_mdb2[i])
      interaction[,2] = tfs_mdb2[[i]]
      interaction[,3] = 1
      inter_NATACT = rbind(inter_NATACT,interaction)
    }
  }

  colnames(inter_NATACT) = inter_NATACT[1,]
  inter_NATACT = inter_NATACT[-1,]
  inter_NATACT = as.data.frame(inter_NATACT)

  dtotal.unique = rbind(dtotal.unique,inter_NATACT)
  dtotal.unique = unique(dtotal.unique)

  #trimming network down - removing self regulation + genes with that only target/recieve
  for (i in 1:5){
    dtotal.unique <- dtotal.unique[dtotal.unique$Source != dtotal.unique$Target, ]
    dtotal.unique <- dtotal.unique[which(dtotal.unique[,1]%in%unique(dtotal.unique[,2])),]
    dtotal.unique <- dtotal.unique[which(dtotal.unique[,2]%in%unique(dtotal.unique[,1])),]
    dtotal.unique=unique(dtotal.unique)
  }
  return (dtotal.unique)

}
rank_networks = function(network_list, core_tfs, tf_list, num_added){
  grn_scores = c()
  q_scores = c()
  orders = c()
  tfs_added = c()
  for (i in 1:length(network_list)){
    network = network_list[[i]]
    extra_tfs = unique(network[,1][!network[,1] %in% core_tfs])
    temp_tfs = tf_list[!duplicated(tf_list[,1]), ]
    temp_tfs = temp_tfs[temp_tfs[,1] %in% extra_tfs,  ,drop = FALSE]
    if (length(temp_tfs[,4]) == num_added){
      q_scores = append(q_scores, -mean(0.5/log(as.numeric(temp_tfs[,4]))))
      grn_scores = append(grn_scores, grn_score(network))
      orders = append(orders, i)
      tfs_added = append(tfs_added, list(extra_tfs))
    } else {
      q_scores = append(q_scores, 100)
      grn_scores = append(grn_scores, grn_score(network))
      orders = append(orders, i)
      tfs_added = append(tfs_added, list(extra_tfs))
    }
  }
  final_scores = grn_scores + q_scores
  ret = data.frame(network_score = grn_scores, q_score = q_scores, final_score = final_scores, order_in_list = orders, tfs_added = I(tfs_added))
  ret = ret[order(ret$final_score),]
  return (ret)
}

nnb_clusters = readRDS(file.path(OUT_DIR, "nnb_clusters.rds"))
ncb_clusters = readRDS(file.path(OUT_DIR, "ncb_clusters.rds"))

both_tfs = readRDS(file.path(OUT_DIR, "both_tfs.rds"))
motifRankings = importRankings(file.path(DB_DIR, "mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"))
nnb_not = readRDS(file.path(OUT_DIR, "nnb_not.rds"))
nnb_not
#a = initial_grn(ncb_clusters, intersect(tfs_1, tfs_2), motifRankings)
#b = initial_grn(nnb_clusters, intersect(tfs_1, tfs_2), motifRankings)
#plot_network(a)
#plot_network(b)
#unique(intersect(intersect(a[,1], a[,2]), intersect(b[,1],b[,2])))
both_tfs
to_add_n = unique(nnb_not[,1])
to_add_n

both_tfs = readRDS(file.path(OUT_DIR, "both_tfs.rds"))
both_tfs

cancer_tfs = readRDS(file.path(OUT_DIR, "cancer_tfs.rds"))
naive_tfs = readRDS(file.path(OUT_DIR, "naive_tfs.rds"))


intersect(cancer_tfs, naive_tfs)

ncb_not = readRDS(file.path(OUT_DIR, "ncb_not.rds"))
to_add_c = unique(ncb_not[,1])

no_add_n = initial_grn(nnb_clusters, both_tfs, motifRankings)
no_add_c = initial_grn(ncb_clusters, both_tfs, motifRankings)
both_tfs_2 = intersect(no_add_n[,1], no_add_c[,1])
core_tfs = both_tfs_2
core_tfs
no_add_n_2 = initial_grn(nnb_clusters, both_tfs_2, motifRankings)
no_add_c_2 = initial_grn(ncb_clusters, both_tfs_2, motifRankings)

to_add_n = append(to_add_n, setdiff(no_add_n[,1], both_tfs_2))
to_add_c = append(to_add_c, setdiff(no_add_c[,1], both_tfs_2))
both_tfs_2
to_add_n
to_add_c

generate_all = function(num_add, to_add, both_tfs, clusters, num_cores){
  all_grns = list()
  motifRankings = importRankings(file.path(DB_DIR, "mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"))
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths(c("~/Rlibs", .libPaths())))
  clusterEvalQ(cl, c(library(NetAct), library(RcisTarget)))
  all_combs = combn(to_add, num_add, simplify = FALSE)
  tf_list = list()
  for (i in 1:length(all_combs)){
    tf_list[[i]] = c(both_tfs, all_combs[[i]])
  }
  all_grns = foreach(tfs = tf_list, .packages = c("NetAct", "RcisTarget", "dplyr", "tidyverse"), .export = c("initial_grn", "genenames", "targetnames", "getnrow", "updatalist", "tflinks")) %dopar% {initial_grn(clusters, tfs, motifRankings)}
  stopCluster(cl)
  good_grns = list()
  for (i in 1:length(all_grns)){
    if (all(tf_list[[i]] %in% all_grns[[i]][,1])){
      good_grns[[length(good_grns) + 1]] = all_grns[[i]]
    }
  }
  return (good_grns)
}
naive_combinations = combn(to_add_n, 1, simplify = FALSE)
for (i in 1:length(naive_combinations)){
  print(naive_combinations[[i]])
}
naive_one_networks = generate_all(1, to_add_n, both_tfs_2, nnb_clusters, 7)
saveRDS(naive_one_networks, file.path(OUT_DIR, "naive_one_networks.rds"))
naive_two_networks = generate_all(2, to_add_n, both_tfs_2, nnb_clusters, 21)
saveRDS(naive_two_networks, file.path(OUT_DIR, "naive_two_networks.rds"))
naive_three_networks = generate_all(3, to_add_n, both_tfs_2, nnb_clusters, 20)
saveRDS(naive_three_networks, file.path(OUT_DIR, "naive_three_networks.rds"))

cancer_one_networks = generate_all(1, to_add_c, both_tfs_2, ncb_clusters, 7)
saveRDS(cancer_one_networks, file.path(OUT_DIR, "cancer_one_networks.rds"))
cancer_two_networks = generate_all(2, to_add_c, both_tfs_2, ncb_clusters, 21)
saveRDS(cancer_two_networks, file.path(OUT_DIR, "cancer_two_networks.rds"))
cancer_three_networks = generate_all(3, to_add_c, both_tfs_2, ncb_clusters, 20)
saveRDS(cancer_three_networks, file.path(OUT_DIR, "cancer_three_networks.rds"))

naive_one_scores = rank_networks(naive_one_networks, both_tfs_2, nnb_not, 1)
naive_two_scores = rank_networks(naive_two_networks, both_tfs_2, nnb_not, 2)
naive_three_scores = rank_networks(naive_three_networks, both_tfs_2, nnb_not, 3)

cancer_one_scores = rank_networks(cancer_one_networks, both_tfs_2, ncb_not, 1)
cancer_two_scores = rank_networks(cancer_two_networks, both_tfs_2, ncb_not, 2)
cancer_three_scores = rank_networks(cancer_three_networks, both_tfs_2, ncb_not, 3)
naive_three_scores
naive_one_scores


cancer_three_networks = readRDS(file.path(OUT_DIR, "cancer_three_networks.rds"))
naive_three_networks = readRDS(file.path(OUT_DIR, "naive_three_networks.rds"))


#save.csv(naive_three_networks[[10]])
to_add_c
write.csv(naive_three_networks[[8]], file.path(OUT_DIR, 'naive_initial_network.csv'))
write.csv(cancer_three_networks[[10]], file.path(OUT_DIR, 'cancer_initial_network.csv'))
