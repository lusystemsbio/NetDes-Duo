library(readr)
library(ggplot2)
library(visNetwork)

DATA_DIR <- "../data"
OUT_DIR  <- "../plots"

# figure 4A
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

naive_genes = c("Stat4", "Nfkb1", "Rel")
cancer_genes = c("Sp3", "Relb", "Jund")

delete_genes_from_tflinks <- function(tflinks, genes) {
  tflinks[!(tflinks$Source %in% genes | tflinks$Target %in% genes), , drop = FALSE]
}

naive_signed <- read_csv(file.path(DATA_DIR, "naive_final_network_signed.csv"), show_col_types = FALSE)[, c("Source","Target", "Interaction")]
cancer_signed <- read_csv(file.path(DATA_DIR, "cancer_final_network_signed.csv"), show_col_types = FALSE)[, c("Source","Target", "Interaction")]

naive_signed  <- delete_genes_from_tflinks(naive_signed,  naive_genes)
cancer_signed <- delete_genes_from_tflinks(cancer_signed, cancer_genes)
print(cancer_signed,n=50)

plot_network_both_signed(cancer_signed, naive_signed)


# figure 4B
tfs <- c("Cebpa","Cebpb","Ctnnb1","Egr1","Ets1","Ets2",
         "Fos","Hif1a","Jun","Nfkb2","Per1","Stat1")

# Edge sets, sign ignored: just the (source -> target) ordered pairs
naive_edges  <- unique(paste(naive_signed$Source,  naive_signed$Target,  sep = "->"))
cancer_edges <- unique(paste(cancer_signed$Source, cancer_signed$Target, sep = "->"))

# Full 12x12 grid of ordered pairs
grid <- expand.grid(source = tfs, target = tfs,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$key  <- paste(grid$source, grid$target, sep = "->")
grid$in_n <- grid$key %in% naive_edges
grid$in_c <- grid$key %in% cancer_edges

grid$status <- with(grid, ifelse( in_n &  in_c, "Both",
                                  ifelse( in_n & !in_c, "Naive only",
                                          ifelse(!in_n &  in_c, "Tumor-Bearing only", "Neither"))))
grid$status <- factor(grid$status,
                      levels = c("Both","Naive only","Tumor-Bearing only","Neither"))

# Axis ordering: source left->right; target top->bottom (Cebpa at top)
grid$source <- factor(grid$source, levels = tfs)
grid$target <- factor(grid$target, levels = rev(tfs))

p_grid = ggplot(grid, aes(source, target, fill = status)) +
  geom_tile(color = "grey80", linewidth = 0.4) +
  scale_fill_manual(values = c("Both"       = "#0000FF",
                               "Naive only" = "#00640080",
                               "Tumor-Bearing only" = "#64000080",
                               "Neither"    = "white"),
                    breaks = c("Both","Naive only","Tumor-Bearing only"),
                    name = NULL, drop = FALSE)  +
  coord_equal() +
  labs(x = "Source", y = "Target") +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 18),
        axis.text.y  = element_text(color = "black", size = 14),
        panel.grid   = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = "top")
ggsave(file.path(OUT_DIR, "fig4B.pdf"), p_grid, device = cairo_pdf, width = 8, height = 6)
