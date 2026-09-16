"""Figure 5 - Python panel (5A); 5B/5C are in fig5.R."""

import os
import sys

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plots")

import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import math
import joblib
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, FormatStrFormatter


# getting naive data
tf11 = pd.read_csv(os.path.join(DATA_DIR, 'naive_expression_data.csv'))
genes = pd.read_csv(os.path.join(DATA_DIR, 'naive_genes_expressed.csv'))
genes = genes['x'].to_list()

tf11 = tf11.transpose()
tf11 = tf11[1:100000]
tf11.index = genes

# restore the average pseudotime values to before log-ing
naive_tf12expression_all = tf11.apply(lambda x: (math.e**x - 1))
naive_tf12expression_all[naive_tf12expression_all <= 0.001] = 0.001

# scale expression
naive_tf12_scaled_all = naive_tf12expression_all.div(
    naive_tf12expression_all.max(axis=1),
    axis=0
)

naive_tf12expression_logtarget_all = naive_tf12_scaled_all.apply(lambda x: np.log2(x))
naive_tf12expression_logtarget_all = naive_tf12expression_logtarget_all.dropna()
naive_tf12expression_all = naive_tf12expression_all.dropna()

network0 = pd.read_csv(os.path.join(DATA_DIR, 'naive_initial_network.csv'))
network1 = network0[["Source", "Target", "Interaction"]]
naive_gene_list = list(pd.unique(network1['Target']))
network1.columns = ['Source', 'Target', 'Interaction']

naive_expression = naive_tf12expression_all.loc[naive_gene_list]
naive_expression_log = naive_tf12expression_logtarget_all.loc[naive_gene_list]

# getting cancer data
tf11 = pd.read_csv(os.path.join(DATA_DIR, 'cancer_expression_data.csv'))
genes = pd.read_csv(os.path.join(DATA_DIR, 'cancer_genes_expressed.csv'))
genes = genes['x'].to_list()

tf11 = tf11.transpose()
tf11 = tf11[1:100000]
tf11.index = genes

# restore the average pseudotime values to before log-ing
cancer_tf12expression_all = tf11.apply(lambda x: (math.e**x - 1))
cancer_tf12expression_all[cancer_tf12expression_all <= 0.001] = 0.001

# scale expression
cancer_tf12_scaled_all = cancer_tf12expression_all.div(
    cancer_tf12expression_all.max(axis=1),
    axis=0
)

cancer_tf12expression_logtarget_all = cancer_tf12_scaled_all.apply(lambda x: np.log2(x))
cancer_tf12expression_logtarget_all = cancer_tf12expression_logtarget_all.dropna()
cancer_tf12expression_all = cancer_tf12expression_all.dropna()

network0 = pd.read_csv(os.path.join(DATA_DIR, 'cancer_initial_network.csv'))
network1 = network0[["Source", "Target", "Interaction"]]
cancer_gene_list = list(pd.unique(network1['Target']))
network1.columns = ['Source', 'Target', 'Interaction']

cancer_expression = cancer_tf12expression_all.loc[cancer_gene_list]
cancer_expression_log = cancer_tf12expression_logtarget_all.loc[cancer_gene_list]


# figure 5A (naive)
# -----------------------------
# Naive baseline cloud + full real pseudotime trajectory
# -----------------------------
naive_baseline_cloud = pd.read_csv(os.path.join(DATA_DIR, "naive_baseline_cloud.csv"), index_col = 0)

# baseline cloud: runs x genes
X_cloud = naive_baseline_cloud.T.values

# fit PCA on the simulated baseline cloud
pca_model = PCA(n_components=2, random_state=0)
pcs_cloud = pca_model.fit_transform(X_cloud)

explained = pca_model.explained_variance_ratio_ * 100

pcs_cloud_df = pd.DataFrame(
    pcs_cloud,
    index=naive_baseline_cloud.columns,
    columns=["PC1", "PC2"]
)

# real pseudotime data in same scaled expression space
naive_scaled = naive_expression.div(naive_expression.max(axis=1), axis=0)
naive_scaled = naive_scaled.loc[naive_gene_list]

# pseudotime points: timepoints x genes
X_real = naive_scaled.T.values

# project real pseudotime trajectory into the baseline PCA axes
pcs_real = pca_model.transform(X_real)

pcs_real_df = pd.DataFrame(
    pcs_real,
    columns=["PC1", "PC2"]
)

pcs_real_df["pseudotime_index"] = np.arange(len(pcs_real_df))

# plot
plt.figure(figsize=(6, 5))

xy = np.vstack([pcs_cloud_df["PC1"].values, pcs_cloud_df["PC2"].values])
dens = gaussian_kde(xy)(xy)
order = dens.argsort()          # densest points drawn last
# log norm, scaled to this condition's own range: the naive cloud spans 360x in
# density, so on a linear scale its sparse arm collapses into one colour
sc = plt.scatter(
    pcs_cloud_df["PC1"].values[order],
    pcs_cloud_df["PC2"].values[order],
    c=dens[order],
    cmap="viridis",
    norm=LogNorm(vmin=dens.min(), vmax=dens.max()),
    s=12,
    alpha=0.35,
    label="Simulated States"
)


#plt.scatter(
#    pcs_real_df["PC1"],
#    pcs_real_df["PC2"],
#    c=pcs_real_df["pseudotime_index"],
#    s=25,
#    cmap="viridis",
#    label="Real pseudotime points"
#)

plt.scatter(
    pcs_real_df["PC1"].iloc[0],
    pcs_real_df["PC2"].iloc[0],
    s=120,
    marker="o",
    color="tab:orange",
    edgecolor="black",
    linewidth=1.2,
    label="Observed Start"
)

plt.scatter(
    pcs_real_df["PC1"].iloc[-1],
    pcs_real_df["PC2"].iloc[-1],
    s=120,
    marker="s",
    color="tab:green",
    edgecolor="black",
    linewidth=1.2,
    label="Observed End"
)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.xlabel(f"PC1 ({explained[0]:.1f}%)", fontsize = 14)
plt.ylabel(f"PC2 ({explained[1]:.1f}%)", fontsize = 14)
cb = plt.colorbar(sc, ticks=LogLocator(subs=(1, 2, 5)))
cb.ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
cb.set_label("Density", fontsize=14)
cb.ax.tick_params(labelsize=12)
plt.legend(fontsize=10)
plt.tight_layout()

plt.savefig(os.path.join(OUT_DIR, 'fig5A_naive.pdf'))
plt.show()


# figure 5A (tumor-bearing)
# -----------------------------
# Cancer baseline cloud + full real pseudotime trajectory
# -----------------------------
cancer_baseline_cloud = pd.read_csv(os.path.join(DATA_DIR, "cancer_baseline_cloud.csv"), index_col = 0)

# baseline cloud: runs x genes
X_cloud = cancer_baseline_cloud.T.values

# fit PCA on the simulated baseline cloud
pca_model = PCA(n_components=2, random_state=0)
pcs_cloud = pca_model.fit_transform(X_cloud)

explained = pca_model.explained_variance_ratio_ * 100

pcs_cloud_df = pd.DataFrame(
    pcs_cloud,
    index=cancer_baseline_cloud.columns,
    columns=["PC1", "PC2"]
)

# real pseudotime data in same scaled expression space
cancer_scaled = cancer_expression.div(cancer_expression.max(axis=1), axis=0)
cancer_scaled = cancer_scaled.loc[cancer_gene_list]

# pseudotime points: timepoints x genes
X_real = cancer_scaled.T.values

# project real pseudotime trajectory into the baseline PCA axes
pcs_real = pca_model.transform(X_real)

pcs_real_df = pd.DataFrame(
    pcs_real,
    columns=["PC1", "PC2"]
)

pcs_real_df["pseudotime_index"] = np.arange(len(pcs_real_df))

# plot
plt.figure(figsize=(6, 5))

xy = np.vstack([pcs_cloud_df["PC1"].values, pcs_cloud_df["PC2"].values])
dens = gaussian_kde(xy)(xy)
order = dens.argsort()          # densest points drawn last
# log norm, scaled to this condition's own range: the naive cloud spans 360x in
# density, so on a linear scale its sparse arm collapses into one colour
sc = plt.scatter(
    pcs_cloud_df["PC1"].values[order],
    pcs_cloud_df["PC2"].values[order],
    c=dens[order],
    cmap="viridis",
    norm=LogNorm(vmin=dens.min(), vmax=dens.max()),
    s=12,
    alpha=0.35,
    label="Simulated States"
)


#plt.scatter(
#    pcs_real_df["PC1"],
#    pcs_real_df["PC2"],
#    c=pcs_real_df["pseudotime_index"],
#    s=25,
#    cmap="viridis",
#    label="Pseudotime Points"
#)

plt.scatter(
    pcs_real_df["PC1"].iloc[0],
    pcs_real_df["PC2"].iloc[0],
    s=120,
    marker="o",
    color="tab:orange",
    edgecolor="black",
    linewidth=1.2,
    label="Observed Start"
)

plt.scatter(
    pcs_real_df["PC1"].iloc[-1],
    pcs_real_df["PC2"].iloc[-1],
    s=120,
    marker="s",
    color="tab:green",
    edgecolor="black",
    linewidth=1.2,
    label="Observed End"
)


plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.xlabel(f"PC1 ({explained[0]:.1f}%)", fontsize = 14)
plt.ylabel(f"PC2 ({explained[1]:.1f}%)", fontsize = 14)
cb = plt.colorbar(sc, ticks=LogLocator(subs=(1, 2, 5)))
cb.ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
cb.set_label("Density", fontsize=14)
cb.ax.tick_params(labelsize=12)
plt.legend(fontsize=10)

plt.savefig(os.path.join(OUT_DIR, 'fig5A_cancer.pdf'))

plt.tight_layout()
plt.show()
