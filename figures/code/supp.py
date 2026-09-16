"""Supplemental Python panels (S4, S6, S8, S9, S12); the rest are in supp.R."""

import os
import sys

NETDESDUO_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plots")

sys.path.insert(0, NETDESDUO_ROOT)
import NetDesDuo
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import random
import importlib
import math
import joblib
from sklearn.metrics import mean_squared_error
importlib.reload(NetDesDuo)
import math
import matplotlib.pyplot as plt


# figure S4
n_sweep = pd.read_csv(os.path.join(DATA_DIR, "figS4_n_sweep.csv"))
n_vals          = list(n_sweep["n"])
naive_sizes     = list(n_sweep["naive_size"])
cancer_sizes    = list(n_sweep["cancer_size"])
naive_mse_vals  = np.array(n_sweep["naive_mse"])
cancer_mse_vals = np.array(n_sweep["cancer_mse"])

# --- Plot MSE vs N ---
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(n_vals, naive_mse_vals,  color='#377eb8', linewidth=2, marker='o', label='Naive')
ax.plot(n_vals, cancer_mse_vals, color='#e41a1c', linewidth=2, marker='o', label='Tumor-Bearing')

ax.set_xlabel('n', fontsize=14)
ax.set_ylabel('Average MSE ', fontsize=12)
ax.set_title('Average MSE vs n', fontsize=16, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=12)
ax.tick_params(axis='both', labelsize=12)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "figS4_mse_vs_n.pdf"))
plt.show()

# --- Plot MSE vs Size---
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(naive_sizes, naive_mse_vals,  color='#377eb8', linewidth=2, marker='o', label='Naive')
ax.plot(cancer_sizes, cancer_mse_vals, color='#e41a1c', linewidth=2, marker='o', label='Tumor-Bearing')

ax.set_xlabel('Network Size', fontsize=14)
ax.set_ylabel('Average MSE', fontsize=12)
ax.set_title('Average MSE vs Size', fontsize=16, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=12)
ax.tick_params(axis='both', labelsize=12)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "figS4_mse_vs_size.pdf"))
plt.show()

k_sweep = pd.read_csv(os.path.join(DATA_DIR, "figS4_k_sweep.csv"))
k_vals          = list(k_sweep["k"])
naive_mse_vals  = np.array(k_sweep["naive_mse"])
cancer_mse_vals = np.array(k_sweep["cancer_mse"])
jacc_vals       = list(k_sweep["jaccard"])

# --- Plot MSE vs K---
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(k_vals, naive_mse_vals/(0.01 * naive_mse_vals[0]),  color='#377eb8', linewidth=2, marker='o', label='Naive')
ax.plot(k_vals, cancer_mse_vals/(0.01 * cancer_mse_vals[0]), color='#e41a1c', linewidth=2, marker='o', label='Tumor-Bearing')

ax.set_xlabel('k', fontsize=14)
ax.set_ylabel('Percent Incrase in Average MSE From k = 0', fontsize=12)
ax.set_title('Average MSE Increase vs k', fontsize=16, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=12)
ax.tick_params(axis='both', labelsize=12)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "figS4_mse_vs_k.pdf"))
plt.show()

# --- Jaccard vs k ---
plt.figure()
plt.plot(k_vals, jacc_vals, marker="o")
plt.xlabel("k")
plt.ylabel("Jaccard Index")
plt.title(f"Naive/Cancer Similarity vs k")
plt.grid(True)
plt.savefig(os.path.join(OUT_DIR, "figS4_jaccard_vs_k.pdf"), dpi=300)
plt.show()


# figure S6
for prefix in ("naive", "cancer"):
    tf11 = pd.read_csv(os.path.join(DATA_DIR, f'{prefix}_expression_data.csv'))
    genes = pd.read_csv(os.path.join(DATA_DIR, f'{prefix}_genes_expressed.csv'))
    genes = genes['x'].to_list()
    tf11 = tf11.transpose()
    tf11 = tf11[1:100000]
    tf11.index = genes

    # restore the average pseudotime values to before log-ing
    tf12expression = tf11.apply(lambda x: (math.e**x - 1))
    tf12expression[tf12expression <= 0.001] = 0.001
    # scale expression
    tf12_scaled = tf12expression.div(tf12expression.max(axis=1), axis=0)
    tf12expression_logtarget = tf12_scaled.apply(lambda x: np.log2(x))
    tf12expression_logtarget = tf12expression_logtarget.dropna()
    tf12expression = tf12expression.dropna()

    network0 = pd.read_csv(os.path.join(DATA_DIR, f'{prefix}_initial_network.csv'))
    network1 = network0[["Source", "Target", "Interaction"]]
    gene_list = list(pd.unique(network1['Target']))
    network1.columns = ['Source', 'Target', 'Interaction']

    ob_newall = joblib.load(os.path.join(DATA_DIR, f"{prefix}_ob_newall_signed.joblib"))
    res_final = joblib.load(os.path.join(DATA_DIR, f"{prefix}_res_final_signed.joblib"))

    pseudotime_pick = pd.read_csv(os.path.join(DATA_DIR, f'{prefix}_pseudotime_pick.csv'))[['x']]
    NetDesDuo.plot_network_results(
        gene_list=gene_list,
        network=network1,
        tfexpression=tf12expression,
        ob_newall=ob_newall,
        res_final=res_final,
        tfexpression_logtarget=tf12expression_logtarget,
        pseudotime_pick=pseudotime_pick,
        output_file=os.path.join(OUT_DIR, f"figS6_{prefix}.pdf")
    )


# figure S8
tf11 = pd.read_csv(os.path.join(DATA_DIR, 'naive_expression_data.csv'))
genes = pd.read_csv(os.path.join(DATA_DIR, 'naive_genes_expressed.csv'))
genes = genes['x'].to_list()
tf11 = tf11.transpose()
tf11 = tf11[1:100000]
tf11.index = genes

tf12expression = tf11.apply(lambda x: (math.e**x - 1))
tf12expression[tf12expression <= 0.001] = 0.001
tf12_scaled = tf12expression.div(tf12expression.max(axis=1), axis=0)
tf12expression_logtarget = tf12_scaled.apply(lambda x: np.log2(x))
tf12expression_logtarget = tf12expression_logtarget.dropna()
tf12expression = tf12expression.dropna()

network0 = pd.read_csv(os.path.join(DATA_DIR, 'naive_initial_network.csv'))
network1 = network0[["Source", "Target", "Interaction"]]
gene_list = list(pd.unique(network1['Target']))
network1.columns = ['Source', 'Target', 'Interaction']

ob_newall = joblib.load(os.path.join(DATA_DIR, "naive_ob_newall_signed.joblib"))
res_final = joblib.load(os.path.join(DATA_DIR, "naive_res_final_signed.joblib"))
pseudotime_pick = pd.read_csv(os.path.join(DATA_DIR, 'naive_pseudotime_pick.csv'))[['x']]

drivers = gene_list
tf12expression_logtarget2, tf12expression_target2, gene_position,gene_list2 = NetDesDuo.process_gene_expression(tf12expression_logtarget,
                                                                                                             tf12expression,
                                                                                                             gene_list,
                                                                                                             ob_newall)
ob_newall2, res_final2,gene_position2=NetDesDuo.Fit_withoutinput(gene_list=gene_list2,
                                                              res_final=res_final,
                                                              ob_newall=ob_newall,
                                                              gene_position=gene_position,
                                                              tfexpression_logtarget=tf12expression_logtarget)

cache_naive = joblib.load(os.path.join(DATA_DIR, "naive_driving_cache_Egr1_Ets1.joblib"))
def drive (pair, tf12expression_target2=tf12expression_target2, cache_naive=cache_naive):
    i = pair[0]
    j = pair[1]
    if i == j:
        file_name = os.path.join(OUT_DIR, "figS8_" + gene_list[i] + ".pdf")
        print(file_name)
        dri_genes = [gene_list[i]]
    else:
        sorted_name = sorted((gene_list[i], gene_list[j]))
        file_name = os.path.join(OUT_DIR, "figS8.pdf")
        print(file_name)
        dri_genes = list(sorted_name)

    df3, df4,cache_naive,figs = NetDesDuo.driving_results(
        dri_genes= dri_genes,
        gene_list=gene_list,
        tfexpression_target2=tf12expression_target2,
        network=network1,
        tfexpression=tf12expression,
        pseudotime_pick=pseudotime_pick,
        res_final=res_final2,
        ob_newall=ob_newall2,
        gene_position=gene_position2,
        cached = True,
        cache = cache_naive)
    figs.savefig(file_name, dpi=300, bbox_inches='tight')
    plt.close(figs)
    return cache_naive

drive((gene_list.index("Egr1"), gene_list.index("Ets1")))


# figure S9
tf11 = pd.read_csv(os.path.join(DATA_DIR, 'cancer_expression_data.csv'))
genes = pd.read_csv(os.path.join(DATA_DIR, 'cancer_genes_expressed.csv'))
genes = genes['x'].to_list()
tf11 = tf11.transpose()
tf11 = tf11[1:100000]
tf11.index = genes

tf12expression = tf11.apply(lambda x: (math.e**x - 1))
tf12expression[tf12expression <= 0.001] = 0.001
tf12_scaled = tf12expression.div(tf12expression.max(axis=1), axis=0)
tf12expression_logtarget = tf12_scaled.apply(lambda x: np.log2(x))
tf12expression_logtarget = tf12expression_logtarget.dropna()
tf12expression = tf12expression.dropna()

network0 = pd.read_csv(os.path.join(DATA_DIR, 'cancer_initial_network.csv'))
network1 = network0[["Source", "Target", "Interaction"]]
gene_list = list(pd.unique(network1['Target']))
network1.columns = ['Source', 'Target', 'Interaction']

ob_newall = joblib.load(os.path.join(DATA_DIR, "cancer_ob_newall_signed.joblib"))
res_final = joblib.load(os.path.join(DATA_DIR, "cancer_res_final_signed.joblib"))
pseudotime_pick = pd.read_csv(os.path.join(DATA_DIR, 'cancer_pseudotime_pick.csv'))[['x']]

drivers = gene_list
tf12expression_logtarget2, tf12expression_target2, gene_position,gene_list2 = NetDesDuo.process_gene_expression(tf12expression_logtarget,
                                                                                                             tf12expression,
                                                                                                             gene_list,
                                                                                                             ob_newall)
ob_newall2, res_final2,gene_position2=NetDesDuo.Fit_withoutinput(gene_list=gene_list2,
                                                              res_final=res_final,
                                                              ob_newall=ob_newall,
                                                              gene_position=gene_position,
                                                              tfexpression_logtarget=tf12expression_logtarget)

cache_cancer = joblib.load(os.path.join(DATA_DIR, "cancer_driving_cache_Cebpb_Jun.joblib"))
def drive (pair, tf12expression_target2=tf12expression_target2, cache_cancer=cache_cancer):
    i = pair[0]
    j = pair[1]
    if i == j:
        file_name = os.path.join(OUT_DIR, "figS9_" + gene_list[i] + ".pdf")
        print(file_name)
        dri_genes = [gene_list[i]]
    else:
        sorted_name = sorted((gene_list[i], gene_list[j]))
        file_name = os.path.join(OUT_DIR, "figS9.pdf")
        print(file_name)
        dri_genes = list(sorted_name)

    df3, df4,cache_cancer,figs = NetDesDuo.driving_results(
        dri_genes= dri_genes,
        gene_list=gene_list,
        tfexpression_target2=tf12expression_target2,
        network=network1,
        tfexpression=tf12expression,
        pseudotime_pick=pseudotime_pick,
        res_final=res_final2,
        ob_newall=ob_newall2,
        gene_position=gene_position2,
        cached = True,
        cache = cache_cancer)
    figs.savefig(file_name, dpi=300, bbox_inches='tight')
    plt.close(figs)
    return cache_cancer

drive((gene_list.index("Cebpb"), gene_list.index("Jun")))

# figure S12 - only the 12 TFs present in both initial networks
core = set.intersection(*(set(pd.read_csv(os.path.join(DATA_DIR, f"{c}_initial_network.csv"))["Target"])
                          for c in ("naive", "cancer")))
for cond, label in (("naive", "Naive"), ("cancer", "Tumor-Bearing")):
    df = pd.read_csv(os.path.join(DATA_DIR, f"figS_pert_rank_{cond}.csv"))
    df = df[df["perturbation"].str.split("_").map(lambda g: set(g) <= core)].reset_index(drop=True)
    df["rank"] = df.index + 1
    df["Perturbation"] = df["perturbation"].str.replace("_", " + ", regex=False)
    fig, axes = plt.subplots(1, 2, figsize=(9, 6.6))
    fig.suptitle(f"{label}: Perturbation Strength Ranking", fontsize=17, fontweight="bold", y=0.97)
    for ax, part, heading in zip(
            axes, [df.nsmallest(20, "rank"), df.nlargest(20, "rank").sort_values("rank", ascending=False)],
            ["Largest Forward Shift", "Smallest Forward/Largest Backward Shift"]):
        ax.axis("off"); ax.set_title(heading, fontsize=12.5, pad=14, color="#333333")
        cells = [[int(r["rank"]), r["Perturbation"], f'{r["mean_shift"]:+.4f}'] for _, r in part.iterrows()]
        tbl = ax.table(cellText=cells, colLabels=["Rank", "Perturbation", "Principal Curve Shift"],
                       colWidths=[.147, .367, .416], cellLoc="center", loc="upper center")
        tbl.auto_set_font_size(False); tbl.set_fontsize(11); tbl.scale(1, 1.5)
        for (row, col), cell in tbl.get_celld().items():
            cell.set_linewidth(0.8); cell.set_edgecolor("#222222"); t = cell.get_text()
            cell.visible_edges = "TB" if row == 0 else ("B" if row == len(cells) else "")
            cell.set_facecolor("white" if row == 0 or row % 2 else "#f7f7f7")
            if row == 0: t.set_fontweight("bold")
            if col: t.set_ha("left" if col == 1 else "right"); cell.PAD = 0.04
            if col == 2 and row: t.set_fontweight("bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(OUT_DIR, f"figS12_pert_rank_{cond}.pdf"), bbox_inches="tight")
    plt.close(fig)
