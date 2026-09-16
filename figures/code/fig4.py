"""Figure 4 - Python panels (4C, 4D); 4A/4B are in fig4.R."""

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


# figure 4C (crop of the S6 grid)
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
    ids = list(range(len(gene_list)))

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
        output_file=os.path.join(OUT_DIR, f"fig4C_{prefix}.pdf")
    )


# figure 4D
from matplotlib.lines import Line2D


naive_net = pd.read_csv(os.path.join(DATA_DIR, "naive_combined_signed.csv"))
cancer_net = pd.read_csv(os.path.join(DATA_DIR, "cancer_combined_signed.csv"))

naive_tmp = naive_net[["Source", "Target", "Interaction"]].rename(
    columns={"Interaction": "interaction_naive"}
)

cancer_tmp = cancer_net[["Source", "Target", "Interaction"]].rename(
    columns={"Interaction": "interaction_cancer"}
)

# Match rows where Source and Target are the same in both networks
merged = naive_tmp.merge(
    cancer_tmp,
    on=["Source", "Target"],
    how="inner"
)

# Remove rows where cancer interaction is zero, since ratio would be undefined
merged = merged[merged["interaction_cancer"] != 0].copy()

# Ratio: naive / cancer
merged["ratio_naive_over_cancer"] = (
    merged["interaction_naive"] / merged["interaction_cancer"]
)

# Sort smallest to largest
merged_sorted = merged.sort_values("ratio_naive_over_cancer").reset_index(drop=True)

# Optional edge label
merged_sorted["edge"] = merged_sorted["Source"] + "→" + merged_sorted["Target"]

# Color rule:
# green if both interactions are on same side of 1
# red otherwise
merged_sorted["color"] = np.where(
    (
        ((merged_sorted["interaction_naive"] > 1) & (merged_sorted["interaction_cancer"] > 1)) |
        ((merged_sorted["interaction_naive"] < 1) & (merged_sorted["interaction_cancer"] < 1))
    ),
    "green",
    "red"
)

# Plot
plt.figure(figsize=(8, 3))

plt.scatter(
    range(len(merged_sorted)),
    merged_sorted["ratio_naive_over_cancer"],
    c=merged_sorted["color"],
    alpha=0.8
)

plt.axhline(1, linestyle="--", linewidth=1, color="black")
plt.axhline(2, linestyle="--", linewidth=1, color="red")
plt.axhline(0.5, linestyle="--", linewidth=1, color="red")


plt.yscale("log")

plt.xlabel("Rank", fontsize = 16)
plt.ylabel("Ratio", fontsize = 16)
plt.title("Interaction Strength Ratios", fontsize=  18)


legend_elements = [
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Concordant edge",
        markerfacecolor="green",
        markersize=8
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Discordant edge",
        markerfacecolor="red",
        markersize=8
    )
]
plt.legend(
    handles=legend_elements,
    loc="upper center",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False
)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4D.pdf"))
plt.show()
