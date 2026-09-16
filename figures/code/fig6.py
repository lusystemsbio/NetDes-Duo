"""Figure 6 - Python panels (6C, 6D); 6A/6B are in fig6.R."""

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


# figure 6C naive (crop of the S8 grid)
tf11 = pd.read_csv(os.path.join(DATA_DIR, 'naive_expression_data.csv'))
genes = pd.read_csv(os.path.join(DATA_DIR, 'naive_genes_expressed.csv'))
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
        file_name = os.path.join(OUT_DIR, "fig6C_naive_" + gene_list[i] + ".pdf")
        print(file_name)
        dri_genes = [gene_list[i]]
    else:
        sorted_name = sorted((gene_list[i], gene_list[j]))
        file_name = os.path.join(OUT_DIR, "fig6C_naive.pdf")
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


# figure 6C tumor-bearing (crop of the S9 grid)
tf11 = pd.read_csv(os.path.join(DATA_DIR, 'cancer_expression_data.csv'))
genes = pd.read_csv(os.path.join(DATA_DIR, 'cancer_genes_expressed.csv'))
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
        file_name = os.path.join(OUT_DIR, "fig6C_cancer_" + gene_list[i] + ".pdf")
        print(file_name)
        dri_genes = [gene_list[i]]
    else:
        sorted_name = sorted((gene_list[i], gene_list[j]))
        file_name = os.path.join(OUT_DIR, "fig6C_cancer.pdf")
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


# figure 6D
NAIVE_F2   = joblib.load(os.path.join(DATA_DIR, "naive_f2.joblib"))
NAIVE_DIF2 = joblib.load(os.path.join(DATA_DIR, "naive_dif2.joblib"))
NAIVE_F1   = joblib.load(os.path.join(DATA_DIR, "naive_f1.joblib"))
NAIVE_DIF1 = joblib.load(os.path.join(DATA_DIR, "naive_dif1.joblib"))

CANCER_F2   = joblib.load(os.path.join(DATA_DIR, "cancer_f2.joblib"))
CANCER_DIF2 = joblib.load(os.path.join(DATA_DIR, "cancer_dif2.joblib"))
CANCER_F1   = joblib.load(os.path.join(DATA_DIR, "cancer_f1.joblib"))
CANCER_DIF1 = joblib.load(os.path.join(DATA_DIR, "cancer_dif1.joblib"))

# =========================================================
# HELPERS
# =========================================================
def prep_pairs(f2, dif2):
    """
    Pair deletions: merge f2 and dif2 on edge identity.
    Keep only positive MSE values <= 10000.
    """
    f2_tmp = f2[["gene1", "gene2", "position1", "position2", "MSE"]].rename(
        columns={"MSE": "mse_f"}
    )

    dif2_tmp = dif2[["gene1", "gene2", "position1", "position2", "MSE"]].rename(
        columns={"MSE": "mse_dif"}
    )

    out = f2_tmp.merge(
        dif2_tmp,
        on=["gene1", "gene2", "position1", "position2"],
        how="inner"
    )

    out = out[
        (out["mse_f"] > 0) &
        (out["mse_dif"] > 0) &
        (out["mse_f"] <= 10000) &
        (out["mse_dif"] <= 10000)
    ].copy()

    out["kind"] = "Pair"
    out["label"] = out["gene1"].astype(str) + " + " + out["gene2"].astype(str)
    return out[["mse_f", "mse_dif", "kind", "label"]]


def prep_singles(f1, dif1):
    """
    Single deletions: merge f1 and dif1 on gene identity.
    Keep only positive MSE values <= 10000.
    """
    f1_tmp = f1[["gene", "MSE"]].rename(columns={"MSE": "mse_f"})
    dif1_tmp = dif1[["gene", "MSE"]].rename(columns={"MSE": "mse_dif"})

    out = f1_tmp.merge(dif1_tmp, on="gene", how="inner")

    out = out[
        (out["mse_f"] > 0) &
        (out["mse_dif"] > 0) &
        (out["mse_f"] <= 10000) &
        (out["mse_dif"] <= 10000)
    ].copy()

    out["kind"] = "Single"
    out["label"] = out["gene"].astype(str)
    return out[["mse_f", "mse_dif", "kind", "label"]]


def add_condition(df, condition_name):
    df = df.copy()
    df["condition"] = condition_name
    return df


# =========================================================
# PREP DATA
# =========================================================
naive_pairs   = add_condition(prep_pairs(NAIVE_F2, NAIVE_DIF2), "Naive")
naive_singles = add_condition(prep_singles(NAIVE_F1, NAIVE_DIF1), "Naive")

cancer_pairs   = add_condition(prep_pairs(CANCER_F2, CANCER_DIF2), "Cancer")
cancer_singles = add_condition(prep_singles(CANCER_F1, CANCER_DIF1), "Cancer")

naive_all  = pd.concat([naive_pairs, naive_singles], ignore_index=True)
cancer_all = pd.concat([cancer_pairs, cancer_singles], ignore_index=True)
all_data   = pd.concat([naive_all, cancer_all], ignore_index=True)

# =========================================================
# GLOBAL AXIS LIMITS
# Use BOTH conditions together, then apply same limits to both panels
# Round outward to exact powers of 10
# =========================================================
x_min_raw = all_data["mse_f"].min()
x_max_raw = all_data["mse_f"].max()
y_min_raw = all_data["mse_dif"].min()
y_max_raw = all_data["mse_dif"].max()

x_lo = 10 ** np.floor(np.log10(x_min_raw))
x_hi = 10 ** np.ceil(np.log10(x_max_raw))
y_lo = 10 ** np.floor(np.log10(y_min_raw))
y_hi = 10 ** np.ceil(np.log10(y_max_raw))

print("X limits:", x_lo, "to", x_hi)
print("Y limits:", y_lo, "to", y_hi)

# =========================================================
# PLOT STYLING
# "Double as big"
# =========================================================
TITLE_FS  = 32
LABEL_FS  = 28
TICK_FS   = 24
LEGEND_FS = 24
ANNOT_FS  = 20

PAIR_SIZE   = 110
SINGLE_SIZE = 150

PAIR_COLOR   = "tab:blue"
SINGLE_COLOR = "tab:orange"

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

for ax, condition_df, title in zip(
    axes,
    [naive_all, cancer_all],
    ["Naive", "Tumor-Bearing"]
):
    pair_df = condition_df[condition_df["kind"] == "Pair"]
    single_df = condition_df[condition_df["kind"] == "Single"]

    ax.scatter(
        pair_df["mse_f"],
        pair_df["mse_dif"],
        s=PAIR_SIZE,
        color=PAIR_COLOR,
        marker="o",
        alpha=0.8,
        label="Pair"
    )

    ax.scatter(
        single_df["mse_f"],
        single_df["mse_dif"],
        s=SINGLE_SIZE,
        color=SINGLE_COLOR,
        marker="^",
        alpha=0.9,
        label="Single"
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    ax.set_title(title, fontsize=TITLE_FS)

    ax.set_xlabel("MSE: Forward vs Experimental", fontsize=LABEL_FS)
    ax.set_ylabel("MSE: Forward vs Backward", fontsize=LABEL_FS)

    ax.tick_params(axis="both", which="major", labelsize=TICK_FS, length=9, width=1.5)
    ax.tick_params(axis="both", which="minor", length=5, width=1.2)

    # label the 2 best pairs and best single (lowest mse_f), placed off the markers
    sel = pd.concat([pair_df.assign(_m=pair_df["mse_f"].astype(float)).nsmallest(2, "_m"),
                     single_df.assign(_m=single_df["mse_f"].astype(float)).nsmallest(1, "_m")])
    fig.canvas.draw()
    pts = ax.transData.transform(condition_df[["mse_f", "mse_dif"]].astype(float).values)
    bb, placed = ax.get_window_extent(), []
    for _, r in sel.iterrows():
        anc = ax.transData.transform((float(r["mse_f"]), float(r["mse_dif"])))
        probe = ax.annotate(r["label"], (.5, .5), xycoords="axes fraction",
                            fontsize=ANNOT_FS, fontweight="bold")
        e = probe.get_window_extent(fig.canvas.get_renderer()); w, h = e.width, e.height
        probe.remove()
        gx, gy = (v.ravel() for v in np.meshgrid(
            np.linspace(bb.x0 + 4, bb.x1 - w - 4, 70), np.linspace(bb.y0 + 4, bb.y1 - h - 4, 54)))
        d = np.hypot(gx + w / 2 - anc[0], gy + h / 2 - anc[1])
        cover = ((pts[:, :1] > gx - 10) & (pts[:, :1] < gx + w + 10) &
                 (pts[:, 1:] > gy - 10) & (pts[:, 1:] < gy + h + 10)).sum(0)
        score = 500 * cover + 0.9 * d + np.where((d < 26) | (d > 300), 1e6, 0)
        for p in placed:
            score += 900 * ((gx - 20 < p[2]) & (gx + w + 20 > p[0]) &
                            (gy - 20 < p[3]) & (gy + h + 20 > p[1]))
        k = score.argmin()
        placed.append((gx[k] - 10, gy[k] - 10, gx[k] + w + 10, gy[k] + h + 10))
        ax.annotate(r["label"], (float(r["mse_f"]), float(r["mse_dif"])),
                    xytext=((gx[k] - anc[0]) * 72 / fig.dpi, (gy[k] - anc[1]) * 72 / fig.dpi),
                    textcoords="offset points", fontsize=ANNOT_FS, fontweight="bold",
                    ha="left", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=1.2, color="black", shrinkA=2, shrinkB=4))

# One legend on the right
handles, labels = axes[0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="center left",
    bbox_to_anchor=(0.91, 0.5),
    frameon=False,
    fontsize=LEGEND_FS,
    markerscale=1.4
)

fig.subplots_adjust(right=0.88, wspace=0.32)
plt.savefig(os.path.join(OUT_DIR, "fig6D.pdf"), bbox_inches="tight")
plt.show()
