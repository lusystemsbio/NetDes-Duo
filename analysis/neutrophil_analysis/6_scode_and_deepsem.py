"""Step 6 - SCODE, DeepSEM and base NetDes baselines for the overlap comparison.
Outputs are committed to ../neutrophil_data, so step 7 runs without either tool.
"""

import os, shutil, subprocess, sys, time
import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
import NetDesDuo

DATA_DIR     = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            "..", "neutrophil_data"))
# The four constants below point outside the repo - edit for your setup.
WORK_DIR     = "/projects/lulab/alex/neutrophil_case/network_optimization"
SCODE_R      = "/projects/lulab/alex/Ode_case/SCODE/SCODE.R"
DEEPSEM_REPO = "/projects/lulab/alex/DeepSEM"
RSCRIPT      = "/shared/EL9/explorer/R/4.4.1/bin/Rscript"
TMP_ROOT     = os.path.join(DATA_DIR, "_tmp_baselines")

D, MAXITE, N_EPOCHS = 4, 100, 120   # SCODE latent dim / iterations, DeepSEM epochs
CUTS = [round(0.5 * i, 2) for i in range(1, 201)]   # NetDes MSE cuts to sweep


def load(cond):
    """That condition's binned expression, restricted to its initial-network genes."""
    expr = pd.read_csv(os.path.join(WORK_DIR, f"{cond}_results/{cond}_expression_data.csv"))
    genes = pd.read_csv(os.path.join(WORK_DIR, f"{cond}_results/{cond}_genes_expressed.csv"))["x"].tolist()
    expr = expr.transpose()[1:100000]
    expr.index = genes
    keep = pd.read_csv(os.path.join(DATA_DIR, f"{cond}_initial_network.csv"))["Target"].unique()
    expr = expr.loc[[g for g in keep if g in expr.index]]
    pt = pd.read_csv(os.path.join(DATA_DIR, f"{cond}_pseudotime_pick.csv"))["x"].values
    return expr.astype(float), pt


def run_scode(cond, expr, pt):
    tmp = os.path.join(TMP_ROOT, f"scode_{cond}"); os.makedirs(tmp, exist_ok=True)
    fe, ft, out = (os.path.join(tmp, n) for n in ("expr.txt", "time.txt", "out"))
    np.savetxt(fe, expr.values, delimiter="\t", fmt="%.6f")
    np.savetxt(ft, np.column_stack([np.arange(len(pt)), pt]), delimiter="\t", fmt="%.6f")
    g, c = expr.shape
    r = subprocess.run([RSCRIPT, SCODE_R, fe, ft, out, str(g), str(D), str(c), str(MAXITE)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"SCODE {cond}: {r.stderr[-400:]}")
    A = np.loadtxt(os.path.join(out, "A.txt"))
    genes = list(expr.index)
    return pd.DataFrame([{"Source": genes[j], "Target": genes[i], "weight": float(A[i, j])}
                         for i in range(g) for j in range(g) if i != j])


def run_deepsem(cond, expr):
    tmp = os.path.join(TMP_ROOT, f"deepsem_{cond}"); os.makedirs(tmp, exist_ok=True)
    csv, prefix = os.path.join(tmp, "data.csv"), os.path.join(tmp, "result")
    expr.T.to_csv(csv, index=False)
    r = subprocess.run([sys.executable, "main.py", "--task", "non_celltype_GRN", "--data_file", csv,
                        "--save_name", prefix, "--setting", "test", "--n_epochs", str(N_EPOCHS)],
                       capture_output=True, text=True, cwd=DEEPSEM_REPO)
    if r.returncode != 0:
        raise RuntimeError(f"DeepSEM {cond}: {r.stderr[-400:]}")
    grn = pd.read_csv(os.path.join(prefix, "GRN_inference_result.tsv"), sep="\t")
    grn.columns = ["Source", "Target", "weight"][:len(grn.columns)]
    return grn[grn.Source != grn.Target]


def run_netdes(cond, core, target):
    """delete_int at the MSE cut whose core edge count lands closest to target."""
    net = pd.read_csv(os.path.join(DATA_DIR, f"{cond}_initial_network.csv"))[["Source", "Target", "Interaction"]]
    genes = list(pd.unique(net["Target"]))
    itest = joblib.load(os.path.join(WORK_DIR, f"{cond}_results/{cond}_int_test.joblib"))
    consis = joblib.load(os.path.join(WORK_DIR, f"{cond}_results/{cond}_consis_res.joblib"))
    best = None
    for cut in CUTS:
        g = NetDesDuo.ob_to_network(genes, NetDesDuo.delete_int(genes, net, itest, consis, cut)[2])
        n = sum(1 for s, t in zip(g.Source, g.Target) if {s, t} <= core)
        if best is None or abs(n - target) < best[0]:
            best = (abs(n - target), cut, n, g)
    _, cut, n, g = best
    print(f"  netdes   cut {cut:5} -> {len(g):3} edges, {n} in core (target {target})")
    return g


if __name__ == "__main__":
    init = {c: pd.read_csv(os.path.join(DATA_DIR, f"{c}_initial_network.csv")) for c in ("naive", "cancer")}
    core = set(init["naive"]["Target"]) & set(init["cancer"]["Target"])
    target = {c: sum(1 for s, t in zip(*[pd.read_csv(os.path.join(DATA_DIR, f"{c}_final_network_signed.csv"))[k]
                                         for k in ("Source", "Target")]) if {s, t} <= core)
              for c in ("naive", "cancer")}

    for cond in ("naive", "cancer"):
        expr, pt = load(cond)
        print(f"{cond}: {expr.shape[0]} genes x {expr.shape[1]} bins")
        for name, fn in (("scode", lambda: run_scode(cond, expr, pt)),
                         ("deepsem", lambda: run_deepsem(cond, expr))):
            t0 = time.time(); grn = fn()
            grn = grn.reindex(grn.weight.abs().sort_values(ascending=False).index)
            path = os.path.join(DATA_DIR, f"{cond}_{name}_weights.csv")
            grn.to_csv(path, index=False)
            print(f"  {name:8} {len(grn):4} weights -> {os.path.basename(path)}  ({time.time()-t0:.0f}s)")
        run_netdes(cond, core, target[cond]).to_csv(
            os.path.join(DATA_DIR, f"{cond}_netdes_network.csv"), index=False)
    shutil.rmtree(TMP_ROOT, ignore_errors=True)
