"""Step 3 - SCODE and DeepSEM on the synthetic benchmark, one run per
(seed, network). Outputs are committed under ../synthetic_data, so step 4
runs without either tool. Run as a script: the worker pools start on import.
"""

import os

BASE_DIR     = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            "..", "synthetic_data"))
SCODE_R      = "/projects/lulab/alex/Ode_case/SCODE/SCODE.R"
DEEPSEM_REPO = "/projects/lulab/alex/DeepSEM"

SEEDS = [125,177,629,753,1419, 1434, 1477,1594,1612,2120,
         2569, 3219, 3646, 4119, 4250, 4510, 5042, 5236, 6072, 6534]


# SCODE
#
# Per (seed, data_idx):
#   1. Load data_{idx}.csv (7 genes x 1000 timepoints)
#   2. Write expr.txt + time.txt -> call Rscript SCODE.R
#   3. Read A.txt (7x7) -> edge list CSV -> copy into all 6 edge folders

import os, shutil, time, subprocess
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# ── Config ────────────────────────────────────────
TMP_DIR    = os.path.join(BASE_DIR, "pipeline_logs", "scode")

DATA_IDS   = [1, 2]
EDGES      = [27, 30, 33, 36, 39, 42]
G, C, D, I = 7, 1000, 4, 100
N_WORKERS  = 8

os.makedirs(TMP_DIR, exist_ok=True)

# ── Helpers ───────────────────────────────────────
def case_dir(seed):          return os.path.join(BASE_DIR, f"ODE_Case_{seed}")
def edge_dir(seed, e):       return os.path.join(case_dir(seed), f"total_{e}_edges")
def out_path(seed, e, didx): return os.path.join(edge_dir(seed, e), f"SCODE_weights_{didx}.csv")

def cleanup(*paths):
    for p in paths:
        if os.path.isdir(p):  shutil.rmtree(p, ignore_errors=True)
        elif os.path.isfile(p): os.remove(p)

# ── Core ──────────────────────────────────────────
def run_scode(seed, didx):
    """Run SCODE for one (seed, data_idx). Returns (seed, didx, status, seconds)."""

    #if all(os.path.isfile(out_path(seed, e, didx)) for e in EDGES):
    #    return (seed, didx, "skipped", 0.0)

    df_path = os.path.join(case_dir(seed), f"data_{didx}.csv")
    #if not os.path.isfile(df_path):
    #    return (seed, didx, f"FAIL: missing {df_path}", 0.0)

    data = pd.read_csv(df_path, index_col=0)
    genes = list(data.index)
    #if data.shape != (G, C):
    #    return (seed, didx, f"FAIL: shape {data.shape}", 0.0)

    # Temp files (unique per worker)
    tag        = f"{seed}_{didx}"
    tmp_expr   = os.path.join(TMP_DIR, f"{tag}_expr.txt")
    tmp_time   = os.path.join(TMP_DIR, f"{tag}_time.txt")
    tmp_outdir = os.path.join(TMP_DIR, f"{tag}_out")
    tmp_csv    = os.path.join(TMP_DIR, f"{tag}_grn.csv")
    os.makedirs(tmp_outdir, exist_ok=True)

    # expr.txt: TAB-separated, not space
    np.savetxt(tmp_expr, data.values, delimiter="\t", fmt="%.6f")

    # time.txt: two columns (index, pseudotime), TAB-separated
    times = np.column_stack([np.arange(C), np.linspace(0, 1, C)])
    np.savetxt(tmp_time, times, delimiter="\t", fmt="%.6f")
    # Run SCODE
    cmd = ["/shared/EL9/explorer/R/4.4.1/bin/Rscript", SCODE_R, tmp_expr, tmp_time, tmp_outdir,
           str(G), str(D), str(C), str(I)]
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        cleanup(tmp_expr, tmp_time, tmp_outdir)
        return (seed, didx, f"FAIL: {result.stderr[-300:]}", elapsed)

    # Read A.txt → edge list
    a_path = os.path.join(tmp_outdir, "A.txt")
    if not os.path.isfile(a_path):
        cleanup(tmp_expr, tmp_time, tmp_outdir)
        return (seed, didx, "FAIL: A.txt missing", elapsed)

    A = np.loadtxt(a_path)
    edges = [{"source": genes[j], "target": genes[i], "weight": float(A[i, j])}
             for i in range(G) for j in range(G) if i != j]
    grn = pd.DataFrame(edges).sort_values("weight", key=abs, ascending=False)
    grn.to_csv(tmp_csv, index=False)

    # Copy to all edge folders
    copied = 0
    for e in EDGES:
        dest = out_path(seed, e, didx)
        if os.path.isdir(edge_dir(seed, e)) and not os.path.isfile(dest):
            shutil.copy2(tmp_csv, dest)
            copied += 1

    cleanup(tmp_expr, tmp_time, tmp_outdir, tmp_csv)
    return (seed, didx, f"ok ({copied}/6)", elapsed)

# ── Run ───────────────────────────────────────────
valid = [s for s in SEEDS if os.path.isdir(case_dir(s))]
tasks = [(s, d) for s in valid for d in DATA_IDS]
print(f"SCODE: {len(tasks)} tasks, {len(valid)} seeds, {N_WORKERS} workers\n")

with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
    futures = {pool.submit(run_scode, s, d): (s, d) for s, d in tasks}
    for i, fut in enumerate(as_completed(futures), 1):
        seed, didx, status, sec = fut.result()
        icon = "✓" if status.startswith("ok") else "–" if status == "skipped" else "✗"
        print(f"[{icon}] seed={seed} data={didx}  {status}  ({sec:.0f}s)  [{i}/{len(tasks)}]")

print("\nDone.")


# DeepSEM

import os, subprocess, shutil, time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

TMP_ROOT     = os.path.join(BASE_DIR, "_deepsem_tmp")

DATA_INDICES = [1, 2]
EDGE_COUNTS  = [27, 30, 33, 36, 39, 42]
G, C         = 7, 1000
N_EPOCHS     = 120
N_WORKERS    = 2   # single GPU — keep low to avoid VRAM contention
# ──────────────────────────────────────────────────

os.makedirs(TMP_ROOT, exist_ok=True)

def out_path(seed, e, d):
    return os.path.join(BASE_DIR, f"ODE_Case_{seed}",
                        f"total_{e}_edges", f"DeepSEM_weights_{d}.csv")

def parse_columns(df):
    """Map DeepSEM output columns → (source, target, weight)."""
    col_map = {}
    for c in df.columns:
        cl = c.strip().lower().replace(" ", "")
        if   cl in ("gene1", "tf", "source"):                        col_map[c] = "source"
        elif cl in ("gene2", "target"):                              col_map[c] = "target"
        elif any(k in cl for k in ("weight", "importance", "score")): col_map[c] = "weight"
    if set(col_map.values()) != {"source", "target", "weight"}:
        raise ValueError(f"Can't parse columns: {list(df.columns)}")
    return df.rename(columns=col_map)[["source", "target", "weight"]]

def run_deepsem(seed, idx):
    tag = f"seed={seed} data={idx}"

    if all(os.path.isfile(out_path(seed, e, idx)) for e in EDGE_COUNTS):
        return tag, "skipped", 0.0

    src = os.path.join(BASE_DIR, f"ODE_Case_{seed}", f"data_{idx}.csv")
    if not os.path.isfile(src):
        return tag, f"FAIL: missing {src}", 0.0

    data = pd.read_csv(src, index_col=0)
    if data.shape != (G, C):
        return tag, f"FAIL: shape {data.shape}", 0.0

    # Transpose (G×C) → (C×G), write temp input
    tmp = os.path.join(TMP_ROOT, f"{seed}_{idx}")
    os.makedirs(tmp, exist_ok=True)
    tmp_csv = os.path.join(tmp, "data.csv")
    data.T.to_csv(tmp_csv, index=False)

    save_prefix = os.path.join(tmp, "result")
    cmd = ["python", "main.py",
           "--task", "non_celltype_GRN",
           "--data_file", tmp_csv,
           "--save_name", save_prefix,
           "--setting", "test",
           "--n_epochs", str(N_EPOCHS)]

    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=DEEPSEM_REPO)
    elapsed = time.time() - t0

    if proc.returncode != 0:
        shutil.rmtree(tmp, ignore_errors=True)
        return tag, f"FAIL: {proc.stderr[-500:]}", elapsed

    # Find output TSV — DeepSEM writes to save_name/GRN_inference_result.tsv
    tsv = os.path.join(save_prefix, "GRN_inference_result.tsv")
    if not os.path.isfile(tsv):
        # Fallback: search recursively
        found = []
        for root, dirs, files in os.walk(tmp):
            found += [os.path.join(root, f) for f in files if f.endswith(".tsv")]
        if not found:
            shutil.rmtree(tmp, ignore_errors=True)
            return tag, "FAIL: no .tsv output", elapsed
        tsv = found[0]

    # Parse → standardize → distribute
    try:
        grn = parse_columns(pd.read_csv(tsv, sep="\t"))
    except ValueError as e:
        shutil.rmtree(tmp, ignore_errors=True)
        return tag, f"FAIL: {e}", elapsed

    grn = grn.sort_values("weight", key=abs, ascending=False)
    result_csv = os.path.join(tmp, "grn.csv")
    grn.to_csv(result_csv, index=False)

    copied = 0
    for e in EDGE_COUNTS:
        dest = out_path(seed, e, idx)
        if os.path.isdir(os.path.dirname(dest)) and not os.path.isfile(dest):
            shutil.copy2(result_csv, dest)
            copied += 1

    shutil.rmtree(tmp, ignore_errors=True)
    return tag, f"ok ({copied}/6 copied)", elapsed

# ── RUN ───────────────────────────────────────────
valid = [s for s in SEEDS if os.path.isdir(os.path.join(BASE_DIR, f"ODE_Case_{s}"))]
tasks = [(s, d) for s in valid for d in DATA_INDICES]
print(f"DeepSEM: {len(tasks)} tasks, {N_WORKERS} workers, {len(valid)}/{len(SEEDS)} seeds valid")

with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
    futs = {pool.submit(run_deepsem, s, d): (s, d) for s, d in tasks}
    for i, fut in enumerate(as_completed(futs), 1):
        tag, status, elapsed = fut.result()
        sym = "✓" if status.startswith("ok") else "–" if status.startswith("skip") else "✗"
        print(f"[{sym}] {tag}  {status}  ({elapsed:.0f}s)  [{i}/{len(tasks)}]")

shutil.rmtree(TMP_ROOT, ignore_errors=True)
print("Done.")
