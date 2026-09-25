# NetDes-Duo

**Net**work inference and optimization using **D**ynamical **e**quation **s**imulations for **Du**al-condition **o**verlap (NetDes-Duo)

NetDes-Duo is a computational method for optimizing the gene regulatory networks
(GRNs) of core transcription factors (TFs) in two related conditions at once,
based on gene expression time trajectories. It is an extension of
[NetDes](https://github.com/lusystemsbio/NetDes), which optimizes one network at
a time. 

Below is the code reproducing the benchmarking, neutrophil application, and figures for NetDes-Duo 

Paper: https://doi.org/10.64898/2026.09.16.752194

## Layout

| Folder | Content | Figure |
|---|---|---|
| `NetDesDuo.py` | the method/all functions required for application of method | — |
| `analysis/synthetic_data/` | Synthetic Networks: ODE trajectories, ground truth, initial networks | — |
| `analysis/synthetic_benchmarking/` | Synthetic benchmarking with NetDes-Duo + comparison methods | Fig.2 |
| `analysis/neutrophil_analysis/` | Naive vs Tumour-bearing neutrophil analysis| Fig.3–6 |
| `analysis/neutrophil_data/` | small committed intermediate results | — |
| `figures/` | code and inputs required to plot the figures | Fig.2-6 |


## Environment

Python 3.10+ (numpy, pandas, scipy, scikit-learn, matplotlib, joblib) and R 4.4 (Seurat, destiny, RcisTarget, NetAct, GENIE3, ppcor, PRROC, princurve, ggplot2,
visNetwork). SCODE and DeepSEM are third-party tools with their own enviroments. 

## Data

Raw data is not committed. Access at GSE217143 (samples GSM6705648, GSM6705649, and GSM6705650). 
Smaller intermediate results are commited, but larger intermediate results are not.
