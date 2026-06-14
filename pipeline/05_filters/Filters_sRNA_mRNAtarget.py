#!/usr/bin/env python3
"""
sRNA-mRNA Target Filtering and Integration Pipeline

Generated from sRNA_mRNAtarget.ipynb. Edit the USER CONFIGURATION block below,
then run:  python Filters_sRNA_mRNAtarget.py
"""

# # sRNA ↔ mRNA target filtering pipeline

# ## 1) User‑defined variables
# 
# Edit **only this cell** (paths and parameters).

# =========================
# USER CONFIGURATION
# =========================

# (A) Prediction file (CSV)
PREDICOES_CSV = "Prediction_test.csv"

# (B) DESeq2 result files (TSV) — one per strain/condition
DEG_FILES = [
    "DEG_result_A.tsv",
    "DEG_result_B.tsv",
]

# (C) Adjusted p‑value cutoff to call a DEG
PADJ_CUTOFF = 0.05

# (D) Cross‑strain consistency:
# how many strains must agree (>0 up, <0 down)
MIN_STRAINS_CONSISTENT = 2

# (E) Energy / probability outlier filters
FILTERS_OUTLIER = {
    "E_intaRNA_max": -2.44,
    "E_Rnaplex_max": -32.6,
    "E_TargetRNA3_max": -5.13,
    "Probability_TargetRNA3_min": 0.06,
    "Probability_sRNARFTarget_min": 0.40
}

# (F) STRING module nodes and edges
MODULE_NODES_FILE = "CytoscapeInput-nodes.txt"   # expected: nodeName, nodeAttr
MODULE_EDGES_FILE = "CytoscapeInput-edges.txt"   # expected: fromNode, toNode, weight

# (G) Output
OUTPUT_CSV = "filtered_weight.csv"

# ## 2) Imports and helper functions

import os
import numpy as np
import pandas as pd

def _check_exists(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

def load_predictions(csv_path):
    _check_exists(csv_path)
    df = pd.read_csv(csv_path)
    required = [
        "sRNA","Target",
        "E_intaRNA","E_Rnaplex","E_TargetRNA3",
        "Probability_TargetRNA3","Probability_sRNARFTarget"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {csv_path}: {missing}")
    return df

def filter_outliers(df, f):
    mask = (
        (df["E_intaRNA"] <= f["E_intaRNA_max"]) &
        (df["E_Rnaplex"] <= f["E_Rnaplex_max"]) &
        (df["E_TargetRNA3"] <= f["E_TargetRNA3_max"]) &
        (df["Probability_TargetRNA3"] >= f["Probability_TargetRNA3_min"]) &
        (df["Probability_sRNARFTarget"] >= f["Probability_sRNARFTarget_min"])
    )
    return df.loc[mask].copy()

def load_deg_tables(files, padj):
    out = {}
    for fp in files:
        _check_exists(fp)
        df = pd.read_csv(fp, sep="\t")
        if not {"gene_id","log2FoldChange","padj"}.issubset(df.columns):
            raise ValueError(f"{fp} must contain gene_id, log2FoldChange, padj")
        out[fp] = df.loc[df["padj"]<=padj,["gene_id","log2FoldChange"]].copy()
    return out

def merge_deg_all(d):
    merged=None
    for fp,df in d.items():
        name = fp.replace("biofilm_vs_plank_","").replace(".deseq2.results.tsv","")
        tmp = df.rename(columns={"log2FoldChange":name})
        merged = tmp if merged is None else pd.merge(merged,tmp,on="gene_id",how="outer")
    return merged

def compute_consistent_status(df, min_strains):
    num = df.select_dtypes(include=[np.number])
    def row_status(r):
        if (r>0).sum()>=min_strains: return "upregulated"
        if (r<0).sum()>=min_strains: return "downregulated"
        return np.nan
    out=df.copy()
    out["regulation_status"]=num.apply(row_status,axis=1)
    return out.dropna(subset=["regulation_status"])

def annotate_deg(pred,deg):
    m = deg.set_index("gene_id")["regulation_status"].to_dict()
    pred=pred.copy()
    pred["sRNA_DEG"]=pred["sRNA"].map(m)
    pred["Target_DEG"]=pred["Target"].map(m)
    return pred

def load_modules(fp):
    _check_exists(fp)
    df = pd.read_csv(fp, sep="\t")
    # Fix specific bad column name from R export
    df = df.rename(columns={"nodeAttr[nodesPresent, ]": "nodeAttr"})
    if not {"nodeName","nodeAttr"}.issubset(df.columns):
        raise ValueError("nodes file must contain nodeName,nodeAttr")
    return df

def annotate_modules(pred,mod):
    m=mod.set_index("nodeName")["nodeAttr"].to_dict()
    pred=pred.copy()
    pred["sRNA_Module"]=pred["sRNA"].map(m)
    pred["Target_Module"]=pred["Target"].map(m)
    return pred

def load_edges(fp):
    _check_exists(fp)
    df=pd.read_csv(fp,sep="\t")
    if not {"fromNode","toNode","weight"}.issubset(df.columns):
        raise ValueError("edges file must contain fromNode,toNode,weight")
    df["pair"]=df.apply(lambda r: tuple(sorted([r["fromNode"],r["toNode"]])),axis=1)
    return df[["pair","weight"]]

def pick_best_weight(df,edges):
    tmp=df.copy()
    tmp["pair"]=tmp.apply(lambda r: tuple(sorted([r["sRNA"],r["Target"]])),axis=1)
    m=tmp.merge(edges,on="pair",how="left")
    m=m.sort_values(["Target","weight"],ascending=[True,False])
    m=m[~m["weight"].isna()].drop_duplicates("Target")
    return m.drop(columns=["pair"])

# ## 3) Run pipeline

# 3.1 Load predictions
pred = load_predictions(PREDICOES_CSV)
print("Total predictions:",len(pred))
# pred.head()   # (notebook preview; no-op in script)

# 3.2 Apply energy/probability filters
pred_f = filter_outliers(pred,FILTERS_OUTLIER)
print("After outlier filter:",len(pred_f))
# pred_f.head()   # (notebook preview; no-op in script)

# 3.3 Load DEGs and compute cross‑strain consistency
deg_raw = load_deg_tables(DEG_FILES,PADJ_CUTOFF)
deg_all = merge_deg_all(deg_raw)
deg_cons = compute_consistent_status(deg_all,MIN_STRAINS_CONSISTENT)
print("Consistent DEGs:",len(deg_cons))
# deg_cons.head()   # (notebook preview; no-op in script)

# 3.4 Annotate predictions with DEG status (sRNA + target)
pred_deg = annotate_deg(pred_f,deg_cons)
pred_deg = pred_deg.dropna(subset=["sRNA_DEG","Target_DEG"])
print("After DEG filter:",len(pred_deg))
# pred_deg.head()   # (notebook preview; no-op in script)

# 3.5 Keep only interactions inside same STRING module
mods = load_modules(MODULE_NODES_FILE)
pred_mod = annotate_modules(pred_deg,mods)
filtered = pred_mod[pred_mod["sRNA_Module"]==pred_mod["Target_Module"]]
print("After module filter:",len(filtered))
# filtered.head()   # (notebook preview; no-op in script)

# 3.6 Quick stats
print("Unique sRNAs:",filtered["sRNA"].nunique())
print("Unique targets:",filtered["Target"].nunique())
print("Mean sRNAs per target:",filtered.groupby("Target")["sRNA"].nunique().mean())

# 3.7 Attach edge weights and keep best sRNA per target
edges = load_edges(MODULE_EDGES_FILE)
final_df = pick_best_weight(filtered,edges)
print("Final interactions:",len(final_df))
# final_df.head()   # (notebook preview; no-op in script)

# 3.8 Save output
final_df.to_csv(OUTPUT_CSV,index=False)
print("Saved:",OUTPUT_CSV)

