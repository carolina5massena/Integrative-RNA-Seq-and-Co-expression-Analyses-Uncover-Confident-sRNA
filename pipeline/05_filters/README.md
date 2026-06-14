# sRNA–mRNA Target Filtering and Integration Pipeline

## Overview

The pipeline produces a **high-confidence, module-consistent, cross-strain-supported interaction table**.

The workflow includes:

1. Energy and probability filtering of predicted interactions  
2. Cross-strain DEG consistency filtering  
3. Annotation of sRNA and target DEG status  
4. Network weight integration  
5. Selection of the strongest sRNA per target  

All configurable parameters are centralized in the **USER CONFIGURATION** section at the top of the script.

This module ships in two equivalent forms:

- **`Filters_sRNA_mRNAtarget.py`** — a standalone script (edit the config block, then run).
- **`sRNA_mRNAtarget.ipynb`** — the same logic as a Jupyter notebook for interactive use.

The same step is also available, fully configurable, as the **05 Filters** tab in the GUI.

---

# Inputs

The script requires the following input files.

---

## 1) Interaction Prediction File (CSV)

Defined in the script:

```python
PREDICOES_CSV = "Prediction_test.csv"
```

### Required Columns

| Column | Description |
|--------|------------|
| `sRNA` | sRNA identifier |
| `Target` | Target gene identifier |
| `E_intaRNA` | IntaRNA interaction energy |
| `E_Rnaplex` | RNAplex interaction energy |
| `E_TargetRNA3` | TargetRNA3 interaction energy |
| `Probability_TargetRNA3` | TargetRNA3 probability score |
| `Probability_sRNARFTarget` | sRNARFTarget probability score |

These columns are mandatory.

---

## 2) DESeq2 Result Files (TSV)

Defined in the script:

```python
DEG_FILES = [
    "DEG_result_A.tsv",
    "DEG_result_B.tsv",
    # ... one entry per comparison
]
```

Each file must contain:

| Column | Description |
|--------|------------|
| `gene_id` | Gene identifier |
| `log2FoldChange` | Expression fold change |
| `padj` | Adjusted p-value |

Filtering applied:

- `padj <= PADJ_CUTOFF`
- Cross-strain consistency:  
  At least `MIN_STRAINS_CONSISTENT` strains must agree on direction  
  (>0 = upregulated, <0 = downregulated)

---

## 3) WGCNA nodes Module Nodes File (TSV)

Defined in the script:

```python
MODULE_NODES_FILE = "CytoscapeInput-nodes.txt"
```

Required columns:

| Column | Description |
|--------|------------|
| `nodeName` | Gene name |
| `nodeAttr` | Module identifier |

---

## 4) WGCNA edges Module Edges File (TSV)

Defined in the script:

```python
MODULE_EDGES_FILE = "CytoscapeInput-edges.txt"
```

Required columns:

| Column | Description |
|--------|------------|
| `fromNode` | Gene A |
| `toNode` | Gene B |
| `weight` | Edge weight |

---

# Filtering Logic

## Step 1 — Energy and Probability Filters

Thresholds defined in:

```python
FILTERS_OUTLIER = {
    "E_intaRNA_max": -2.44,
    "E_Rnaplex_max": -32.6,
    "E_TargetRNA3_max": -5.13,
    "Probability_TargetRNA3_min": 0.06,
    "Probability_sRNARFTarget_min": 0.40
}
```

Only interactions passing all thresholds are retained.

---

## Step 2 — Cross-Strain DEG Consistency

Parameters:

```python
PADJ_CUTOFF = 0.05
MIN_STRAINS_CONSISTENT = 2
```

A gene is considered:

- **Upregulated** if ≥ MIN_STRAINS_CONSISTENT strains have log2FC > 0  
- **Downregulated** if ≥ MIN_STRAINS_CONSISTENT strains have log2FC < 0  

Both sRNA and target must have consistent DEG status.

---

## Step 3 — WGCNA Module Consistency

Only interactions where:

```
sRNA_Module == Target_Module
```

are retained.

---

## Step 4 — Network Weight Selection

- For each target, only the sRNA with the highest weight is retained.

---

# Output

Defined in the script:

```python
OUTPUT_CSV = "filtered_weight.csv"
```

The final output contains high-confidence interactions including:

- sRNA  
- Target  
- Energy scores  
- Probability scores  
- DEG status  
- Module annotation  

Each target appears only once (best-weight interaction retained).

---

# Requirements

## Software

- Python ≥ 3.9

## Python Packages

Install once:

```bash
pip install pandas numpy
```

---

# Execution

Edit the **USER CONFIGURATION** block (input paths and parameters), then run
from the terminal:

```bash
python Filters_sRNA_mRNAtarget.py
```

Or open `sRNA_mRNAtarget.ipynb` and run the cells in order. Either way, make sure
the input files referenced in the config (prediction CSV, DEG TSVs, and the
Cytoscape nodes/edges files) are reachable from the paths you set.
