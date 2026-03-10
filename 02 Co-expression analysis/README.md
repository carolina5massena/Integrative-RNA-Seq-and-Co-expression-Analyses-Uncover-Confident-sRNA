# WGCNA Gene Co-expression Network Analysis

## Overview

This script performs a complete **Weighted Gene Co-expression Network Analysis (WGCNA)** starting from:

1. A gene-level RNA-seq expression matrix  
2. A sample trait table  

The workflow includes:

- Sample clustering and outlier detection  
- Soft-threshold power selection  
- Network construction (Adjacency and TOM)  
- Dynamic module detection and merging  
- Module–trait correlation analysis  
- Hub gene identification  
- Cytoscape network export  

All user-editable parameters are centralized in the **USER INPUTS** section at the top of the script.

---

# Inputs

The script requires **two mandatory input files**.

## 1) Expression Matrix (TSV)

Defined in the script as:

```r
expr_counts_tsv <- "salmon.merged.gene_counts.tsv"
```

### Required Format

- Tab-separated file (`.tsv`)
- First column: gene identifier
- Remaining columns: numeric sample expression values
- No missing values allowed

### Required Columns

| Column | Description |
|--------|------------|
| `gene_id_col` (default: `gene_name`) | Gene identifier column |
| Sample columns | Numeric expression values |

### Example

```
gene_name	Sample1	Sample2	Sample3
GeneA	    10	    15	    8
GeneB	    5	    7	    6
```

### Important Requirements

- All sample columns must be numeric.
- Duplicate gene IDs are automatically aggregated (summed).
- Sample column names must exactly match the sample IDs in the trait table.

---

## 2) Trait Table (CSV)

Defined in the script as:

```r
traits_csv <- "table_treatment.csv"
```

### Required Format

- Comma-separated file (`.csv`)
- Must contain sample identifiers and at least one categorical trait column.

### Required Columns

| Parameter | Default | Description |
|-----------|----------|------------|
| `trait_sample_col` | `rownames` | Sample ID column |
| `treatment_col` | `treatment` | Categorical variable for module–trait correlation |

### Example

```
rownames,treatment
Sample1,Control
Sample2,Treated
Sample3,Treated
```

### Critical Matching Rule

Sample IDs in the trait file **must exactly match** the expression matrix sample column names.

---

# Outputs

All outputs are written to the directory defined by:

```r
out_dir <- "."
```

The script generates the following files:

## Quality Control

- `sample_dendrogram_outlier_check.png`
- `sample_dendrogram_with_traits.png`

## Network Construction

- `soft_thresholding_power_selection.png`
- `gene_dendrogram_with_modules.png`

## Module–Trait Analysis

- `Module_Trait_Heatmap.pdf`

## Gene and Network Exports

- `gene_module_membership.csv`
- `top_hub_genes_per_module.csv`
- `CytoscapeInput-edges.txt`
- `CytoscapeInput-nodes.txt`

---

# Required Parameters

Located in the `USER INPUTS` section of the script.

## Input File Paths

```r
expr_counts_tsv <- "your_expression_file.tsv"
traits_csv      <- "your_trait_file.csv"
```

## Column Names (change only if necessary)

```r
gene_id_col      <- "gene_name"
trait_sample_col <- "rownames"
treatment_col    <- "treatment"
```

## Soft-Threshold Power (MANDATORY)

```r
softPower <- 10
```

To select `softPower`:

1. Run the script once.
2. Inspect `soft_thresholding_power_selection.png`.
3. Choose a power where the Scale-Free Topology Fit (R²) is acceptable (commonly ≥ 0.8).
4. Update `softPower`.
5. Re-run the analysis.

---

# Requirements

## Software

- R ≥ 4.1

## Required R Packages

Install before running:

```r
install.packages(c("vctrs", "tibble", "fastDummies", "flashClust"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("WGCNA", ask = FALSE, update = FALSE)
```

The script internally enables multithreading using `allowWGCNAThreads()`.

---

# Execution

Run from terminal:

```bash
Rscript WGCNA.R
```

Or inside R:

```r
source("WGCNA.R")
```
