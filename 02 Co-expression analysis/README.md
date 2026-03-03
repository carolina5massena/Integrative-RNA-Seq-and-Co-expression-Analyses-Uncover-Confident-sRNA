# WGCNA Network Analysis Module

## Overview

This module performs a complete **Weighted Gene Co-expression Network Analysis (WGCNA)** workflow starting from a gene-level expression matrix and a sample trait table.

The pipeline includes:

- Sample clustering and outlier detection
- Soft-threshold power selection
- Network construction and TOM calculation
- Dynamic module detection
- Module merging
- Module–trait correlation analysis
- Hub gene identification
- Cytoscape network export

All parameters are centralized in the **USER INPUTS** section at the top of the script to ensure transparency and reproducibility.

---

# Requirements

## Software

- **R ≥ 4.1**

## Required R Packages

Install before running:

```r
install.packages(c(
  "tidyverse",
  "WGCNA",
  "flashClust",
  "fastDummies"
))
