# Position Categorization Module

## Overview

This module classifies small RNA (sRNA) genomic coordinates relative to annotated gene intervals in a strand-aware manner.

Each sRNA is assigned to one of the following categories:

- **Intragenic**
- **Antisense**
- **5' UTR**
- **3' UTR**
- **Intergenic**

The classification procedure is deterministic and based on configurable overlap thresholds and strand orientation.  
This implementation was developed to support the reproducible analyses described in the associated manuscript.

---

## Methodological Description

The classification is applied hierarchically as follows:

<img width="678" height="615" alt="image" src="https://github.com/user-attachments/assets/dbc7baf3-39d3-4336-bcd1-ba38865d0308" />

### 1. Intragenic / Antisense

An sRNA is classified as overlapping a gene if:

- It is fully contained within the gene interval, **or**
- The overlap between the sRNA and gene is ≥ 90% of the sRNA length.

Strand interpretation:

- Same strand → **Intragenic**
- Opposite strand → **Antisense**

---

### 2. UTR Classification (Strand-Aware)

Unclassified sRNAs are evaluated against upstream and downstream windows defined around each gene (default: ±150 bp).

If the overlap between the sRNA and the UTR window meets the defined threshold:

- For genes on the `+` strand:
  - Upstream → **5' UTR**
  - Downstream → **3' UTR**
- For genes on the `-` strand:
  - Upstream → **3' UTR**
  - Downstream → **5' UTR**

---

### 3. Intergenic

Any sRNA not assigned in previous passes is classified as:

**Intergenic**

---

## Input Data

### sRNA Table

CSV or TSV file containing the following columns:

| Column   | Description                        |
|----------|------------------------------------|
| exon_id  | Unique identifier for the sRNA     |
| start    | Genomic start coordinate           |
| end      | Genomic end coordinate             |
| strand   | Strand (`+` or `-`)                |

---

### Gene Annotation

Either:

- A processed gene interval table with the same columns as above, or  
- A GTF file from which gene intervals are extracted.

---

## Output

The module generates a table (Excel or CSV) containing:

| Column    | Description                              |
|-----------|------------------------------------------|
| exon_id   | sRNA identifier                          |
| start     | Genomic start coordinate                 |
| end       | Genomic end coordinate                   |
| strand    | Strand                                   |
| location  | Assigned category                        |
| gene_id   | Associated gene (if applicable)          |

---

# Requirements

## Software

- Python ≥ 3.9
- Jupyter Notebook (if using the notebook version)

## Python Dependencies

Install using:

```bash
pip install pandas numpy openpyxl
