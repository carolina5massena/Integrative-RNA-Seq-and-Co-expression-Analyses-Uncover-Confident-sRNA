# Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA

A reproducible pipeline that identifies and prioritizes biologically consistent
**sRNA–mRNA regulatory interactions** by combining sequence-based target
prediction with transcriptomic and network-level evidence. A desktop GUI wraps
the analysis steps so the full workflow can be run without writing code.

> **New here?** Jump straight to the [**Quick Start**](#quick-start) to install,
> then follow the [**Tutorial — Running the First Test**](#tutorial--running-the-first-test)
> to run the whole pipeline end-to-end on the bundled example data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Required External Tools](#required-external-tools)
- [Software Requirements](#software-requirements)
- [Quick Start](#quick-start)
- [The GUI](#the-gui)
- [Tutorial — Running the First Test](#tutorial--running-the-first-test)
- [Pipeline Modules (reference)](#pipeline-modules-reference)
- [Auxiliary RNA-Seq Commands](#auxiliary-rna-seq-commands)
- [Troubleshooting](#troubleshooting)
- [Contact & Citation](#contact--citation)

---

## Overview

The pipeline integrates several lines of evidence to move from raw sequencing
data to a short list of high-confidence regulatory interactions:

1. **Position Classification** — classify each sRNA relative to annotated genes.
2. **RNA-seq processing** — quantify expression with `nf-core/rnaseq`
   *(see [Auxiliary Commands](#auxiliary-rna-seq-commands))*.
3. **Co-expression analysis (WGCNA)** — build a gene network and detect modules.
4. **Differential expression (DESeq2)** — find differentially expressed genes
   *(also available via [Auxiliary Commands](#auxiliary-rna-seq-commands) using `nf-core/differentialabundance`)*.
5. **Target prediction** — run external sRNA–mRNA prediction tools.
6. **Filtering & integration** — combine all evidence into a final interaction table.

By requiring agreement across sequence prediction, differential expression, and
network module membership, the workflow yields interactions supported by multiple
independent signals.

<img width="606" height="340" alt="Pipeline overview" src="https://github.com/user-attachments/assets/8f1a853f-b106-4f29-a660-7e94f5b9712e" />

*Processes drawn with a dotted outline in the figure are not included in this repository.*

---

## Repository Structure

```
.
├── README.md                      ← you are here (start point + tutorial)
├── INSTALL.md                     ← detailed install & troubleshooting
├── requirements.txt               ← Python dependencies
├── install_packages.R             ← R/Bioconductor dependencies
├── setup.bat / setup.sh           ← one-click environment setup
├── run.bat   / run.sh             ← launch the GUI
├── AUXILIARY_COMMANDS.sh          ← SRA download + nf-core RNA-seq / DESeq2 commands
│
├── GUI/                           ← PySide6 desktop application
│   ├── app.py                     ← the interface
│   ├── pipeline_logic.py          ← step implementations
│   └── README.md                  ← GUI details
│
├── pipeline/                      ← standalone (command-line) analysis modules
│   ├── 01_position_classification/   ← sRNA position categorization (Python)
│   ├── 02_coexpression_wgcna/        ← WGCNA network analysis (R)
│   └── 05_filters/                   ← evidence integration & filtering (Python)
│
├── data/                          ← example inputs for the tutorial
├── Results/                       ← example outputs produced by the tutorial
└── docs/                          ← tutorial as PDF / DOCX (archived copies)
```

Each module under `pipeline/` contains its own `README.md` with a detailed
description of inputs, outputs, and parameters — see
[Pipeline Modules](#pipeline-modules-reference) for links. The `pipeline/`
modules and the GUI are two ways to run the same analyses: use the GUI for a
guided, no-code workflow, or the standalone scripts for command-line / scripted
runs.

---

## Required External Tools

Target prediction is performed by external programs, which must be installed
separately. Run any subset; the [Filters](#pipeline-modules-reference) step
combines whichever outputs are available.

| Tool | Purpose | Link |
|------|---------|------|
| **IntaRNA** (`IntaRNAsTar` personality) | sRNA–mRNA interaction prediction, optimized for genome-wide sRNA-target screens | https://github.com/BackofenLab/IntaRNA |
| **RNAplex** (ViennaRNA Package) | Fast interaction prediction | https://www.tbi.univie.ac.at/RNA/ |
| **TargetRNA3** | Target prediction with probability score | https://cs.wellesley.edu/~btjaden/TargetRNA3 |
| **sRNARFTarget** | Machine-learning target prediction | https://github.com/BioinformaticsLabAtMUN/sRNARFTarget |

---

## Software Requirements

The pipeline uses two runtimes:

- **Python ≥ 3.10** — the GUI and the Python analysis steps.
- **R ≥ 4.1** — the WGCNA and DESeq2 analyses.

Install those two once per machine. The setup scripts below handle every Python
and R package and locate `Rscript` automatically. See [`INSTALL.md`](INSTALL.md)
for full details and troubleshooting.

---

## Quick Start

1. Install [Python](https://www.python.org/downloads/) (on Windows, tick
   *"Add Python to PATH"*) and [R](https://cran.r-project.org/). You do **not**
   need to add R to PATH — setup finds it for you.

2. Run the one-click setup:

   **Windows** — double-click `setup.bat`.
   **macOS / Linux** — `bash setup.sh`

3. Launch the application:

   **Windows** — double-click `run.bat`.
   **macOS / Linux** — `bash run.sh`

Setup creates a local `.venv`, installs the Python packages from
`requirements.txt`, locates `Rscript`, and installs the R packages via
`install_packages.R`. It is safe to re-run.

<details>
<summary><b>Manual install (if you prefer not to use the scripts)</b></summary>

```bash
# Python (PySide6, pandas, numpy, openpyxl)
pip install -r requirements.txt

# R: CRAN (tibble, fastDummies, flashClust) + Bioconductor (WGCNA, DESeq2)
Rscript install_packages.R

# Launch
python GUI/app.py
```
</details>

---

## The GUI

A PySide6 desktop application (`GUI/app.py`) wraps the analysis steps behind a
tabbed interface. Every input is editable in the UI, each tab is self-contained,
and the path to `Rscript` is detected automatically and pre-filled.

| Tab | What it does | Runtime |
|-----|--------------|---------|
| **01 Position Classification** | Classifies sRNAs from a GTF file using overlap threshold and UTR window. | Python |
| **02 Co-expression analysis** | Configures and runs WGCNA; generates `WGCNA_configured.R`. | R |
| **03 Differential Expression** | Runs DESeq2 on a count matrix + metadata; generates `DESeq2_configured.R` and one DEG file per comparison. | R |
| **04 Combine Predictions** | Merges the four target-prediction outputs into one `Prediction_combined.csv`. | Python |
| **05 Filters** | Integrates predictions, DEGs, and WGCNA modules into the final table. | Python |

Each analysis tab has its own **Run** button; work runs in a background thread
and streams a log. See [`GUI/README.md`](GUI/README.md) for full details.

---

## Tutorial — Running the First Test

This walkthrough runs the entire pipeline on the bundled example data in
`data/`. The expected console output is shown for each step so you can
confirm your installation is working. *(A PDF/DOCX copy lives in [`docs/`](docs).)*

### A. Installation and setup

1. Run `setup.bat` (or `setup.sh`) to install the dependencies into a local
   virtual environment.
2. Once installation completes, run `run.bat` (or `run.sh`) to launch the
   application.

### 1. Position Classification

1. Select the GTF annotation file `Test_annotation.gtf` from the `data` folder.
3. Click **Run Position Classification**.

**Expected output:**

```
Found 681 sRNAs and 3013 genes.
Category counts:
location
Intergenic    374
Antisense     219
3' UTR         31
5' UTR         30
Intragenic     27
Saved sRNA_annotation.xlsx in the desired output directory.
```

### 2. Co-expression Analysis (WGCNA)

1. Select `Count_test.tsv` as the expression matrix (TSV).
2. Select `table_treatment_WGCNA.csv` as the trait table (CSV).
4. Adjust the Rscript executable path if required.
5. Click **Run WGCNA**.

The following files are generated in the results folder:

```
sample_dendrogram_outlier_check.png
sample_dendrogram_with_traits.png
soft_thresholding_power_selection.png
gene_dendrogram_with_modules.png
Module_Trait_Heatmap.pdf          ← module–trait heatmap for comparison
gene_module_membership.csv
top_hub_genes_per_module.csv
CytoscapeInput-edges.txt
CytoscapeInput-nodes.txt
```

### 3. Differential Expression (DESeq2)

1. Set the count matrix TSV to `Count_test.tsv`.
2. Set the sample metadata CSV to `metadata_DEG.csv`.
3. For this test, set the **Comparison** column to `comparison`. This splits the
   metadata into groups (e.g., A, B, C, D) and runs DESeq2 once per group,
   producing one DEG file per comparison (`DEG_result_A.tsv`,
   `DEG_result_B.tsv`, …).
5. Adjust the Rscript executable path if required.
6. Click **Run DESeq2**.

**Expected output:**

```
Shared dispersion estimated from 8 samples over 3211 genes.
comparison A -> DEG_result_A.tsv [own dispersion] sig padj<0.05: 591
comparison B -> DEG_result_B.tsv [own dispersion] sig padj<0.05: 1804
```

### 4. Combine Predictions

The four target-prediction programs each produce a separate output file (see
[Auxiliary RNA-Seq Commands](#auxiliary-rna-seq-commands) for how to run them).
The **04 Combine Predictions** tab merges those files, on each (Target, sRNA)
pair, into the single prediction table that the Filters step consumes.

1. Open the **04 Combine Predictions** tab.
2. Add each program's output file from the `data/` folder:
   - **intaRNA file:** `intaRNA_output.txt`
   - **RNAplex file:** `Rnaplex_output.txt`
   - **TargetRNA3 file:** `TargetRNA3_output.txt`
   - **sRNARFTarget file:** `sRNARFTarget_output.csv`
3. Leave the default options (compute empirical p-values for intaRNA/RNAplex;
   keep only pairs predicted by every program provided).
4. Set the **Output CSV** to `Prediction_combined.csv`.
5. Click **Combine Predictions**.

**Expected output:**

```
intaRNA: 32789 unique target-sRNA pairs
RNAplex: 32789 unique target-sRNA pairs
TargetRNA3: 32789 unique target-sRNA pairs
sRNARFTarget: 32789 unique target-sRNA pairs
Merged (inner join): 32789 target-sRNA pairs
Combined prediction table: 32789 rows, 10 columns.
```

The resulting `Prediction_combined.csv` has the same columns as the bundled
`Prediction_test.csv` and is used as the prediction input in the next step.

### 5. Filtering

1. **Prediction CSV:** select `Prediction_combined.csv` (from step 4) — or the
   bundled `Prediction_test.csv`.
2. **DEG TSV files:** add both `DEG_result_A.tsv` and `DEG_result_B.tsv`.
3. **Cytoscape Node:** add `CytoscapeInput-nodes.txt`.
4. **Cytoscape Edges:** add `CytoscapeInput-edges.txt`.
5. The **Min strains consistent** value defaults to `2` — leave it as is for
   this test.

**Expected output (using `Prediction_combined.csv`):**

```
Total predictions: 32789
After energy/probability filter (>= 5/5 of
  ['E_intaRNA', 'E_Rnaplex', 'E_TargetRNA3',
   'Probability_TargetRNA3', 'Probability_sRNARFTarget']): 28586
Consistent DEGs: 572
After DEG filter: 421
After module filter (same): 233
After requiring network edge: 182
After best_per_target selection: 146
Final interactions: 146
```

*(With the bundled `Prediction_test.csv` the totals differ slightly — 33996
predictions — but the final count is also 146.)*

#### Optional — sRNA location filter

You can restrict the results to sRNAs of a chosen genomic category using the
annotation table from step 1. In the **sRNA location filter** group of the
Filters tab:

1. Tick **Enable location filter**.
2. **Annotation file:** add `sRNA_annotation.xlsx` (produced in step 1).
3. Under **Keep locations**, leave only **Antisense** ticked.

**Expected output (Antisense only):**

```
Total predictions: 32789
After energy/probability filter (>= 5/5 of ['E_intaRNA', 'E_Rnaplex', 'E_TargetRNA3', 'Probability_TargetRNA3', 'Probability_sRNARFTarget']): 28586
Consistent DEGs: 572
After DEG filter: 421
After location filter (Antisense): 101 (from 421; sRNAs not found in the annotation file are dropped)
After module filter (same): 58
After requiring network edge: 48
After best_per_target selection: 46
Final interactions: 46
```

If your numbers match, the pipeline is installed correctly and you are ready to
run it on your own data.

---

## Pipeline Modules (reference)

Detailed documentation for each step lives in its module folder:

| Module | Description | Docs |
|--------|-------------|------|
| **01 Position Classification** | Strand-aware classification of sRNAs as Intragenic, Antisense, 5′ UTR, 3′ UTR, or Intergenic. | [`pipeline/01_position_classification/`](pipeline/01_position_classification/README.md) |
| **02 Co-expression (WGCNA)** | Full WGCNA workflow: clustering, soft-threshold selection, module detection, module–trait correlation, hub genes, Cytoscape export. | [`pipeline/02_coexpression_wgcna/`](pipeline/02_coexpression_wgcna/README.md) |
| **05 Filters** | Energy/probability filtering, cross-strain DEG consistency, module consistency, and best-per-target selection. | [`pipeline/05_filters/`](pipeline/05_filters/README.md) |
| **GUI** | The PySide6 application that drives all of the above. | [`GUI/README.md`](GUI/README.md) |

---

## Auxiliary RNA-Seq Commands

The RNA-seq quantification and the sRNA–mRNA target-prediction programs are run
outside the GUI. Reference commands for all of them are collected in
[`AUXILIARY_COMMANDS.sh`](AUXILIARY_COMMANDS.sh):

- Download raw reads from SRA (`prefetch`)
- Convert SRA to FASTQ (`fasterq-dump`)
- Run `nf-core/rnaseq`
- Run the four target-prediction tools — **IntaRNA**, **RNAplex**,
  **TargetRNA3**, and **sRNARFTarget** — whose outputs are merged into the
  prediction table consumed by [05 Filters](#pipeline-modules-reference)

**Requirements:** SRA Toolkit, Nextflow, Docker, the nf-core pipelines, and the
four prediction tools (see [Required External Tools](#required-external-tools)).

---

## Troubleshooting

**`OSError: [Errno 2]` / "enable long paths" during `pip install` (Windows).**
This is Windows' 260-character path limit being hit while installing PySide6
(deeply nested files), made worse by a long project path. Fix it by enabling
long-path support:

1. Open **PowerShell as Administrator** and run:
   ```powershell
   New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" -Name "LongPathsEnabled" -Value 1 -PropertyType DWORD -Force
   ```
2. **Restart the computer** for the setting to take effect.
3. Delete the partially-created `.venv` folder, then run `setup.bat` again.

   Alternatively, move the project to a short path such as `C:\srna` before
   running setup.

**`FileNotFoundError: [WinError 2]` when running an analysis.** R isn't being
found. Re-run setup, or paste the full path to `Rscript.exe` into the GUI field.

**`there is no package called 'WGCNA'` / `'DESeq2'`.** The R packages didn't
install. Run `Rscript install_packages.R` (Bioconductor installs can take several
minutes the first time).

More cases are covered in [`INSTALL.md`](INSTALL.md#troubleshooting).

---

## Contact & Citation

For issues or reproducibility questions, please open a GitHub issue in this
repository. If you use this pipeline in your work, please cite the associated
manuscript.
