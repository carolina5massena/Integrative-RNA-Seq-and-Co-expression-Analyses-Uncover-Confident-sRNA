# Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA

A reproducible pipeline that identifies and prioritizes biologically consistent
**sRNA–mRNA regulatory interactions** by combining sequence-based target
prediction with transcriptomic and network-level evidence.

Everything runs through a single **desktop GUI** — no scripting required. The GUI
is organized as five tabs that take you from a genome annotation to a short list
of high-confidence interactions, and this README documents every tab and every
option in one place.

> **New here?** Run the [**Quick Start**](#quick-start) to install and launch the
> app, then follow the [**Tutorial — Running the First Test**](#tutorial--running-the-first-test)
> to reproduce a full run on the bundled example data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Required External Tools](#required-external-tools)
- [Software Requirements](#software-requirements)
- [Quick Start](#quick-start)
- [Using the GUI](#using-the-gui)
  - [01 Position Classification](#01-position-classification)
  - [02 Co-expression Analysis (WGCNA)](#02-co-expression-analysis-wgcna)
  - [03 Differential Expression (DESeq2)](#03-differential-expression-deseq2)
  - [04 Combine Predictions](#04-combine-predictions)
  - [05 Filters](#05-filters)
- [Tutorial — Running the First Test](#tutorial--running-the-first-test)
- [Auxiliary RNA-Seq Commands](#auxiliary-rna-seq-commands)
- [Troubleshooting](#troubleshooting)
- [Contact & Citation](#contact--citation)

---

## Overview

The pipeline integrates several lines of evidence to move from raw sequencing
data to a short list of high-confidence regulatory interactions:

1. **Position Classification** — classify each sRNA relative to annotated genes.
2. **RNA-seq processing** — quantify expression with `nf-core/rnaseq`
   *(run outside the GUI; see [Auxiliary Commands](#auxiliary-rna-seq-commands))*.
3. **Co-expression analysis (WGCNA)** — build a gene network and detect modules.
4. **Differential expression (DESeq2)** — find differentially expressed genes.
5. **Target prediction** — run external sRNA–mRNA prediction tools
   *(run outside the GUI; see [Auxiliary Commands](#auxiliary-rna-seq-commands))*.
6. **Combine & filter** — merge the predictions and filter them against the
   expression and network evidence into a final interaction table.

By requiring agreement across sequence prediction, differential expression, and
network module membership, the workflow yields interactions supported by multiple
independent signals.

The GUI covers steps 1, 3, 4 and 6 directly, and configures + runs the WGCNA
analysis for step 3. Steps that depend on heavy external tools — RNA-seq
quantification and the four target-prediction programs — are run from the command
line using the reference commands in
[`AUXILIARY_COMMANDS.sh`](AUXILIARY_COMMANDS.sh), and their outputs are then
loaded back into the GUI.

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/20f6f1f8-3c42-4e23-9db3-f3debeb7b761" />

*Processes drawn with a dotted outline in the figure are not included in this repository.*

---

## Repository Structure

```
.
├── README.md                      ← you are here (full guide + tutorial)
├── INSTALL.md                     ← detailed install & troubleshooting
├── requirements.txt               ← Python dependencies
├── install_packages.R             ← R/Bioconductor dependencies
├── setup.bat / setup.sh           ← one-click environment setup
├── run.bat   / run.sh             ← launch the GUI
├── AUXILIARY_COMMANDS.sh          ← SRA download, nf-core RNA-seq, target prediction
│
├── GUI/                           ← the desktop application (everything runs here)
│   ├── app.py                     ← the PySide6 interface (the five tabs)
│   └── pipeline_logic.py          ← the analysis implementations
│
├── data/                          ← example inputs for the tutorial
├── Results/                       ← example outputs produced by the tutorial
└── docs/                          ← tutorial as PDF / DOCX (archived copies)
```

The entire workflow is driven by the GUI in `GUI/app.py`; `pipeline_logic.py`
holds the analysis code and the embedded R templates it generates and runs. There
are no separate command-line modules to maintain — every step and option lives in
the GUI and is documented in [Using the GUI](#using-the-gui).

---

## Required External Tools

Target prediction is performed by external programs, installed separately. Run
any subset; the [Combine Predictions](#04-combine-predictions) step merges
whichever outputs are available. See
[Auxiliary RNA-Seq Commands](#auxiliary-rna-seq-commands) for the exact commands.

| Tool | Purpose | Link |
|------|---------|------|
| **IntaRNA** (`IntaRNAsTar` personality) | sRNA–mRNA interaction prediction, optimized for genome-wide sRNA-target screens | https://github.com/BackofenLab/IntaRNA |
| **RNAplex** (ViennaRNA Package) | Fast interaction prediction | https://www.tbi.univie.ac.at/RNA/ |
| **TargetRNA3** | Target prediction with probability score | https://cs.wellesley.edu/~btjaden/TargetRNA3 |
| **sRNARFTarget** | Machine-learning target prediction | https://github.com/BioinformaticsLabAtMUN/sRNARFTarget |

RNA-seq quantification additionally uses the SRA Toolkit, Nextflow, Docker, and
the `nf-core/rnaseq` pipeline.

---

## Software Requirements

The pipeline uses two runtimes:

- **Python ≥ 3.10** — the GUI and the Python analysis steps (Position
  Classification, Combine Predictions, Filters).
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

## Using the GUI

The application is a tabbed interface. Every input is editable in the UI, each
tab is self-contained, and the path to `Rscript` is detected automatically and
pre-filled where R is needed. Each analysis tab has its own **Run** button; work
runs in a background thread and streams a log, and a dialog reports
success/failure.

| Tab | What it does | Runtime |
|-----|--------------|---------|
| [**01 Position Classification**](#01-position-classification) | Classifies each sRNA relative to annotated genes. | Python |
| [**02 Co-expression Analysis**](#02-co-expression-analysis-wgcna) | Configures and runs WGCNA; builds the co-expression network. | R |
| [**03 Differential Expression**](#03-differential-expression-deseq2) | Runs DESeq2 on a count matrix + metadata. | R |
| [**04 Combine Predictions**](#04-combine-predictions) | Merges the four target-prediction outputs into one table. | Python |
| [**05 Filters**](#05-filters) | Integrates predictions, DEGs, location, and network evidence into the final table. | Python |

A typical run goes top to bottom: classify sRNAs (01), build the network (02),
call DEGs (03), assemble the prediction table (04), then filter everything into
the final list (05).

---

### 01 Position Classification

Classifies each small RNA's genomic coordinates relative to annotated gene
intervals, in a **strand-aware** manner, and assigns one category per sRNA:
**Intragenic**, **Antisense**, **5′ UTR**, **3′ UTR**, or **Intergenic**.

**How classification works (applied hierarchically):**

1. **Intragenic / Antisense** — an sRNA overlaps a gene if it is fully contained
   in the gene interval, or the overlap is ≥ the *overlap threshold* (default 90%
   of the sRNA length). Same strand ⇒ **Intragenic**; opposite strand ⇒
   **Antisense**.
2. **UTR (strand-aware)** — otherwise the sRNA is tested against the windows
   around each gene (default ±150 bp). For `+`-strand genes: upstream ⇒ **5′ UTR**,
   downstream ⇒ **3′ UTR**; for `-`-strand genes the assignment is reversed.
3. **Intergenic** — anything not assigned above.

**Inputs & options**

| Field | Description |
|-------|-------------|
| GTF annotation file | Genome annotation; gene intervals are extracted from it. |
| Feature / ID keys | Which GTF feature and attribute identify genes. |
| sRNA & gene prefixes | Prefixes that distinguish sRNA vs gene records. |
| Rename prefixes | Optionally relabel IDs on output. |
| Overlap threshold | Minimum overlap (fraction of sRNA length) for Intragenic/Antisense (default 0.90). |
| UTR window | Size of the up/downstream UTR window in bp (default 150). |
| Output file | The annotation table to write (default `sRNA_annotation.xlsx`). |

**Output** — a table (`sRNA_annotation.xlsx` by default) with `exon_id`, `start`,
`end`, `strand`, `location` (the assigned category), and `gene_id`. This file
feeds the optional **sRNA location filter** in [05 Filters](#05-filters).

---

### 02 Co-expression Analysis (WGCNA)

Runs a complete **Weighted Gene Co-expression Network Analysis**: sample
clustering and outlier detection, soft-threshold power selection, network
construction (adjacency + TOM), dynamic module detection and merging,
module–trait correlation, hub-gene identification, and a Cytoscape export. The
tab assembles a configured `WGCNA_configured.R` and runs it via `Rscript`.

**Inputs**

- **Expression matrix (TSV)** — first column a gene identifier, remaining columns
  numeric sample values; no missing values. Duplicate gene IDs are summed. Sample
  column names must match the trait table exactly.
- **Trait table (CSV)** — a sample-ID column plus at least one categorical trait
  used for module–trait correlation.

**Key options**

| Field | Description |
|-------|-------------|
| Gene ID / sample / treatment columns | Column names in the two inputs. |
| Soft power | The soft-thresholding power (see below). |
| Min module size | Minimum genes per module. |
| Deep split | Module-splitting sensitivity. |
| Merge cut height | Threshold for merging similar modules. |
| Cytoscape threshold | Edge-weight threshold for the network export. |
| Threads | Number of CPU threads. |
| Network / correlation type | Network and correlation options. |

**Choosing the soft power:** run once, open
`soft_thresholding_power_selection.png`, pick a power where the scale-free
topology fit (R²) is acceptable (commonly ≥ 0.8), set it, and re-run.

**Outputs** (written to the chosen results folder):

```
sample_dendrogram_outlier_check.png
sample_dendrogram_with_traits.png
soft_thresholding_power_selection.png
gene_dendrogram_with_modules.png
Module_Trait_Heatmap.pdf
gene_module_membership.csv
top_hub_genes_per_module.csv
CytoscapeInput-edges.txt        ← fromNode, toNode, weight  (used by 05 Filters)
CytoscapeInput-nodes.txt        ← nodeName, nodeAttr (module)  (used by 05 Filters)
```

---

### 03 Differential Expression (DESeq2)

Runs DESeq2 on a raw count matrix plus sample metadata, producing the DEG tables
used for the cross-strain consistency filter. The tab generates a configured
`DESeq2_configured.R` and runs it via `Rscript`.

**Inputs & options**

| Field | Description |
|-------|-------------|
| Count matrix (TSV) | Raw gene-level counts. |
| Sample metadata (CSV) | Sample annotations. |
| Gene-ID / sample / condition columns | Column names in the inputs. |
| Reference vs. test level | Which condition levels to contrast. |
| Minimum-count filter | Drop low-count genes before testing. |
| **Comparison column** | If set, splits the metadata into groups and runs DESeq2 **once per group**, writing one DEG file per comparison (e.g. `DEG_result_A.tsv`, `DEG_result_B.tsv`, …). |

**Output** — one TSV per comparison containing `gene_id`, `log2FoldChange`,
`padj` (plus `baseMean`, `pvalue`, `lfcSE`, `stat`). These plug straight into the
DEG inputs of [05 Filters](#05-filters).

> RNA-seq quantification (producing the count matrix) and an alternative
> differential-abundance route are run outside the GUI — see
> [Auxiliary RNA-Seq Commands](#auxiliary-rna-seq-commands).

---

### 04 Combine Predictions

Merges the native output of up to four target-prediction programs into a single
prediction table (`Prediction_combined.csv`) with the columns the Filters step
expects. Pairs are merged on each `(Target, sRNA)` combination; every file is
optional — leave a program blank to skip it.

**Inputs**

| Field | Expected format |
|-------|-----------------|
| intaRNA file | `;`-separated: `id1;id2;start1;end1;start2;end2;E` |
| RNAplex file | blocks of `>target`, `>sRNA`, structure line with `(energy)` |
| TargetRNA3 file | tab-separated; must include an `sRNA` column plus `Target`, `Energy`, `P-value`, `Probability` |
| sRNARFTarget file | `sRNA_ID`, `mRNA_ID`, `Prediction_Probability` |

**Options**

| Option | Default | Meaning |
|--------|---------|---------|
| Compute empirical p-values for intaRNA / RNAplex | on | These programs report only an energy; their p-value column is derived as an empirical p-value from the energy distribution. |
| Keep only pairs predicted by every provided program | on | Inner join across the files you supply (untick for an outer join with blanks). |
| Relabel IDs to Gene# / sRNA# | off | Replace names with tidy `Gene1…` / `sRNA1…` labels and write a mapping file. |
| Output CSV | `Prediction_combined.csv` | The merged prediction table. |

**Output** — `Prediction_combined.csv` (same columns as the bundled
`Prediction_test.csv`): `Target`, `sRNA`, `E_intaRNA`, `p_intaRNA`, `E_Rnaplex`,
`p_Rnaplex`, `E_TargetRNA3`, `p_TargetRNA3`, `Probability_TargetRNA3`,
`Probability_sRNARFTarget`.

---

### 05 Filters

The final step. It takes the prediction table and progressively filters it
against differential expression, genomic location, and co-expression-network
evidence. Filtering happens in the order below; **every stage is optional and
configurable**. Defaults shown are the GUI defaults.

#### Inputs

| Input | Required? | Format | Key columns |
|-------|-----------|--------|-------------|
| Prediction table | Optional* | CSV | `sRNA`, `Target` (+ any score columns) |
| DEG result file(s) | Optional | TSV | `gene_id`, `log2FoldChange`, `padj` (+ `baseMean`) |
| WGCNA nodes file | Optional | TSV | `nodeName`, `nodeAttr` |
| WGCNA edges file | Optional* | TSV | `fromNode`, `toNode`, `weight` |
| sRNA annotation (location) | Optional | XLSX/CSV/TSV | sRNA-ID column (e.g. `exon_id`), `location` |

\* Provide **either** a prediction table **or** an edges file. Add multiple DEG
files (comma-separated) to build a cross-strain consensus.

#### Step 0 — Building the pair list

- **Prediction CSV given** → pairs and their scores come from it.
- **Prediction CSV blank** → pairs are built from the edges file: each edge with
  exactly one endpoint whose name starts with the **sRNA prefix** (default
  `sRNA`) becomes an sRNA→target pair; the energy/probability step is then
  skipped (no scores available).

#### Step 1 — Energy / probability filter

Each metric can be individually enabled/disabled, has its own threshold, and is
skipped automatically if its column is absent. Instead of requiring *all* enabled
metrics (AND), you can require only **M of N** to pass.

By default only the four calibrated consensus metrics are enabled (the
four-tool AND rule); the four energy/`Probability_TargetRNA3` metrics are
shipped disabled and can be re-enabled if needed.

| Option | Default | Meaning |
|--------|---------|---------|
| `E_intaRNA max (<=)` | **disabled**, `-2.44` | Keep IntaRNA energy at or below this. |
| `E_Rnaplex max (<=)` | **disabled**, `-32.6` | Keep RNAplex energy at or below this. |
| `E_TargetRNA3 max (<=)` | **disabled**, `-5.13` | Keep TargetRNA3 energy at or below this. |
| `Probability_TargetRNA3 min (>=)` | **disabled**, `0.06` | Keep TargetRNA3 probability at or above this. |
| `Probability_sRNARFTarget min (>=)` | enabled, `0.45` | Keep sRNARFTarget probability at or above this. |
| `p_intaRNA max (<=)` | enabled, `0.41` | Keep IntaRNA p-value at or below this. |
| `p_Rnaplex max (<=)` | enabled, `0.76` | Keep RNAplex p-value at or below this. |
| `p_TargetRNA3 max (<=)` | enabled, `0.89` | Keep TargetRNA3 p-value at or below this. |
| `Min metrics that must pass` | `All enabled (AND)` | Require this many enabled metrics to pass (M-of-N); `0` = all. |

#### Step 2 — DEG consistency & direction

| Option | Default | Meaning |
|--------|---------|---------|
| `padj cutoff` | `0.05` | A gene must have `padj ≤` this in a file to count there. |
| `Min baseMean (>=)` | **enabled**, `10.0` | Minimum DESeq2 `baseMean` for a gene to count as a DEG. |
| `Keep log2FC >=` / `Keep log2FC <=` | **both disabled** | Optional log2FC inclusion band; off = keep any log2FC. |
| `Consensus type` | `presence (ignore log2FC)` | `sign of log2FC` = agree on direction; `presence` = called in enough files regardless of sign. |
| `Consensus: min files (X)` | `5` | Number of DEG files a gene must be called in to reach consensus. **Lower it to match how many DEG files you actually load** (e.g. `2` for the bundled test). |
| `Required sRNA direction` | off (`upregulated`) | Restrict to `upregulated` / `downregulated` sRNAs. |
| `Required target direction` | off (`upregulated`) | Restrict to `upregulated` / `downregulated` targets. |
| `sRNA vs target direction` | off (`same`) | Require the pair to move in the `same` or `opposite` direction. |

Direction is per gene: `log2FC > 0` ⇒ upregulated, `< 0` ⇒ downregulated.

#### Step 3 — sRNA genomic location filter

Restrict results to sRNAs of chosen categories, using the annotation table from
[01 Position Classification](#01-position-classification). sRNAs absent from the
annotation file are dropped.

| Option | Default | Meaning |
|--------|---------|---------|
| `Enable location filter` | off | Turn the location filter on. |
| `Annotation file` | — | The `sRNA_annotation.xlsx`/`.csv` from step 01. |
| `sRNA ID column` | `exon_id` | Column in the annotation holding the sRNA name. |
| `Keep locations` | all ticked | Any of `Intragenic`, `Antisense`, `5' UTR`, `3' UTR`, `Intergenic`. |

#### Step 4 — WGCNA module relationship

| Option | Default | Meaning |
|--------|---------|---------|
| `sRNA / target module` | `same` | Keep pairs in the **same** module, in **different** modules, or **any** (skip). |

Requires the nodes file; pairs where either gene lacks a module are dropped when
the mode is `same`/`different`.

#### Step 5 — Network edge weight & selection

| Option | Default | Meaning |
|--------|---------|---------|
| `Interaction must have a network edge` | on | Drop pairs with no co-expression edge/weight. |
| `Min edge weight (>=)` | **disabled**, `0.0` | Drop pairs whose edge weight is below this. |
| `Keep best per` | **enabled**, `best_per_target` | `best_per_target` (one sRNA per target), `best_per_srna` (one target per sRNA). Untick to keep all. |
| `Top N per group` | `0` (single best) | With selection on, keep the top N partners per group by edge weight. |

"Best" = the pair with the highest network edge `weight`. Because
`best_per_target` keeps only the single strongest sRNA per unique target, this is
the step that collapses many surviving pairs down to one representative per
target.

#### Output

Default file: `filtered_weight.csv` (editable). Tick **Generate all intermediate
tables** to also save the result of every filtering step next to the output CSV,
named with a step suffix so each stage of the funnel can be inspected on its own:

| File | Content |
|------|---------|
| `filtered_weight_step1_base.csv` | Starting sRNA–target pairs. |
| `filtered_weight_step2_energy.csv` | After the energy/probability consensus. |
| `filtered_weight_step3_deg.csv` | After the DEG filter. |
| `filtered_weight_step4_wgcna.csv` | After same-module + require-edge (+ min-weight), **before** best-pair selection. |
| `filtered_weight_step5_bestpair.csv` | After `best_per_target` / `best_per_srna` prioritization. |
| `filtered_weight_step6_position.csv` | After the genomic-location filter (final). |

A step that is skipped still writes a table (identical to the previous step) so
the full funnel is always captured.

Each retained interaction carries the original `sRNA`/`Target`, whichever score
columns were present, and the annotations added during filtering:

| Column(s) | Added by |
|-----------|----------|
| `sRNA`, `Target` | input |
| `E_*`, `Probability_*` | input prediction table |
| `sRNA_DEG`, `Target_DEG` | Step 2 (`upregulated` / `downregulated`) |
| `sRNA_Location` | Step 3 (only if the location filter ran) |
| `sRNA_Module`, `Target_Module` | Step 4 (only if the module filter ran) |
| `weight` | Step 5 (co-expression edge weight) |

> **Notes:** empirical p-values for IntaRNA/RNAplex are computed in
> [04 Combine Predictions](#04-combine-predictions), not here; the location
> filter drops sRNAs not found in the annotation; and no multiple-testing
> correction is applied across interactions (filtering is threshold-based).

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
2. Click **Run Position Classification**.

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
3. Adjust the Rscript executable path if required.
4. Click **Run WGCNA**.

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
4. Adjust the Rscript executable path if required.
5. Click **Run DESeq2**.

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
5. The **Consensus: min files (X)** value now defaults to `5` — the calibrated
   whole-study setting (six RNA-seq datasets). **This small test ships only two
   DEG files, so change it to `2`.** If you leave it at `5`, no gene can be
   present in five files and the DEG filter empties the table.
6. *(Optional)* Tick **Generate all intermediate tables** to also write one CSV
   per filtering step (`filtered_weight_step1_base.csv` through
   `_step6_position.csv`) next to the Output CSV.

**Expected output (using `Prediction_combined.csv`, `min files (X) = 2`, no
location filter):**

```
Total pairs: 32789
[1] ENERGY / PROBABILITY  (>= 4/4 of Probability_sRNARFTarget,
    p_intaRNA, p_Rnaplex, p_TargetRNA3): removed 21639, kept 11150
[2] DEG  (padj <= 0.05, baseMean >= 10, presence in >= 2 files): 129
[3] WGCNA  module (same): 77; require network edge: 64
    best_per_target (single best per target): 59
Final interactions: 59
```

*(With the bundled `Prediction_test.csv` the starting total is 33996 instead of
32789, so the intermediate counts differ slightly.)*

#### Optional — sRNA location filter

You can restrict the results to sRNAs of a chosen genomic category using the
annotation table from step 1. In the **sRNA location filter** group of the
Filters tab:

1. Tick **Enable location filter**.
2. **Annotation file:** add `sRNA_annotation.xlsx` (produced in step 1).
3. Under **Keep locations**, leave only **Antisense** ticked.

**Expected output (Antisense only, `min files (X) = 2`):**

```
Total pairs: 32789
[1] ENERGY / PROBABILITY  (>= 4/4 metrics): 11150
[2] DEG  (presence in >= 2 files): 129
[3] WGCNA  module (same): 77; require network edge: 64
    best_per_target (single best per target): 59
[4] POSITION CLASSIFICATION  keep Antisense: removed 37 -> 22
Final interactions: 22
```

The genomic-location filter always runs **last**, so it is applied after the
`best_per_target` selection (59 pairs), keeping the 22 whose sRNA is `Antisense`.
sRNAs not found in the annotation file are dropped.

If your numbers match, the pipeline is installed correctly and you are ready to
run it on your own data.

#### Large-scale run — reproducing the published results

The bundled `data/` set is a small subset for smoke-testing. To reproduce the
values reported in the article this program accompanies, run the Filters step on
the full-scale inputs instead:

1. Download the full input dataset (prediction table, six DEG files, WGCNA
   nodes/edges) from:
   <https://drive.google.com/drive/folders/166prrrQZ9MAkCUC_nlXPUnniiOHwWW2K?usp=sharing>
2. In **05 Filters**, keep the calibrated defaults. Because the full study uses
   **six** RNA-seq datasets, the default **Consensus: min files (X) = 5** is the
   correct value here (do **not** lower it to `2` as in the small test).
3. Tick **Generate all intermediate tables** to write one CSV per step. A
   reference copy of the expected step tables is available at:
   <https://drive.google.com/drive/folders/1_W_4aiV7zdxE4kJcBU1Uvc5Cw-NFOXQq?usp=sharing>

With those defaults (four-tool energy/probability AND consensus; DEG `padj ≤
0.05`, `baseMean ≥ 10`, no log2FC band, presence in `≥ 5` of the six files; same
WGCNA module; require a network edge; `best_per_target`; no location filter), the
run reproduces the following funnel:

| Step | File | Pairs kept |
|------|------|-----------:|
| Base pairs (prediction table) | `filtered_weight_step1_base.csv` | 1,731,060 |
| Energy/probability (≥ 4/4 metrics, AND) | `filtered_weight_step2_energy.csv` | 949,716 |
| DEG (padj ≤ 0.05, baseMean ≥ 10, presence in ≥ 5 files) | `filtered_weight_step3_deg.csv` | 149,545 |
| WGCNA same-module + require network edge | `filtered_weight_step4_wgcna.csv` | 8,618 |
| `best_per_target` (single best per target) | `filtered_weight_step5_bestpair.csv` | 1,290 |
| Position filter (OFF here) | `filtered_weight_step6_position.csv` | 1,290 |

**Final interactions: 1,290** (`filtered_weight.csv`).

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
  prediction table by [04 Combine Predictions](#04-combine-predictions)

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
