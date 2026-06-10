# sRNA Pipeline GUI (PySide6)

A desktop application that wraps the three analysis steps of this repository
behind a tabbed interface. Every user input is editable in the UI.

Each tab is self-contained: all of its input files and parameters are configured
on that tab.

## Tabs

1. **01 Position Classification** — GTF file, feature/ID keys, sRNA & gene
   prefixes, rename prefixes, overlap threshold, UTR window, output file. Pure
   Python (no R needed).
2. **02 Co-expression analysis** — all WGCNA `USER INPUTS` (input files, column
   names, soft power, module size, deep split, merge cut height, Cytoscape
   threshold, threads, network/correlation type). Generates a configured
   `WGCNA_configured.R` and runs it via `Rscript`.
3. **Differential Expression** — runs DESeq2 (R) on a raw count matrix + sample
   metadata. Set the gene-ID/sample/condition columns, reference vs. test level,
   and a minimum-count filter. Writes a TSV with `gene_id`, `log2FoldChange`,
   `padj` (plus baseMean/pvalue/lfcSE/stat) that plugs straight into the Filters
   tab. Generates a configured `DESeq2_configured.R` and runs it via `Rscript`.
4. **03 Filters** — prediction CSV, one or more DEG TSVs (comma-separated),
   WGCNA nodes/edges files, padj cutoff, min strains, and all five
   energy/probability thresholds. Pure Python.

Each analysis tab has its own **Run** button. Work runs in a background thread and
streams a log; a dialog reports success/failure.

## Install

```bash
pip install PySide6 pandas numpy openpyxl
```

Step 02 also requires R (≥ 4.1) with `WGCNA`, `tibble`, `fastDummies`, `flashClust`
installed; the Differential Expression tab requires R with Bioconductor `DESeq2`.
In both cases `Rscript` must be on your PATH (or point to it in the tab).

## Run

```bash
cd GUI
python app.py
```

## Files

- `app.py` — PySide6 interface.
- `pipeline_logic.py` — parameterized implementations of the three steps.
