# Installation & Setup

This pipeline has two runtimes: **Python** (the GUI) and **R** (the WGCNA and
DESeq2 analyses). The scripts below install both and wire up the paths
automatically, so the same steps work on any computer.

## Prerequisites (install these once, per machine)

1. **Python 3.10+** — https://www.python.org/downloads/
   On Windows, tick **"Add Python to PATH"** during install.
2. **R 4.1+** — https://cran.r-project.org/
   You do **not** need to add R to PATH — the setup finds it for you.

## One-click setup

### Windows
Double-click **`setup.bat`** (or run it in a terminal). It will:
- create a local virtual environment `.venv`,
- install the Python packages from `requirements.txt`,
- locate `Rscript.exe` automatically and install the R packages.

Then launch the app by double-clicking **`run.bat`**.

### macOS / Linux
```bash
bash setup.sh
bash run.sh
```

That's it. Re-running setup is safe — it skips anything already installed.

## What gets installed

**Python** (`requirements.txt`): PySide6, pandas, numpy.

**R** (`install_packages.R`):
- CRAN: tibble, fastDummies, flashClust
- Bioconductor: WGCNA, DESeq2

## How the R path is configured automatically


1. the `RSCRIPT` environment variable, if set to a valid file;
2. `Rscript` on the system PATH;
3. common install locations (newest R version wins on Windows, e.g.
   `C:\Program Files\R\R-4.6.0\bin\Rscript.exe`);
4. falling back to the bare name `Rscript`.

The resolved path is pre-filled in the **"Rscript executable:"** field of both
the WGCNA and DESeq2 tabs. You can still override it manually (Browse button) if
R is installed somewhere unusual, or set the `RSCRIPT` environment variable to
force a specific R version.

## Manual install (if you prefer not to use the scripts)

```bash
# Python
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Windows
# .venv/bin/python -m pip install -r requirements.txt      # macOS/Linux

# R packages
Rscript install_packages.R

# Run
.venv/Scripts/python GUI/app.py    # Windows
# .venv/bin/python GUI/app.py        # macOS/Linux
```

## Troubleshooting

- **`FileNotFoundError: [WinError 2]` when running an analysis** — R isn't being
  found. Run `setup.bat` again, or paste the full path to `Rscript.exe` into the
  GUI field. Find it with PowerShell:
  `Get-ChildItem "C:\Program Files\R" -Recurse -Filter Rscript.exe`
- **`there is no package called 'WGCNA'` / `'DESeq2'`** — the R packages didn't
  install. Run `Rscript install_packages.R` and watch for errors (Bioconductor
  installs can take several minutes the first time).
- **`python` not recognized** — Python isn't on PATH. Reinstall it with the
  "Add to PATH" option, or use the `py` launcher.
