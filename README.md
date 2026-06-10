# Integrative RNA-Seq and Co-expression Analyses Uncover Confident sRNA

## Overview

The workflow combines sequence-based predictions with transcriptomic and network-level evidence to produce high-confidence regulatory interactions.

A more complete description of each of the modes and their use is available in the readme file within each module.

<img width="606" height="340" alt="image" src="https://github.com/user-attachments/assets/8f1a853f-b106-4f29-a660-7e94f5b9712e" />

*Note: The processes shown in the figure with a dotted outline were not included in this repository.*

---
#Test Data and Expected Results

A tutorial explaining how to use this program, along with the expected results, is provided in the PDF file named “Run first test.”

---

# Optional External Tools

## Interaction Prediction

- **IntaRNA**  
  https://github.com/BackofenLab/IntaRNA  

- **RNAplex (ViennaRNA Package)**  
  https://www.tbi.univie.ac.at/RNA/  

- **TargetRNA3**  
  https://cs.wellesley.edu/~btjaden/TargetRNA3 

- **sRNARFTarget**  
  https://github.com/BioinformaticsLabAtMUN/sRNARFTarget


---


# Software Requirements

The pipeline uses two runtimes: **Python ≥ 3.10** (for the GUI and data steps)
and **R ≥ 4.1** (for the WGCNA and DESeq2 analyses). Install those two once per
machine, then use the setup scripts below — they handle every package and
configure the R path automatically. See [`INSTALL.md`](INSTALL.md) for full
details and troubleshooting.

## Quick start (recommended)

Install [Python](https://www.python.org/downloads/) (tick *"Add to PATH"* on
Windows) and [R](https://cran.r-project.org/), then run the one-click setup:

**Windows** — double-click `setup.bat`, then `run.bat` to launch the GUI.

**macOS / Linux**
```
bash setup.sh
bash run.sh
```

The setup creates a local `.venv`, installs the Python packages from
`requirements.txt`, locates `Rscript` automatically, and installs the R packages
via `install_packages.R`. It is safe to re-run.

## Manual install

```
# Python (PySide6, pandas, numpy)
pip install -r requirements.txt

# R: CRAN (tibble, fastDummies, flashClust) + Bioconductor (WGCNA, DESeq2)
Rscript install_packages.R
```

## GUI

A PySide6 GUI (`GUI/app.py`) wraps the WGCNA and DESeq2 steps. The path to
`Rscript` is detected automatically and pre-filled in each tab — no manual
configuration is needed on a typical install. See [`GUI/README.md`](GUI/README.md).

## Troubleshooting

**`OSError: [Errno 2] No such file or directory` / "enable long paths" during
`pip install` (Windows).** This is Windows' 260-character path limit being hit
while installing PySide6 (which has deeply nested files), made worse if the
project lives in a long folder path. Enable long-path support, then re-run setup:

1. Open **PowerShell as Administrator** and run:
   ```powershell
   New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" -Name "LongPathsEnabled" -Value 1 -PropertyType DWORD -Force
   ```
2. **Restart the computer** for the setting to take effect.
3. Delete the partially-created `.venv` folder in the project directory.
4. Run `setup.bat` again.

   Alternatively, avoid the limit by moving the project to a short path such as
   `C:\srna` (so files end up at `C:\srna\GUI\app.py`) before running setup.

```


# Contact

For issues or reproducibility questions, please open a GitHub issue in this repository.
