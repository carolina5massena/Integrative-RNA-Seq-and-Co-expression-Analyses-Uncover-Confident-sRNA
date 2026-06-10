#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# install_packages.R
# Installs every R package the pipeline needs, idempotently.
# Safe to re-run: already-installed packages are skipped.
#
# Run directly:   Rscript install_packages.R
# (the setup scripts call this for you)
# ---------------------------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Packages pulled from CRAN
cran_pkgs <- c("tibble", "fastDummies", "flashClust")

# Packages pulled from Bioconductor.
# impute, preprocessCore and GO.db are WGCNA dependencies that are not always
# pulled in automatically, so we install them explicitly (before WGCNA).
bioc_pkgs <- c("impute", "preprocessCore", "GO.db", "WGCNA", "DESeq2")

message("== Checking CRAN packages ==")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing CRAN package: %s", pkg))
    install.packages(pkg)
  } else {
    message(sprintf("Already installed: %s", pkg))
  }
}

message("== Checking Bioconductor packages ==")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager")
  install.packages("BiocManager")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing Bioconductor package: %s", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    message(sprintf("Already installed: %s", pkg))
  }
}

# ---- Verify everything loads -------------------------------------------------
message("== Verifying installation ==")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing <- all_pkgs[!vapply(all_pkgs,
                            function(p) requireNamespace(p, quietly = TRUE),
                            logical(1))]
if (length(missing) > 0) {
  stop(sprintf("FAILED to install: %s", paste(missing, collapse = ", ")))
}
message("All R packages installed and importable. R: ",
        R.version.string)
