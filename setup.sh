#!/usr/bin/env bash
# ===========================================================================
#  setup.sh  -  One-click setup for macOS / Linux
#  - creates a local Python virtual environment (.venv)
#  - installs Python dependencies (requirements.txt)
#  - locates Rscript and installs the R/Bioconductor packages
#  Usage:  bash setup.sh
# ===========================================================================
set -e
cd "$(dirname "$0")"

echo
echo "=== [1/4] Locating Python ==========================================="
PY=""
for c in python3 python; do
  if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
  echo "ERROR: Python 3 not found. Install it from https://www.python.org/downloads/"
  exit 1
fi
echo "Using: $PY ($($PY --version))"

echo
echo "=== [2/4] Creating virtual environment (.venv) ====================="
if [ ! -x ".venv/bin/python" ]; then
  "$PY" -m venv .venv
else
  echo ".venv already exists - reusing it."
fi
VENV_PY=".venv/bin/python"

echo
echo "=== [3/4] Installing Python packages ==============================="
"$VENV_PY" -m pip install --upgrade pip
"$VENV_PY" -m pip install -r requirements.txt

echo
echo "=== [4/4] Locating Rscript and installing R packages =============="
RSCRIPT="$(command -v Rscript || true)"
if [ -z "$RSCRIPT" ]; then
  for cand in \
    /Library/Frameworks/R.framework/Resources/bin/Rscript \
    /opt/homebrew/bin/Rscript \
    /usr/local/bin/Rscript \
    /usr/bin/Rscript; do
    if [ -x "$cand" ]; then RSCRIPT="$cand"; break; fi
  done
fi
if [ -z "$RSCRIPT" ]; then
  echo "WARNING: Rscript not found. Install R from https://cran.r-project.org/"
  echo "Then run:  Rscript install_packages.R"
else
  echo "Using Rscript: $RSCRIPT"
  "$RSCRIPT" install_packages.R || echo "WARNING: R package install reported an error."
fi

echo
echo "=== Setup complete ================================================="
echo "Launch the GUI with:   bash run.sh"
