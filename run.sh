#!/usr/bin/env bash
# Launch the pipeline GUI using the local virtual environment.
set -e
cd "$(dirname "$0")"
if [ ! -x ".venv/bin/python" ]; then
  echo "Virtual environment not found. Run:  bash setup.sh"
  exit 1
fi
exec ".venv/bin/python" "GUI/app.py"
