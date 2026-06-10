@echo off
REM ===========================================================================
REM  setup.bat  -  One-click setup for Windows
REM  - creates a local Python virtual environment (.venv)
REM  - installs the Python dependencies (requirements.txt)
REM  - locates Rscript automatically and installs the R/Bioconductor packages
REM  Just double-click this file, or run it from a terminal.
REM ===========================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo.
echo === [1/4] Locating Python ===========================================
set "PY="
where py        >nul 2>nul && set "PY=py"
if not defined PY ( where python >nul 2>nul && set "PY=python" )
if not defined PY (
  echo ERROR: Python was not found on PATH.
  echo Install Python 3.10+ from https://www.python.org/downloads/  ^(tick "Add to PATH"^)
  pause & exit /b 1
)
echo Using Python launcher: !PY!

echo.
echo === [2/4] Creating virtual environment ^(.venv^) =====================
if not exist ".venv\Scripts\python.exe" (
  !PY! -m venv .venv || ( echo ERROR: failed to create venv & pause & exit /b 1 )
) else (
  echo .venv already exists - reusing it.
)
set "VENV_PY=.venv\Scripts\python.exe"

echo.
echo === [3/4] Installing Python packages ================================
"!VENV_PY!" -m pip install --upgrade pip
"!VENV_PY!" -m pip install -r requirements.txt || ( echo ERROR: pip install failed & pause & exit /b 1 )

echo.
echo === [4/4] Locating Rscript and installing R packages ===============
set "RSCRIPT="
for /f "delims=" %%i in ('where Rscript 2^>nul') do set "RSCRIPT=%%i"
if not defined RSCRIPT (
  REM Search standard install dirs, newest version first ^(/o-n = sort by name desc^)
  for %%B in ("C:\Program Files\R" "C:\Program Files (x86)\R" "%LOCALAPPDATA%\Programs\R") do (
    if exist "%%~B" for /f "delims=" %%i in ('dir /b /o-n "%%~B" 2^>nul') do (
      if not defined RSCRIPT if exist "%%~B\%%i\bin\Rscript.exe" set "RSCRIPT=%%~B\%%i\bin\Rscript.exe"
    )
  )
)
if not defined RSCRIPT (
  echo WARNING: Rscript not found. Install R from https://cran.r-project.org/bin/windows/base/
  echo Then re-run this script, or run:  Rscript install_packages.R
) else (
  echo Using Rscript: !RSCRIPT!
  "!RSCRIPT!" install_packages.R || echo WARNING: R package install reported an error - see messages above.
)

echo.
echo === Setup complete =================================================
echo Launch the GUI with:   run.bat
echo.
pause
