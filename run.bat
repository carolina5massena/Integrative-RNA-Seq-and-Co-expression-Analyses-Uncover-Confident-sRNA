@echo off
REM Launch the pipeline GUI using the local virtual environment.
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Virtual environment not found. Run setup.bat first.
  pause & exit /b 1
)
".venv\Scripts\python.exe" "GUI\app.py"
if errorlevel 1 pause
