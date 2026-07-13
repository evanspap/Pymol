@echo off
REM RMSD Calculator Script Runner
REM This script runs the RMSD calculation with proper Python environment

echo ============================================================
echo RMSD Calculator for Multi-Frame PDB Files
echo ============================================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo Error: Python is not installed or not in PATH
    echo Please install Python 3.7+ from https://www.python.org
    pause
    exit /b 1
)

REM Check if NumPy is installed
python -c "import numpy" >nul 2>&1
if errorlevel 1 (
    echo NumPy not found. Installing...
    python -m pip install numpy
)

echo.
echo Running RMSD calculation...
echo.

REM Run the script
python "%~dp0calculate_rmsd.py"

echo.
echo Script completed. Results saved to rmsd_results.csv
pause
