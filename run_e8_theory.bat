@echo off
REM E8 Theory of Everything - Windows Launcher
REM This script sets UTF-8 encoding before running the simulation

echo ========================================
echo E8 Theory of Everything v2.1
echo Windows UTF-8 Launcher
echo ========================================
echo.

REM Set Python to use UTF-8 encoding for input/output
set PYTHONIOENCODING=utf-8

REM Run the simulation with all arguments passed through
python run_unified_theory.py %*

REM Check if successful
if %ERRORLEVEL% == 0 (
    echo.
    echo ========================================
    echo Simulation completed successfully!
    echo ========================================
) else (
    echo.
    echo ========================================
    echo ERROR: Simulation failed with code %ERRORLEVEL%
    echo ========================================
)
