@echo off
echo Checking Python environment...

REM Activate conda environment if it exists
call conda activate ils 2>nul
if %errorlevel% equ 0 (
    echo Successfully activated 'ils' conda environment
) else (
    echo Warning: Could not activate 'ils' conda environment
    echo Make sure you have the correct environment activated
    pause
    exit /b 1
)

echo Starting ILS UI...
cd frontend
streamlit run ils_ui.py
pause
