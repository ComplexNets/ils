@echo off
echo Checking Python environment...

REM Activate conda environment using full path
call C:\Users\X1\miniconda3\Scripts\activate.bat ils 2>nul
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
