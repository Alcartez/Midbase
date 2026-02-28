@echo off
REM XY Counsel - Complete DE Recalculation and Export Pipeline
REM Run this script to recalculate all differential expression statistics

echo ================================================
echo XY Counsel - DE Recalculation Pipeline
echo ================================================
echo.

set PYTHON_EXE=C:\Users\user\miniconda3\envs\xy_counsel\python.exe
set PROJECT_DIR=d:\Projects\Merge_Midbase_Serenova\XY_Counsel

cd /d %PROJECT_DIR%

echo [Step 1/3] Recalculating within-study DE (3 studies with controls)...
echo.
%PYTHON_EXE% recalculate_de_manual.py
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Within-study DE calculation failed!
    pause
    exit /b 1
)

echo.
echo.
echo [Step 2/3] Calculating cross-study DE (5 comparisons with pooled controls)...
echo.
%PYTHON_EXE% cross_study_de.py
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Cross-study DE calculation failed!
    pause
    exit /b 1
)

echo.
echo.
echo [Step 3/3] Exporting DE results to CSV files...
echo.
%PYTHON_EXE% export_de_results.py
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: DE export failed!
    pause
    exit /b 1
)

echo.
echo ================================================
echo SUCCESS! All DE calculations complete.
echo ================================================
echo.
echo Results saved to: %PROJECT_DIR%\de_results\
echo Database updated: %PROJECT_DIR%\midbase.db
echo.
echo You can now launch the Streamlit platform with proper p-values:
echo   streamlit run Home.py
echo.
pause
