@echo off
chcp 65001 >nul
cd /d "%~dp0"
title Wald Brandenburg - Zusammenfassung
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
echo Fasse alle Kacheln zusammen und zeichne die Brandenburg-Karte (ca. 1 Minute) ...
"%PY%" bb_wald.py --out ergebnisse --karte
echo.
pause
