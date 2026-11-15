@echo off
chcp 65001 >nul
cd /d "%~dp0"
title Wald Brandenburg - Teil 1 (West)

if not exist "%~dp0bb_wald.py" (
  echo bb_wald.py liegt nicht neben dieser Datei.
  echo Bitte START.bat im entpackten Ordner Downloads\bb_wald starten, nicht direkt aus der ZIP.
  pause
  exit /b 1
)

set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"

"%PY%" -c "import sys" >nul 2>&1
if errorlevel 1 (
  echo Python fehlt - wird jetzt installiert. Windows fragt eventuell nach einer Bestaetigung.
  winget install -e --id Python.Python.3.12 --architecture x64 --scope user --accept-package-agreements --accept-source-agreements
  set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
)

"%PY%" -c "import sys" >nul 2>&1
if errorlevel 1 (
  echo.
  echo Python konnte nicht automatisch installiert werden.
  echo Bitte von python.org den "Windows installer (64-bit)" installieren und START.bat nochmal doppelklicken.
  pause
  exit /b 1
)

"%PY%" -c "import numpy, scipy, matplotlib, laspy, lazrs" >nul 2>&1
if errorlevel 1 (
  echo Installiere Bibliotheken ...
  "%PY%" -m pip install --quiet --no-warn-script-location numpy scipy matplotlib "laspy[lazrs]"
)

echo.
echo Teil 1: Kacheln Ost 250 bis 436 km. Abbrechen mit Strg+C, Neustart setzt fort.
"%PY%" bb_wald.py --out ergebnisse --bbox 250 5600 436 6000
if errorlevel 1 (
  echo.
  echo Da ist etwas schiefgegangen. Bitte einen Screenshot von diesem Fenster in den Chat schicken.
  pause
  exit /b 1
)

echo.
echo Fertig. Ergebnisse liegen im Ordner "ergebnisse".
if exist ergebnisse\uebersicht.png start "" "ergebnisse\uebersicht.png"
pause
