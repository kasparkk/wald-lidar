@echo off
chcp 65001 >nul
cd /d "%~dp0"
title Wald Brandenburg - Baumart und Vitalitaet fuer Hammer
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -c "import tifffile, sklearn" >nul 2>&1
if errorlevel 1 (
  echo Installiere zusaetzliche Bibliotheken ...
  "%PY%" -m pip install --quiet --no-warn-script-location tifffile scikit-learn
)
echo.
echo Baumart und Vitalitaet fuer 64 Kacheln um Hammer.
echo Laedt ca. 5,5 GB Luftbilder (werden nach jeder Kachel geloescht). Dauer ca. 1-3 Stunden.
echo.
"%PY%" bb_extra.py --out ergebnisse --bbox 405 5774 412 5781
if errorlevel 1 (
  echo.
  echo Da ist etwas schiefgegangen. Bitte einen Screenshot von diesem Fenster in den Chat schicken.
  pause
  exit /b 1
)
echo.
echo Fertig. Im Viewer "Neu einlesen" klicken.
pause
