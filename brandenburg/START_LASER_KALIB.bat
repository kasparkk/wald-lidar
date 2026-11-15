@echo off
chcp 65001 >nul
cd /d "%~dp0"
title Wald Brandenburg - Laser-Stichprobe
if not exist "%~dp0bb_laser_kalib.py" (
  echo bb_laser_kalib.py liegt nicht neben dieser Datei.
  pause
  exit /b 1
)
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -c "import numpy, scipy, laspy, lazrs" >nul 2>&1
if errorlevel 1 (
  echo Installiere Bibliotheken ...
  "%PY%" -m pip install --quiet --no-warn-script-location numpy scipy "laspy[lazrs]"
)
echo.
echo Laser-Stichprobe: 50 Waldkacheln, je ca. 135 MB Download, wird sofort wieder geloescht.
echo Dauer grob 1 bis 2 Stunden. Abbrechen mit Strg+C, Neustart setzt fort.
echo.
"%PY%" bb_laser_kalib.py --out ergebnisse --n 50 --workers 2
if errorlevel 1 (
  echo.
  echo Da ist etwas schiefgegangen. Bitte einen Screenshot von diesem Fenster in den Chat schicken.
  pause
  exit /b 1
)
echo.
echo Fertig. Ergebnis: ergebnisse\laser_kalib\kalibrierung.json
pause
