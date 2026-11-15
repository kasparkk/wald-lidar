@echo off
chcp 65001 >nul
cd /d "%~dp0"
title Wald Brandenburg - Test um Hammer
set "PY=python"

%PY% -c "import sys" >nul 2>&1
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

echo.
echo Installiere Bibliotheken ...
"%PY%" -m pip install --quiet --upgrade pip
"%PY%" -m pip install --quiet numpy scipy matplotlib "laspy[lazrs]"

echo.
echo Starte Test: 49 Kacheln um Loepten/Hammer, ca. 10 Minuten ...
"%PY%" bb_wald.py --out ergebnisse --bbox 405 5774 412 5781

echo.
echo Fertig. Bitte kacheln.csv und uebersicht.png aus dem Ordner "ergebnisse" in den Chat ziehen.
if exist ergebnisse\uebersicht.png start "" "ergebnisse\uebersicht.png"
explorer ergebnisse
pause
