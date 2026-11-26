@echo off
chcp 65001 >nul
cd /d "%~dp0"
title Wald Brandenburg - RCT-Test (holen)
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -c "import numpy, scipy, laspy, lazrs" >nul 2>&1
if errorlevel 1 (
  echo Installiere Bibliotheken ...
  "%PY%" -m pip install --quiet --no-warn-script-location numpy scipy "laspy[lazrs]"
)
"%PY%" rct_holen.py
if errorlevel 1 (
  echo.
  echo Da ist etwas schiefgegangen. Bitte einen Screenshot von diesem Fenster in den Chat schicken.
)
pause
