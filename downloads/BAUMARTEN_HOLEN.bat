@echo off
cd /d "%~dp0"
mkdir baumarten 2>nul
echo Lade Baumartenkarte Deutschland (Koeltzow, FU Berlin, Zenodo 15655451, ca. 690 MB) ...
curl.exe -L --fail -C - -o baumarten\Baumartenkarte_Deutschland_masked.tif "https://zenodo.org/records/15655451/files/Baumartenkarte_Deutschland_masked.tif?download=1"
dir baumarten
echo Fertig.
pause
