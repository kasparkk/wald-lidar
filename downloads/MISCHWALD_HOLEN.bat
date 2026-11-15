@echo off
cd /d "%~dp0"
mkdir mischwald 2>nul
cd mischwald
echo Lade 4 Laserscan-Kacheln mit Kiefer-Laubholz-Mischung (ca. 600 MB) ...
for %%T in (33325-5766 33465-5780 33417-5805 33337-5726) do curl.exe -L --fail -C - -o als_%%T.zip https://data.geobasis-bb.de/geobasis/daten/als/laz/als_%%T.zip
dir
echo Fertig.
pause
