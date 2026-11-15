@echo off
cd /d "%~dp0"
mkdir grumsin 2>nul
cd grumsin
echo Lade Laserscans und Luftbilder (Grumsin, ca. 500 MB) ...
for %%T in (33425-5871 33424-5870) do (
  curl.exe -L --fail -o als_%%T.zip https://data.geobasis-bb.de/geobasis/daten/als/laz/als_%%T.zip
  curl.exe -L --fail -o dop_%%T.zip https://data.geobasis-bb.de/geobasis/daten/dop/rgbi_tif/dop_%%T.zip
)
dir
echo Fertig.
pause
