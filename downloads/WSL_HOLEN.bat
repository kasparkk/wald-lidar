@echo off
cd /d "%~dp0"
mkdir wsl 2>nul
cd wsl
echo Lade WSL-Baumdaten (Tabelle + 3 Punktwolken-Pakete) ...
curl.exe -L --fail -o tree_parameters.csv "https://www.envidat.ch/dataset/c82b9485-0b40-41f4-a798-4e4e725cafdc/resource/450e6bd6-f0c2-4c31-b410-52013c50e988/download/tree_parameters_tls_ad_wl.csv"
for %%F in (site1_WL_group1 site1_WL_group2 site2_AD) do curl.exe -L --fail -C - -o %%F.zip https://zhc-c.s3.cloud.switch.ch/pointclouds/20230516_TLSPointClouds_volumeEstimation/data/%%F.zip
dir
echo Fertig.
pause
