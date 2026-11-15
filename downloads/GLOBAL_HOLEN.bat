@echo off
cd /d "%~dp0"
mkdir wsl 2>nul
curl.exe -L --fail -o wsl\demol_global_AGB_TLS.xlsx "https://zenodo.org/records/5236762/files/demol_et_al_AGB_with_TLS.xlsx?download=1"
dir wsl
pause
