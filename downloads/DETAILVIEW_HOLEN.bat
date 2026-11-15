@echo off
cd /d %USERPROFILE%\Downloads
if not exist detailview mkdir detailview
cd detailview
echo Lade DetailView-Gewichte (FreiData, Uni Freiburg, je ca. 240 MB) ...
curl -L -C - --retry 5 -o model_ft_202412171652_3 "https://freidata.uni-freiburg.de/records/f850a-bb152/files/model_ft_202412171652_3?download=1"
curl -L -C - --retry 5 -o model_202305171452_60 "https://freidata.uni-freiburg.de/records/xw42t-6mt03/files/model_202305171452_60?download=1"
dir
echo Fertig. Fenster kann geschlossen werden.
pause
