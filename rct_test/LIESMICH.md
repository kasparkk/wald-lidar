# Test: RayCloudTools (3Dtrees.earth) gegen nDOM und eigene Laserzählung

Ordner `rct_test` in den Ordner `wald-lidar` legen (neben `brandenburg`), sonst finden die Skripte ihre Funktionen nicht.

## 1. Kacheln holen
`START_RCT_HOLEN.bat` doppelklicken. Lädt drei Kacheln aus dem Raum Hammer/Löpten (je ca. 135 MB):

| Kachel | warum |
|---|---|
| 33411-5775 | dicht, geschlossen, fast reine Kiefer (Kronenschluss 0,90) |
| 33408-5778 | lockerer, etwas Laubholz (0,67, Laub ca. 14 %) |
| 33409-5779 | hoch (Median 22,4 m), offen, etwas Laubholz |

## 2. Bei 3Dtrees hochladen
Die drei Dateien aus `upload/` bei 3Dtrees hochladen und mit RayCloudTools segmentieren lassen.
Ergebnisse (LAZ mit Spalte `PredInstance`) herunterladen und in `ergebnis/` legen. Dateinamen egal.

## 3. Vergleichen
`START_RCT_VERGLEICH.bat` doppelklicken. Ergebnis: `ergebnis/vergleich.csv`, je Kachel:

- Bäume je ha nach RCT, eigener Laserzählung und nDOM (nur Wald, nur ab 6 m)
- `faktor_rct_zu_ndom`: direkt vergleichbar mit dem Kalibrierfaktor 1,84 aus `bb_laser_kalib.py`
- `ndom_in_rct_gefunden`: Anteil nDOM-Bäume, die RCT auch hat (Spitzen höchstens 2,5 m auseinander)
- Median Baumhöhe und Kronendurchmesser nach RCT

Dazu `ergebnis/baeume_<kachel>.csv` mit jedem einzelnen RCT-Baum.

## Lesart
- RCT ≈ eigene Laserzählung: die 1,84 sind belastbar, RCT ist dann vor allem für Einzelbaum-Maße nützlich.
- RCT deutlich höher: die eigene Laserzählung unterschätzt auch noch (Zwischenstand, Unterstand). Faktor neu ansetzen.
- Getestet nur mit künstlichen Daten (15.625 Kegelbäume: RCT-Auswertung zählt exakt 15.625, Kronendurchmesser stimmt).
