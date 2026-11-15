# Wald Brandenburg aus dem nDOM

Ein Skript, das für jede 1-km-Kachel Brandenburgs das offene Höhenmodell der LGB herunterlädt,
Einzelbäume und Bestandeswerte berechnet und die Rohdatei gleich wieder löscht.

## Einrichten (einmal)

    pip install numpy scipy matplotlib "laspy[lazrs]"

## Starten

Erst ein kleiner Test, das Gebiet um Löpten/Hammer (49 Kacheln, ca. 10 Minuten):

    python bb_wald.py --out ergebnisse --bbox 405 5774 412 5781

Dann ganz Brandenburg:

    python bb_wald.py --out ergebnisse

Abbrechen mit Strg+C ist jederzeit möglich. Beim nächsten Start werden fertige Kacheln übersprungen.

## Was es braucht

| | ganz Brandenburg (ca. 30.500 Kacheln) |
|---|---|
| Download | ca. 450 GB, wird nach jeder Kachel gelöscht |
| Speicher für Ergebnisse | ca. 3 GB |
| Arbeitsspeicher | ca. 1,2 GB pro parallelem Prozess |
| Dauer | grob 1 Tag bei 4 Prozessen und 100 Mbit/s |

Gemessen an einer Testkachel: 8 bis 11 Sekunden Rechenzeit pro Kachel.

## Was herauskommt

- `kacheln.csv`: pro Kachel Waldfläche, Zahl der Bäume, Median-Bestandeshöhe, Kronenschluss, grober Vorrat
- `baeume/<Ost>/<Kachel>.csv.gz`: jeder erkannte Baum mit UTM-33-Koordinaten und Höhe
- `uebersicht.png`: Karte der Bestandeshöhen für ganz Brandenburg

## Was die Zahlen bedeuten

- Lage und Höhe der Bäume sind gemessen.
- Die Baumzahl ist eher unterschätzt: dicht stehende Kronen verschmelzen.
- Der Vorrat ist ein unkalibrierter Faustwert (14 m³/ha je Meter mittlerer Höhe × Kronenschluss).
  Mit Forsteinrichtungs- oder Harvesterdaten lässt er sich später kalibrieren.

Daten: GeoBasis-DE / LGB, Datenlizenz Deutschland Namensnennung 2.0
