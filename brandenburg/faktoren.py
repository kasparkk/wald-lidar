"""Laser-Korrektur der Baumzahl (Stichprobe 50 Kacheln) und Vorrats-Angleichung an die BWI 2022."""
import json
K = json.load(open(__import__('os').path.join(__import__('os').path.dirname(__file__),'results','kalibrierung.json'), encoding='latin-1'))
TAB = K['tabelle']; GESAMT = K['faktor_gesamt']
VORRAT_FAKTOR = round(296 / 214, 3)   # BWI 2022 Brandenburg 296 Vfm/ha  /  unser roher Landesmittelwert 214 m³/ha

def faktor(hoehe, schluss):
    h = min(max(hoehe, 6.01), 59.9); c = min(max(schluss, 0.5), 1.0)
    for t in TAB:
        if t['hoehe_m'][0] <= h < t['hoehe_m'][1] and t['kronenschluss'][0] <= c < t['kronenschluss'][1]:
            return t['faktor']
    return GESAMT

def korrigiere(t):
    f = faktor(t['hoehe_median_m'], t['kronenschluss']) if t['wald_ha'] > 0 else 1.0
    t['faktor_laser'] = f
    t['baeume_geschaetzt'] = int(round(t['baeume'] * f))
    t['vorrat_m3'] = int(round(t['vorrat_grob_m3'] * VORRAT_FAKTOR))
    return t
