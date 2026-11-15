#!/usr/bin/env python3
"""
bb_laser_kalib.py – Wie viele Bäume übersieht das Luftbild-Höhenmodell (nDOM)?

Nimmt eine Zufallsstichprobe von Waldkacheln, für die es offene Laser-Punktwolken der LGB gibt,
und zählt auf jeder Kachel die Bäume zweimal: einmal aus dem nDOM (wie bb_wald.py) und einmal aus dem Laser.
Ergebnis: Korrekturfaktoren nach Bestandeshöhe und Kronenschluss, die dann für ganz Brandenburg gelten.

Pro Kachel: Laser ca. 120 MB + nDOM ca. 15 MB laden -> vergleichen -> beide Dateien löschen.

  python bb_laser_kalib.py --out ergebnisse --n 50
  python bb_laser_kalib.py --out ergebnisse --auswerten        # nur Faktoren neu berechnen

Daten: GeoBasis-DE / LGB, dl-de/by-2-0
"""
import argparse, csv, glob, json, os, random, re, sys, tempfile, time, urllib.request, warnings, zipfile
from multiprocessing import Pool
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)

ALS = "https://data.geobasis-bb.de/geobasis/daten/als/laz/"
NDOM = "https://data.geobasis-bb.de/geobasis/daten/ndom/"
HCLS = [(6, 15), (15, 20), (20, 25), (25, 60)]           # Bestandeshöhe (p95) in m
CCLS = [(0.5, 0.7), (0.7, 0.85), (0.85, 1.01)]           # Kronenschluss


def get(url, path):
    for a in range(4):
        try:
            urllib.request.urlretrieve(url, path); return path
        except Exception:
            if a == 3: raise
            time.sleep(5 * (a + 1))


def listing(url, pat):
    html = urllib.request.urlopen(url, timeout=180).read().decode("utf-8", "ignore")
    return sorted({(int(e), int(n)) for e, n in re.findall(pat, html)})


def detect(chm):
    from scipy import ndimage as nd
    s = nd.gaussian_filter(chm, 1.5)
    pk = np.zeros_like(s, bool)
    for lo, hi, w in [(6, 15, 5), (15, 22, 7), (22, 60, 9)]:
        pk |= (s == nd.maximum_filter(s, size=w)) & (s > lo) & (s <= hi)
    return np.nonzero(pk)


def chm_ndom(path, e, n):
    import laspy
    f = laspy.read(path)
    c = np.clip(((np.asarray(f.x) - e * 1000) / 0.2).astype(np.int32), 0, 4999)
    r = np.clip((((n + 1) * 1000 - np.asarray(f.y)) / 0.2).astype(np.int32), 0, 4999)
    g = np.zeros((5000, 5000), np.float32); g[r, c] = np.clip(np.asarray(f.z), 0, None)
    return g.reshape(2500, 2, 2500, 2).max(axis=(1, 3))


def chm_laser(path, e, n):
    """Kronenhöhenmodell 0,4 m aus der Punktwolke; Boden = Klasse 2. Dazu: Anteil Bodentreffer je 20-m-Zelle."""
    import laspy
    from scipy import ndimage as nd
    x0, y1 = e * 1000, (n + 1) * 1000
    f = laspy.read(path)
    X = np.asarray(f.x); Y = np.asarray(f.y); Z = np.asarray(f.z, dtype=np.float32); C = np.asarray(f.classification)
    del f
    m = (X >= x0) & (X < x0 + 1000) & (Y > y1 - 1000) & (Y <= y1) & (C != 7) & (C != 18)
    X, Y, Z, C = X[m], Y[m], Z[m], C[m]
    gc = (X - x0).astype(np.int32).clip(0, 999); gr = (y1 - Y).astype(np.int32).clip(0, 999)
    dgm = np.full((1000, 1000), np.inf, np.float32)
    b = C == 2
    if b.sum() < 1000:
        raise RuntimeError("zu wenige Bodenpunkte")
    np.minimum.at(dgm, (gr[b], gc[b]), Z[b])
    idx = nd.distance_transform_edt(~np.isfinite(dgm), return_distances=False, return_indices=True)
    dgm = dgm[idx[0], idx[1]]
    H = Z - dgm[gr, gc]
    del Z
    ok = (H > -1) & (H < 60)
    cc = ((X[ok] - x0) / 0.4).astype(np.int32).clip(0, 2499); rr = ((y1 - Y[ok]) / 0.4).astype(np.int32).clip(0, 2499)
    chm = np.zeros((2500, 2500), np.float32)
    np.maximum.at(chm, (rr, cc), np.clip(H[ok], 0, None))
    chm = np.where(chm == 0, nd.maximum_filter(chm, 3), chm)
    # Durchdringung: Anteil Punkte unter 1 m an allen Punkten je 20-m-Zelle
    cell = (rr // 50) * 50 + cc // 50
    tot = np.bincount(cell, minlength=2500); low = np.bincount(cell, weights=(H[ok] < 1), minlength=2500)
    return chm, (low / np.maximum(tot, 1)).reshape(50, 50)


def als_date(meta):
    m = re.search(r"Laserscanbefliegung:\s*(\d{4}-\d{2}-\d{2})", re.sub(r"<[^>]+>", " ", meta))
    return m.group(1) if m else ""


def cells(chm):
    b = chm.reshape(50, 50, 50, 50)
    return np.percentile(b, 95, axis=(1, 3)), (b > 5).mean(axis=(1, 3))


def job(args):
    e, n, out = args
    key = f"33{e:03d}-{n:04d}"
    target = os.path.join(out, "laser_kalib", key + ".csv")
    if os.path.exists(target) or os.path.exists(target + ".kein_wald"):
        return "skip"
    with tempfile.TemporaryDirectory() as tmp:
        try:
            nd_path = get(NDOM + f"NDOM_{key}_LAS12.laz", os.path.join(tmp, "ndom.laz"))
            cb = chm_ndom(nd_path, e, n); os.remove(nd_path)
            pb, cvb = cells(cb)
            if ((pb > 10) & (cvb > 0.5)).mean() < 0.4:            # weniger als 40 % Wald -> überspringen
                os.makedirs(os.path.dirname(target), exist_ok=True); open(target + ".kein_wald", "w").close()
                return "kein Wald"
            z = get(ALS + f"als_{key}.zip", os.path.join(tmp, "als.zip"))
            with zipfile.ZipFile(z) as zf:
                lazn = [x for x in zf.namelist() if x.endswith(".laz")][0]
                meta = "".join(zf.read(x).decode("latin-1") for x in zf.namelist() if x.endswith(".html"))
                zf.extract(lazn, tmp)
            os.remove(z)
            cl, pen = chm_laser(os.path.join(tmp, lazn), e, n)
        except Exception as ex:
            return f"Fehler {key}: {ex}"
    pl, cvl = cells(cl)
    wald = (pb > 10) & (pl > 10) & (cvb > 0.5) & (cvl > 0.5)
    stabil = wald & (np.abs(pb - pl) < 4) & (np.abs(cvb - cvl) < 0.2)     # Bestand zwischen Laser- und Luftbildflug unverändert
    cnt = {}
    for name, chm in (("luftbild", cb), ("laser", cl)):
        rr, cc = detect(chm)
        cnt[name] = np.bincount((rr // 50) * 50 + cc // 50, minlength=2500).reshape(50, 50)
    d = als_date(meta)
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kachel", "zelle", "laser_datum", "p95_luftbild", "kronenschluss_luftbild", "durchdringung_laser",
                    "baeume_luftbild", "baeume_laser"])
        for r, c in zip(*np.nonzero(stabil)):
            w.writerow([key, r * 50 + c, d, f"{pb[r, c]:.1f}", f"{cvb[r, c]:.2f}", f"{pen[r, c]:.3f}",
                        int(cnt["luftbild"][r, c]), int(cnt["laser"][r, c])])
    return f"ok {key}: {int(stabil.sum())} Zellen, Luftbild {int(cnt['luftbild'][stabil].sum())} / Laser {int(cnt['laser'][stabil].sum())} Bäume"


def auswerten(out):
    rows = []
    for p in glob.glob(os.path.join(out, "laser_kalib", "*.csv")):
        with open(p) as fh:
            rows += list(csv.DictReader(fh))
    if not rows:
        print("Noch keine Ergebnisse."); return
    P = np.array([float(r["p95_luftbild"]) for r in rows]); C = np.array([float(r["kronenschluss_luftbild"]) for r in rows])
    D = np.array([float(r["durchdringung_laser"]) for r in rows])
    B = np.array([int(r["baeume_luftbild"]) for r in rows]); L = np.array([int(r["baeume_laser"]) for r in rows])
    K = np.array([r["kachel"] for r in rows])
    tab = []
    for h0, h1 in HCLS:
        for c0, c1 in CCLS:
            m = (P >= h0) & (P < h1) & (C >= c0) & (C < c1)
            nk = len(set(K[m]))
            if m.sum() < 25 or nk < 3:
                continue
            # Faktor als Verhältnis der Summen; Streuung über Kacheln (wie stark schwankt es von Ort zu Ort)
            per = [L[m & (K == k)].sum() / max(1, B[m & (K == k)].sum()) for k in set(K[m])]
            tab.append({"hoehe_m": [h0, h1], "kronenschluss": [c0, c1], "zellen": int(m.sum()), "kacheln": nk,
                        "faktor": round(L[m].sum() / max(1, B[m].sum()), 2),
                        "faktor_streuung_kacheln": [round(float(np.percentile(per, 10)), 2), round(float(np.percentile(per, 90)), 2)],
                        "luftbild_je_ha": round(B[m].mean() / 0.04), "laser_je_ha": round(L[m].mean() / 0.04)})
    res = {"kacheln": len(set(K)), "zellen": len(rows), "flaeche_ha": round(len(rows) * 0.04, 1),
           "faktor_gesamt": round(L.sum() / max(1, B.sum()), 2),
           "laser_datum": sorted({r["laser_datum"] for r in rows}),
           "durchdringung_hinweis": "Anteil Laserpunkte unter 1 m je Zelle (wie offen der Bestand für den Laser war)",
           "faktor_nach_durchdringung": {f"{a}-{b}": round(L[(D >= a) & (D < b)].sum() / max(1, B[(D >= a) & (D < b)].sum()), 2)
                                          for a, b in [(0, .2), (.2, .4), (.4, .6), (.6, 1.01)] if ((D >= a) & (D < b)).sum() >= 25},
           "tabelle": tab}
    with open(os.path.join(out, "laser_kalib", "kalibrierung.json"), "w") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print(f"\n{res['kacheln']} Kacheln, {res['flaeche_ha']} ha stabiler Wald. Laser findet {res['faktor_gesamt']}-mal so viele Bäume wie das Luftbild.")
    print("Bestandeshöhe  Kronenschluss  Faktor  (10-90 % der Kacheln)  Luftbild/ha  Laser/ha")
    for t in tab:
        print(f"  {t['hoehe_m'][0]:>2}-{t['hoehe_m'][1]:<2} m      {t['kronenschluss'][0]:.2f}-{min(1, t['kronenschluss'][1]):.2f}     {t['faktor']:>4}   "
              f"({t['faktor_streuung_kacheln'][0]}-{t['faktor_streuung_kacheln'][1]})            {t['luftbild_je_ha']:>4}      {t['laser_je_ha']:>4}")
    print("-> " + os.path.join(out, "laser_kalib", "kalibrierung.json"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="ergebnisse")
    ap.add_argument("--n", type=int, default=50, help="so viele Waldkacheln mit Laserdaten auswerten")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--auswerten", action="store_true")
    a = ap.parse_args()
    if not a.auswerten:
        print("Lade Kachellisten ...", flush=True)
        als = listing(ALS, r"als_33(\d{3})-(\d{4})\.zip")
        ndom = set(listing(NDOM, r"NDOM_33(\d{3})-(\d{4})_LAS12\.laz"))
        cand = [t for t in als if t in ndom]
        random.Random(42).shuffle(cand)            # feste Zufallsreihenfolge: Neustart setzt an derselben Stelle fort
        print(f"{len(cand)} Kacheln mit Laser und Luftbild. Ziel: {a.n} Waldkacheln.", flush=True)
        d = os.path.join(a.out, "laser_kalib"); os.makedirs(d, exist_ok=True)
        fertig = lambda: len(glob.glob(os.path.join(d, "*.csv")))
        i = 0
        with Pool(a.workers) as pool:
            while fertig() < a.n and i < len(cand):
                batch = cand[i:i + a.workers * 2]; i += len(batch)
                for res in pool.imap_unordered(job, [(e, n, a.out) for e, n in batch]):
                    if res != "skip":
                        print(f"[{fertig()}/{a.n}] {res}", flush=True)
    auswerten(a.out)


if __name__ == "__main__":
    main()
