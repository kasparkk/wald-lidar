#!/usr/bin/env python3
"""
bb_extra.py – Baumart (Nadel/Laub) und Vitalität für jeden erkannten Baum aus den
Infrarot-Luftbildern (DOP20 RGBI) der LGB Brandenburg.

Läuft auf Kacheln, die bb_wald.py schon gerechnet hat. Pro Kachel:
Luftbild laden (ca. 86 MB) -> Kronenfarben messen -> Luftbild löschen.
Danach Nadel/Laub-Trennung je Bildflugtag und neue Baumdateien schreiben.

Daten: GeoBasis-DE / LGB, dl-de/by-2-0 – https://data.geobasis-bb.de/geobasis/daten/dop/rgbi_tif/

  python bb_extra.py --out ergebnisse --bbox 405 5774 412 5781
"""
import argparse, csv, glob, gzip, io, json, os, re, sys, tempfile, time, urllib.request, warnings, zipfile
from multiprocessing import Pool
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)

DOP = "https://data.geobasis-bb.de/geobasis/daten/dop/rgbi_tif/dop_33{e:03d}-{n:04d}.zip"


def tree_file(out, e, n):
    return os.path.join(out, "baeume", str(e), f"33{e:03d}-{n:04d}.csv.gz")


def read_trees(path):
    with gzip.open(path, "rt") as fh:
        rows = list(csv.reader(fh))[1:]
    a = np.array([[float(r[0]), float(r[1]), float(r[2])] for r in rows]) if rows else np.zeros((0, 3))
    return a[:, 0], a[:, 1], a[:, 2]


def fetch_dop(e, n, folder, local=None):
    if local:
        tif = glob.glob(os.path.join(local, f"dop_33{e:03d}-{n:04d}*.tif"))
        html = glob.glob(os.path.join(local, f"dop_33{e:03d}-{n:04d}*.html"))
        return tif[0], (open(html[0], encoding="latin-1").read() if html else "")
    zpath = os.path.join(folder, "dop.zip")
    for attempt in range(4):
        try:
            urllib.request.urlretrieve(DOP.format(e=e, n=n), zpath); break
        except Exception:
            if attempt == 3: raise
            time.sleep(5 * (attempt + 1))
    with zipfile.ZipFile(zpath) as z:
        z.extractall(folder)
        meta = "".join(z.read(m).decode("latin-1") for m in z.namelist() if m.endswith(".html"))
    os.remove(zpath)
    return glob.glob(os.path.join(folder, "*.tif"))[0], meta


def features(args):
    e, n, out, local = args
    target = os.path.join(out, "extra", f"33{e:03d}-{n:04d}.npz")
    if os.path.exists(target):
        return "skip"
    tf = tree_file(out, e, n)
    if not os.path.exists(tf):
        return "keine Bäume"
    import tifffile
    from scipy import ndimage as nd
    from scipy.spatial import cKDTree
    E, N, H = read_trees(tf)
    with tempfile.TemporaryDirectory() as tmp:
        try:
            tif, meta = fetch_dop(e, n, tmp, local)
        except Exception as ex:
            return f"Fehler Luftbild 33{e:03d}-{n:04d}: {ex}"
        img = tifffile.imread(tif)
    m = re.search(r"Bildflugdatum:\s*(\d{4}-\d{2}-\d{2})", re.sub(r"<[^>]+>", " ", meta))
    datum = m.group(1) if m else "unbekannt"
    R = img[..., 0].astype(np.float32); G = img[..., 1].astype(np.float32); NIR = img[..., 3].astype(np.float32)
    del img
    ndvi = (NIR - R) / (NIR + R + 1e-6)
    r = np.clip(((n + 1) * 1000 - N) / 0.2, 0, 4999).astype(int)
    c = np.clip((E - e * 1000) / 0.2, 0, 4999).astype(int)
    f_nir = nd.uniform_filter(NIR, 7)[r, c]; f_g = nd.uniform_filter(G, 7)[r, c]
    f_ndvi7 = nd.uniform_filter(ndvi, 7)[r, c]; f_ndvi5 = nd.uniform_filter(ndvi, 5)[r, c]
    os.makedirs(os.path.dirname(target), exist_ok=True)
    np.savez_compressed(target, E=E, N=N, H=H, nir=f_nir.astype(np.float32), g=f_g.astype(np.float32),
                        ndvi7=f_ndvi7.astype(np.float32), ndvi5=f_ndvi5.astype(np.float32), datum=datum)
    return "ok"


def classify_and_write(out, keys):
    from sklearn.mixture import GaussianMixture
    from scipy.spatial import cKDTree
    data = {}
    for k in keys:
        p = os.path.join(out, "extra", k + ".npz")
        if os.path.exists(p):
            d = np.load(p); data[k] = {x: d[x] for x in d.files}; data[k]["datum"] = str(d["datum"])
    by_date = {}
    for k, d in data.items():
        by_date.setdefault(d["datum"], []).append(k)
    report = {}
    for datum, ks in by_date.items():
        X = np.concatenate([np.c_[data[k]["nir"], data[k]["ndvi7"] * 300, data[k]["g"]] for k in ks])
        rng = np.random.default_rng(0)
        Xs = X[rng.choice(len(X), min(len(X), 60000), replace=False)] if len(X) else X
        gm = GaussianMixture(2, random_state=0).fit(Xs)
        laub = int(np.argmax(gm.means_[:, 0]))
        nad = 1 - laub
        # Plausibilität: Laub-Cluster muss deutlich höheres NDVI haben, sonst alles Nadel
        sep = (gm.means_[laub, 1] - gm.means_[nad, 1]) / 300
        ok = sep > 0.06 and gm.weights_[laub] > 0.01
        report[datum] = {"kacheln": len(ks), "laub_ndvi_abstand": round(float(sep), 3), "trennung_ok": bool(ok)}
        for k in ks:
            d = data[k]
            X = np.c_[d["nir"], d["ndvi7"] * 300, d["g"]]
            if len(X) == 0:
                continue
            pr = gm.predict_proba(X)
            isL = (pr.argmax(1) == laub) if ok else np.zeros(len(X), bool)
            sicher = np.round(100 * pr.max(1)).astype(int) if ok else np.full(len(X), 0)
            # Vitalität: NDVI gegen gleichartige Nachbarn im 40-m-Umkreis
            z = np.zeros(len(X)); pts = np.c_[d["E"], d["N"]]
            for cls in (True, False):
                idx = np.where(isL == cls)[0]
                if len(idx) < 5:
                    continue
                t = cKDTree(pts[idx])
                for j, nb in zip(idx, t.query_ball_point(pts[idx], 40)):
                    ref = d["ndvi5"][idx[nb]]; med = np.median(ref)
                    mad = np.median(np.abs(ref - med)) * 1.4826 + 1e-3
                    z[j] = (d["ndvi5"][j] - med) / mad
            dead = (d["ndvi5"] < 0.12) & (d["H"] > 8)
            weak = (~dead) & (z < -3) & (d["H"] > 8)
            vit = np.where(dead, "t", np.where(weak, "a", "v"))
            e, n = int(k[2:5]), int(k[6:])
            with gzip.open(tree_file(out, e, n), "wt", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["ost_utm33", "nord_utm33", "hoehe_m", "art", "sicherheit_pct", "vital", "ndvi"])
                for i in range(len(X)):
                    w.writerow([f"{d['E'][i]:.1f}", f"{d['N'][i]:.1f}", f"{d['H'][i]:.1f}",
                                ("L" if isL[i] else "N") if ok else "?", int(sicher[i]), vit[i], f"{d['ndvi5'][i]:.2f}"])
            js = os.path.join(out, "kacheln", k + ".json")
            if os.path.exists(js):
                with open(js) as fh: s = json.load(fh)
                s.update({"laub_anteil": round(float(isL.mean()), 3) if ok else None,
                          "auffaellig": int(weak.sum()), "abgestorben": int(dead.sum()), "luftbild": datum})
                with open(js, "w") as fh: json.dump(s, fh)
    with open(os.path.join(out, "extra", "bericht.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="ergebnisse")
    ap.add_argument("--bbox", type=int, nargs=4, required=True)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--local", help="Ordner mit bereits entpackten dop_*.tif (Test)")
    a = ap.parse_args()
    e1, n1, e2, n2 = a.bbox
    keys = []
    for fn in sorted(os.listdir(os.path.join(a.out, "kacheln"))):
        m = re.match(r"33(\d{3})-(\d{4})\.json$", fn)
        if m and e1 <= int(m.group(1)) <= e2 and n1 <= int(m.group(2)) <= n2:
            keys.append((int(m.group(1)), int(m.group(2))))
    print(f"{len(keys)} Kacheln: Luftbilder laden und Kronen messen ...", flush=True)
    t0 = time.time()
    with Pool(a.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(features, [(e, n, a.out, a.local) for e, n in keys]), 1):
            if res not in ("ok", "skip"):
                print(res, file=sys.stderr)
            if i % 5 == 0 or i == len(keys):
                el = time.time() - t0
                print(f"{i}/{len(keys)} · noch ca. {(len(keys)-i)*el/i/60:.0f} min", flush=True)
    print("Nadel/Laub trennen und Vitalität berechnen ...", flush=True)
    rep = classify_and_write(a.out, [f"33{e:03d}-{n:04d}" for e, n in keys])
    for datum, r in rep.items():
        print(f"  Bildflug {datum}: {r['kacheln']} Kacheln, Trennung {'ok' if r['trennung_ok'] else 'UNSICHER - alles als ? markiert'}")
    print("Fertig.")


if __name__ == "__main__":
    main()
