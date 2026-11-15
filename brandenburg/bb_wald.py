#!/usr/bin/env python3
"""
bb_wald.py – Einzelbäume und Bestandeswerte für ganz Brandenburg aus dem offenen nDOM der LGB.

Pro 1-km-Kachel: Download (ca. 15 MB) -> Höhenmodell -> Einzelbaumerkennung -> Kennwerte -> Datei löschen.
Es bleiben nur die Ergebnisse übrig (ca. 150 KB pro Kachel), nicht die Rohdaten.

Daten: GeoBasis-DE / LGB, Datenlizenz Deutschland Namensnennung 2.0
       https://data.geobasis-bb.de/geobasis/daten/ndom/

Benutzung:
  pip install numpy scipy "laspy[lazrs]"
  python bb_wald.py --out ergebnisse --bbox 405 5774 412 5781    # Testgebiet um Löpten/Hammer (km, UTM 33)
  python bb_wald.py --out ergebnisse                               # ganz Brandenburg
  python bb_wald.py --out ergebnisse --local NDOM_33408-5778_LAS12.laz   # eine Datei, die schon da ist
  python bb_wald.py --out ergebnisse --karte                       # Übersichtskarte aus den Ergebnissen

Abbrechen und neu starten ist jederzeit möglich: fertige Kacheln werden übersprungen.
"""
import argparse, csv, gzip, json, os, re, sys, tempfile, time, urllib.request
from multiprocessing import Pool
import numpy as np, warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

BASE = "https://data.geobasis-bb.de/geobasis/daten/ndom/"
TILE_RE = re.compile(r"NDOM_33(\d{3})-(\d{4})_LAS12\.laz")
K_VORRAT = 14.0   # m³/ha je m mittlerer Höhe – Ertragstafel-Faustwert Kiefer, NICHT kalibriert


def tile_list():
    html = urllib.request.urlopen(BASE, timeout=120).read().decode("utf-8", "ignore")
    return sorted({(int(e), int(n)) for e, n in TILE_RE.findall(html)})


def download(e, n, folder):
    name = f"NDOM_33{e:03d}-{n:04d}_LAS12.laz"
    path = os.path.join(folder, name)
    for attempt in range(4):
        try:
            urllib.request.urlretrieve(BASE + name, path)
            return path
        except Exception:
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("Download fehlgeschlagen: " + name)


def process(path, e, n):
    import laspy
    from scipy import ndimage as nd
    f = laspy.read(path)
    x0, y1 = e * 1000, (n + 1) * 1000
    c = np.clip(((np.asarray(f.x) - x0) / 0.2).astype(np.int32), 0, 4999)
    r = np.clip(((y1 - np.asarray(f.y)) / 0.2).astype(np.int32), 0, 4999)
    g = np.zeros((5000, 5000), np.float32)
    g[r, c] = np.clip(np.asarray(f.z), 0, None)
    del f

    # Bestandeswerte auf 20-m-Zellen
    b = g.reshape(50, 100, 50, 100)
    p95 = np.percentile(b, 95, axis=(1, 3))
    cover = (b > 5).mean(axis=(1, 3))
    veg = np.where(b > 5, b, np.nan)
    with np.errstate(all="ignore"):
        hmean = np.nan_to_num(np.nanmean(veg, axis=(1, 3)))
    wald = (p95 > 5) & (cover > 0.3)
    vorrat = np.where(wald, K_VORRAT * hmean * cover, 0)

    # Einzelbäume: lokale Maxima im geglätteten 0,4-m-Höhenmodell, Fenster wächst mit der Höhe
    chm = g.reshape(2500, 2, 2500, 2).max(axis=(1, 3))
    s = nd.gaussian_filter(chm, 1.5)
    pk = np.zeros_like(s, bool)
    for lo, hi, w in [(6, 15, 5), (15, 22, 7), (22, 60, 9)]:
        pk |= (s == nd.maximum_filter(s, size=w)) & (s > lo) & (s <= hi)
    rr, cc = np.nonzero(pk)
    h = chm[rr, cc]
    ex = x0 + cc * 0.4 + 0.2
    ny = y1 - rr * 0.4 - 0.2

    stats = {
        "kachel": f"33{e:03d}-{n:04d}", "e_km": e, "n_km": n,
        "wald_ha": round(float(wald.sum() * 0.04), 2),
        "baeume": int(len(h)),
        "hoehe_median_m": round(float(np.median(p95[wald])), 1) if wald.any() else 0.0,
        "kronenschluss": round(float(cover[wald].mean()), 2) if wald.any() else 0.0,
        "vorrat_grob_m3": int(round(float((vorrat * 0.04).sum()))),
    }
    return stats, ex, ny, h


def write_tile(out, stats, ex, ny, h):
    tdir = os.path.join(out, "baeume", str(stats["e_km"]))
    os.makedirs(tdir, exist_ok=True)
    with gzip.open(os.path.join(tdir, stats["kachel"] + ".csv.gz"), "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ost_utm33", "nord_utm33", "hoehe_m"])
        for a, b_, c_ in zip(ex, ny, h):
            w.writerow([f"{a:.1f}", f"{b_:.1f}", f"{c_:.1f}"])
    sdir = os.path.join(out, "kacheln")
    os.makedirs(sdir, exist_ok=True)
    with open(os.path.join(sdir, stats["kachel"] + ".json"), "w") as fh:
        json.dump(stats, fh)


def done(out, e, n):
    return os.path.exists(os.path.join(out, "kacheln", f"33{e:03d}-{n:04d}.json"))


def job(args):
    e, n, out = args
    if done(out, e, n):
        return "skip"
    with tempfile.TemporaryDirectory() as tmp:
        try:
            p = download(e, n, tmp)
            write_tile(out, *process(p, e, n))
            return "ok"
        except Exception as ex:
            return f"Fehler 33{e:03d}-{n:04d}: {ex}"


def merge(out):
    rows = []
    sdir = os.path.join(out, "kacheln")
    for fn in sorted(os.listdir(sdir)) if os.path.isdir(sdir) else []:
        with open(os.path.join(sdir, fn)) as fh:
            rows.append(json.load(fh))
    if not rows:
        return rows
    with open(os.path.join(out, "kacheln.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return rows


def karte(out, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    es = [r["e_km"] for r in rows]; ns = [r["n_km"] for r in rows]
    e0, n0 = min(es), min(ns)
    grid = np.full((max(ns) - n0 + 1, max(es) - e0 + 1), np.nan)
    for r in rows:
        if r["wald_ha"] > 10:
            grid[max(ns) - r["n_km"], r["e_km"] - e0] = r["hoehe_median_m"]
    fig, ax = plt.subplots(figsize=(9, 10))
    im = ax.imshow(grid, cmap="YlGn", vmin=5, vmax=30,
                   extent=[e0, max(es) + 1, n0, max(ns) + 1])
    fig.colorbar(im, ax=ax, label="Median-Bestandeshöhe [m] (Kacheln mit >10 ha Wald)", shrink=.7)
    wald = sum(r["wald_ha"] for r in rows); baeume = sum(r["baeume"] for r in rows)
    ax.set_title(f"Wald in Brandenburg aus dem nDOM – {len(rows):,} Kacheln, "
                 f"{wald/1000:,.0f} Tsd. ha Wald, {baeume/1e6:,.0f} Mio. Bäume erkannt".replace(",", "."))
    ax.set_xlabel("Ost [km, UTM 33]"); ax.set_ylabel("Nord [km]")
    fig.text(.01, .01, "Daten: GeoBasis-DE / LGB, dl-de/by-2-0. Baumzahl unterschätzt, Vorrat unkalibriert.", fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(out, "uebersicht.png"), dpi=130)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="ergebnisse")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--bbox", type=int, nargs=4, metavar=("E_MIN", "N_MIN", "E_MAX", "N_MAX"),
                    help="nur Kacheln in diesem Bereich (km, UTM 33, z. B. 405 5774 412 5781)")
    ap.add_argument("--limit", type=int, help="höchstens so viele Kacheln (zum Testen)")
    ap.add_argument("--local", help="eine bereits heruntergeladene NDOM_*.laz verarbeiten")
    ap.add_argument("--karte", action="store_true", help="nur Zusammenfassung und Übersichtskarte erzeugen")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    if a.local:
        m = TILE_RE.search(os.path.basename(a.local))
        e, n = int(m.group(1)), int(m.group(2))
        t = time.time(); stats, *rest = process(a.local, e, n)
        write_tile(a.out, stats, *rest)
        print(json.dumps(stats, ensure_ascii=False), f"({time.time()-t:.0f} s)")
    elif not a.karte:
        tiles = tile_list()
        if a.bbox:
            e1, n1, e2, n2 = a.bbox
            tiles = [(e, n) for e, n in tiles if e1 <= e <= e2 and n1 <= n <= n2]
        if a.limit:
            tiles = tiles[:a.limit]
        todo = [(e, n, a.out) for e, n in tiles]
        print(f"{len(todo)} Kacheln, {a.workers} parallel. Abbruch mit Strg+C, Neustart setzt fort.")
        t0 = time.time(); fertig = 0
        with Pool(a.workers) as pool:
            for i, res in enumerate(pool.imap_unordered(job, todo), 1):
                if res == "ok":
                    fertig += 1
                elif res != "skip":
                    print(res, file=sys.stderr)
                if i % 25 == 0 or i == len(todo):
                    rate = fertig / max(1e-9, time.time() - t0)
                    rest = (len(todo) - i) / rate / 3600 if rate else float("nan")
                    print(f"{i}/{len(todo)} Kacheln · noch ca. {rest:.1f} h", flush=True)

    rows = merge(a.out)
    if rows:
        karte(a.out, rows)
        print(f"Fertig: {len(rows)} Kacheln -> {a.out}/kacheln.csv, {a.out}/uebersicht.png, {a.out}/baeume/")


if __name__ == "__main__":
    main()
