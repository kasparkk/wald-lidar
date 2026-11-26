"""
Vergleicht die Einzelbäume von RayCloudTools (3Dtrees.earth) mit den eigenen zwei Zählungen:
  ndom   = Kronenmaxima im nDOM der LGB (wie bb_wald.py)
  laser  = Kronenmaxima im Kronenmodell aus der Laser-Punktwolke (wie bb_laser_kalib.py)
  rct    = von RayCloudTools segmentierte Bäume (Extra-Dimension PredInstance)

Eingabe:
  ergebnis/*.laz   die von 3Dtrees zurückgegebenen Dateien (Kachelname im Dateinamen oder wird aus der Lage erkannt)
  upload/, ndom/   von rct_holen.py

Ausgabe:
  ergebnis/baeume_<kachel>.csv   jeder RCT-Baum: Lage, Höhe, Kronenfläche, Punktzahl
  ergebnis/vergleich.csv         eine Zeile je Kachel
Gezählt wird nur im Wald (20-m-Zellen mit p95 > 10 m und Kronenschluss > 0,5 im nDOM), nur Bäume ab 6 m.
"""
import csv, glob, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "brandenburg"))
from bb_laser_kalib import detect, chm_ndom, chm_laser, cells  # noqa: E402

MIN_H = 6.0       # wie detect(): Kronen unter 6 m zählen nicht
MATCH_M = 2.5     # zwei Baumspitzen gelten als derselbe Baum, wenn sie so nah beieinander liegen


def instance_dim(f):
    names = list(f.point_format.dimension_names)
    for n in names:
        if n.lower() == "predinstance":
            return n
    for n in names:
        if "instance" in n.lower():
            return n
    raise RuntimeError("keine Instanz-Spalte gefunden, vorhanden: " + ", ".join(names))


def dgm(X, Y, Z, C, e, n):
    from scipy import ndimage as nd
    x0, y1 = e * 1000, (n + 1) * 1000
    b = (C == 2) & (X >= x0) & (X < x0 + 1000) & (Y > y1 - 1000) & (Y <= y1)
    if b.sum() < 1000:
        return None
    g = np.full((1000, 1000), np.inf, np.float32)
    np.minimum.at(g, ((y1 - Y[b]).astype(np.int32).clip(0, 999), (X[b] - x0).astype(np.int32).clip(0, 999)), Z[b])
    idx = nd.distance_transform_edt(~np.isfinite(g), return_distances=False, return_indices=True)
    return g[idx[0], idx[1]]


def peaks_xy(chm, e, n):
    rr, cc = detect(chm)
    return np.column_stack([e * 1000 + (cc + 0.5) * 0.4, (n + 1) * 1000 - (rr + 0.5) * 0.4])


def in_forest(xy, mask, e, n):
    c = ((xy[:, 0] - e * 1000) // 20).astype(int).clip(0, 49)
    r = ((((n + 1) * 1000 - xy[:, 1]) // 20)).astype(int).clip(0, 49)
    inside = (xy[:, 0] >= e * 1000) & (xy[:, 0] < e * 1000 + 1000) & (xy[:, 1] > n * 1000) & (xy[:, 1] <= n * 1000 + 1000)
    return inside & mask[r, c]


def matched_share(a, b):
    """Anteil der Punkte in a, die in b einen Partner innerhalb MATCH_M haben."""
    from scipy.spatial import cKDTree
    if len(a) == 0 or len(b) == 0:
        return 0.0
    d, _ = cKDTree(b).query(a, distance_upper_bound=MATCH_M)
    return float(np.isfinite(d).mean())


def rct_trees(path, key, e, n):
    import laspy
    from scipy.spatial import ConvexHull
    f = laspy.read(path)
    X = np.asarray(f.x); Y = np.asarray(f.y); Z = np.asarray(f.z, dtype=np.float32)
    I = np.asarray(f[instance_dim(f)]).astype(np.int64)
    C = np.asarray(f.classification) if "classification" in f.point_format.dimension_names else np.zeros(len(X), np.uint8)
    del f
    g = dgm(X, Y, Z, C, e, n)
    if g is None:   # Bodenklasse ging verloren -> Boden aus der Originaldatei
        import laspy as lp
        o = lp.read(os.path.join(HERE, "upload", f"als_{key}.laz"))
        g = dgm(np.asarray(o.x), np.asarray(o.y), np.asarray(o.z, dtype=np.float32), np.asarray(o.classification), e, n)
        del o
    x0, y1 = e * 1000, (n + 1) * 1000
    H = Z - g[(y1 - Y).astype(np.int32).clip(0, 999), (X - x0).astype(np.int32).clip(0, 999)]
    keep = I > 0     # 0 oder -1 = kein Baum
    X, Y, H, I = X[keep], Y[keep], H[keep], I[keep]
    order = np.argsort(I, kind="stable")
    X, Y, H, I = X[order], Y[order], H[order], I[order]
    ids, start = np.unique(I, return_index=True)
    stop = np.r_[start[1:], len(I)]
    out = []
    for t, s, z in zip(ids, start, stop):
        h = H[s:z]
        if len(h) < 10:
            continue
        top = s + int(np.argmax(h))
        crown = (h > 2.0)
        area = 0.0
        if crown.sum() >= 3:
            try:
                area = ConvexHull(np.column_stack([X[s:z][crown], Y[s:z][crown]])).volume
            except Exception:
                area = 0.0
        out.append((int(t), float(X[top]), float(Y[top]), float(h.max()), area, int(z - s)))
    return out


def key_of(path):
    m = re.search(r"(33\d{3})-(\d{4})", os.path.basename(path))
    if m:
        return f"{m.group(1)}-{m.group(2)}"
    import laspy
    with laspy.open(path) as r:
        h = r.header
        return f"33{int(h.mins[0] // 1000):03d}-{int(h.mins[1] // 1000):04d}"


def main():
    res_dir = os.path.join(HERE, "ergebnis")
    files = sorted(p for p in glob.glob(os.path.join(res_dir, "*.laz")) + glob.glob(os.path.join(res_dir, "*.las")))
    if not files:
        print("Keine Dateien in", res_dir, "- bitte die Ergebnisse von 3Dtrees dort ablegen.")
        return 1
    rows = []
    for p in files:
        key = key_of(p)
        e, n = int(key[2:5]), int(key[6:])
        print(f"{key}: lese RCT-Ergebnis ...", flush=True)
        trees = rct_trees(p, key, e, n)
        with open(os.path.join(res_dir, f"baeume_{key}.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["id", "x", "y", "hoehe_m", "kronenflaeche_m2", "punkte"])
            for t in trees:
                w.writerow([t[0], round(t[1], 2), round(t[2], 2), round(t[3], 1), round(t[4], 1), t[5]])
        print(f"{key}: eigene Zählungen ...", flush=True)
        cb = chm_ndom(os.path.join(HERE, "ndom", f"NDOM_{key}.laz"), e, n)
        p95, cov = cells(cb)
        wald = (p95 > 10) & (cov > 0.5)
        nd_xy = peaks_xy(cb, e, n)
        cl, _ = chm_laser(os.path.join(HERE, "upload", f"als_{key}.laz"), e, n)
        la_xy = peaks_xy(cl, e, n)
        T = np.array([[t[1], t[2], t[3], t[4]] for t in trees]) if trees else np.zeros((0, 4))
        T = T[T[:, 2] >= MIN_H]
        rc_xy = T[:, :2]
        fm = in_forest(rc_xy, wald, e, n)
        rc_xy, T = rc_xy[fm], T[fm]
        nd_xy = nd_xy[in_forest(nd_xy, wald, e, n)]
        la_xy = la_xy[in_forest(la_xy, wald, e, n)]
        ha = wald.sum() * 0.04
        r = {
            "kachel": key, "wald_ha": round(ha, 1),
            "rct": len(rc_xy), "laser": len(la_xy), "ndom": len(nd_xy),
            "rct_je_ha": round(len(rc_xy) / max(ha, 1e-9)), "laser_je_ha": round(len(la_xy) / max(ha, 1e-9)),
            "ndom_je_ha": round(len(nd_xy) / max(ha, 1e-9)),
            "faktor_rct_zu_ndom": round(len(rc_xy) / max(len(nd_xy), 1), 2),
            "faktor_laser_zu_ndom": round(len(la_xy) / max(len(nd_xy), 1), 2),
            "ndom_in_rct_gefunden": round(matched_share(nd_xy, rc_xy), 2),
            "rct_in_ndom_gefunden": round(matched_share(rc_xy, nd_xy), 2),
            "laser_in_rct_gefunden": round(matched_share(la_xy, rc_xy), 2),
            "rct_hoehe_median_m": round(float(np.median(T[:, 2])), 1) if len(T) else "",
            "rct_krone_d_median_m": round(float(np.median(2 * np.sqrt(T[:, 3] / np.pi))), 1) if len(T) else "",
        }
        rows.append(r)
        print(f"  Wald {r['wald_ha']} ha | Bäume je ha: RCT {r['rct_je_ha']}, Laser {r['laser_je_ha']}, nDOM {r['ndom_je_ha']}"
              f" | RCT/nDOM {r['faktor_rct_zu_ndom']}, Laser/nDOM {r['faktor_laser_zu_ndom']}", flush=True)
    with open(os.path.join(res_dir, "vergleich.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("\nFertig:", os.path.join(res_dir, "vergleich.csv"))
    print("Zum Vergleich: eigene Laser-Kalibrierung (50 Kacheln) ergab Faktor 1,84 Laser/nDOM.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
