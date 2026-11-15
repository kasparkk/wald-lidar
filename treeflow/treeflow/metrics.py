"""Measurements per tree. Every value says how it was obtained, so nobody mistakes an estimate
for a measurement."""
import numpy as np
from scipy.spatial import ConvexHull


def tree_table(c, hag, tree_id, seeds, chm, lab, cfg, core_bbox=None):
    res = chm['res']
    nmax = max([s['id'] for s in seeds] + [int(lab.max()), int(tree_id.max()), 0])
    areas = np.bincount(lab.ravel(), minlength=nmax + 1) * res * res
    order = np.argsort(tree_id, kind='stable')
    tid_sorted = tree_id[order]
    starts = np.searchsorted(tid_sorted, np.arange(nmax + 2))
    rows = []
    for s in seeds:
        k = s['id']
        idx = order[starts[k]:starts[k + 1]]
        if len(idx) < cfg.get('min_points_per_tree', 15):
            continue
        h = hag[idx]
        stem = s.get('stem')
        px, py = (stem['x'], stem['y']) if stem else (s['x'], s['y'])
        if core_bbox is not None:
            x0, y0, x1, y1 = core_bbox
            if not (x0 <= px < x1 and y0 <= py < y1):
                continue           # tree belongs to the neighbouring tile
        height = float(np.percentile(h, 99.5))
        # crown base: lowest of the points that are clearly away from the stem axis
        dxy = np.hypot(c.x[idx] - px, c.y[idx] - py)
        away = dxy > max(0.6, 3 * (stem['r'] if stem else 0.15))
        crown_pts = idx[away & (h > 1.5)]
        crown_base = float(np.percentile(hag[crown_pts], 5)) if len(crown_pts) >= 10 else None
        area = float(areas[k])
        vol = None
        if len(crown_pts) >= 20:
            try:
                P = np.c_[c.x[crown_pts], c.y[crown_pts], hag[crown_pts]]
                P = P[hag[crown_pts] >= (crown_base or 0)]
                vol = float(ConvexHull(P).volume)
            except Exception:
                vol = None
        dbh = round(2 * stem['r'] * 100, 1) if stem else None
        rows.append(dict(
            tree=k, x=round(px, 2), y=round(py, 2),
            height_m=round(height, 2),
            crown_base_m=None if crown_base is None else round(crown_base, 2),
            crown_area_m2=round(area, 1),
            crown_diameter_m=round(2 * np.sqrt(area / np.pi), 2),
            crown_volume_m3=None if vol is None else round(vol, 1),
            dbh_cm=dbh,
            dbh_source='stem circle fit' if dbh else None,
            stem_fit_rmse_cm=round(stem['rmse'] * 100, 2) if stem else None,
            points=int(len(idx)),
            later_returns_pct=round(100 * float((c.ret[idx] > 1).mean()), 1),
        ))
    return rows


def add_dbh_estimate(rows, a=1.5, b=-6.0):
    """Optional and clearly labelled: DBH from height with a simple linear rule (pine-like default).
    Only fills trees without a measured stem."""
    for r in rows:
        if r['dbh_cm'] is None:
            r['dbh_cm'] = round(max(5.0, a * r['height_m'] + b), 1)
            r['dbh_source'] = f'estimate from height ({a}*h{b:+})'
    return rows
