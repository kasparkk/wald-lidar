"""Compare detected trees with reference trees (field inventory, or segmented trees from
3Dtrees.earth / BioDiv-3DTrees). Reference CSV needs columns x, y and optionally height_m, dbh_cm."""
import csv
import numpy as np
from scipy.optimize import linear_sum_assignment


def load_reference(path):
    rows = list(csv.DictReader(open(path)))
    f = lambda r, k: float(r[k]) if r.get(k) not in (None, '') else np.nan
    return np.array([[f(r, 'x'), f(r, 'y'), f(r, 'height_m'), f(r, 'dbh_cm')] for r in rows])


def compare(rows, ref, max_dist=2.0):
    if not len(rows) or not len(ref):
        return dict(note='nothing to compare')
    det = np.array([[r['x'], r['y'], r['height_m'], r['dbh_cm'] if r.get('dbh_source') == 'stem circle fit' else np.nan]
                    for r in rows])
    D = np.hypot(det[:, None, 0] - ref[None, :, 0], det[:, None, 1] - ref[None, :, 1])
    C = np.where(D <= max_dist, D, 1e6)
    di, ri = linear_sum_assignment(C)
    ok = C[di, ri] < 1e6
    di, ri = di[ok], ri[ok]
    tp = len(di)
    out = dict(reference_trees=len(ref), detected_trees=len(det), matched=tp,
               found_pct=round(100 * tp / len(ref), 1),               # recall
               correct_pct=round(100 * tp / len(det), 1),             # precision
               position_error_m=round(float(np.mean(D[di, ri])), 2) if tp else None)
    dh = det[di, 2] - ref[ri, 2]; dh = dh[np.isfinite(dh)]
    if len(dh):
        out.update(height_bias_m=round(float(dh.mean()), 2), height_rmse_m=round(float(np.sqrt((dh ** 2).mean())), 2))
    dd = det[di, 3] - ref[ri, 3]; dd = dd[np.isfinite(dd)]
    if len(dd):
        out.update(dbh_compared=len(dd), dbh_bias_cm=round(float(dd.mean()), 1),
                   dbh_rmse_cm=round(float(np.sqrt((dd ** 2).mean())), 1))
    return out
