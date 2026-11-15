"""DBH from drone (ULS) points of a segmented tree, with plausibility checks.

The drone has only 4-8 points per tree at breast height, so a slice fit fails. Instead: a tilted, tapering
cylinder fitted to all stem points 0.5-4 m above ground (robust soft-L1 loss, only points within 0.8 m of the
densest stem spot), evaluated at 1.3 m. Tested on 97 Schorfheide trees (BioDiv-3DTrees plot SEW04_GD, beech and
pine, leaf-off) against ground-scan DBH: after the checks 82 trees kept, mean error 2.4 cm (pine 2.2, beech 2.7);
spruce plot AEW33_GD (Swabian Alb): 117 of 221 kept, 2.5 cm with stub trimming (3.4 cm without);
of 15 rejected trees 11 were really more than 5 cm off, 13 thinner than 20 cm. Rule 2 in plausible() was set
after looking at failures, so expect slightly worse results on new data.

  import drone_dbh as D
  r = D.drone_dbh(P, ground_z)                     # P: N x 3 points of one tree
  table = D.plausible(table)                       # DataFrame with columns dbh, height, crown_d, ninl, taper, species
"""
import numpy as np, pandas as pd
from scipy.optimize import least_squares

rng = np.random.default_rng(1)


def circ(xy):
    """Algebraic circle fit (Kasa) with 3 rounds of trimming."""
    for _ in range(3):
        A = np.c_[2 * xy, np.ones(len(xy))]; b = (xy ** 2).sum(1)
        (cx, cy, k), *_ = np.linalg.lstsq(A, b, rcond=None); r = np.sqrt(max(k + cx * cx + cy * cy, 1e-9))
        res = np.abs(np.hypot(xy[:, 0] - cx, xy[:, 1] - cy) - r)
        keep = res < max(2.5 * np.median(res), 0.01)
        if keep.all() or keep.sum() < 4: break
        xy = xy[keep]
    return cx, cy, r, xy


def ransac(xy, it=400, tol=0.04):
    best = None
    for _ in range(it):
        s = xy[rng.choice(len(xy), 3, replace=False)]
        try: cx, cy, r, _ = circ(s)
        except Exception: continue
        if not (0.03 < r < 0.6): continue
        inl = np.abs(np.hypot(xy[:, 0] - cx, xy[:, 1] - cy) - r) < tol
        if best is None or inl.sum() > best.sum(): best = inl
    return circ(xy[best]) if best is not None and best.sum() >= 3 else None


def cyl_fit(S, z13, loss='soft_l1', fs=0.01):
    """S: stem points. Model centre(z) = c0 + a (z - z13), radius(z) = r0 + k (z - z13). Returns radius at z13."""
    r0 = ransac(S[:, :2])
    if r0 is None: return None
    cx, cy, rr, _ = r0
    rr = float(np.clip(rr, 0.021, 0.69))     # start inside the bounds
    def res(p):
        dz = S[:, 2] - z13; X = S[:, 0] - (p[0] + p[2] * dz); Y = S[:, 1] - (p[1] + p[3] * dz)
        return np.hypot(X, Y) - (p[4] + p[5] * dz)
    p = least_squares(res, [cx, cy, 0, 0, rr, -0.005], loss=loss, f_scale=fs,
                      bounds=([-np.inf, -np.inf, -.3, -.3, .02, -.03], [np.inf, np.inf, .3, .3, .7, .005])).x
    inl = np.abs(res(p)) < 0.03
    return dict(r=p[4], taper=p[5], ninl=int(inl.sum()))


TRIM_STUBS_GENERA = ('picea',)   # species that keep dead branches low on the stem


def drone_dbh(P, ground, lo=0.5, hi=4.0, species=''):
    """DBH [cm] of one tree from drone points; None if too few stem points.
    For spruce (dead branch stubs on the lower stem) points more than 3 cm outside the fitted surface are dropped
    and the cylinder refitted twice. On spruce (AEW33_GD, 117 trees kept) this brought the error from 3.4 to 2.5 cm
    and the bias from +2.4 to +0.1 cm; on pine and beech it made things worse (2.3 -> 4.0 cm), and deciding per tree
    from the share of outer points did not work, so it is switched by species."""
    S = P[(P[:, 2] > ground + lo) & (P[:, 2] < ground + hi)]
    if len(S) < 8: return None
    H, xe, ye = np.histogram2d(S[:, 0], S[:, 1], bins=[np.arange(S[:, 0].min(), S[:, 0].max() + .2, .2),
                                                       np.arange(S[:, 1].min(), S[:, 1].max() + .2, .2)])
    i, j = np.unravel_index(H.argmax(), H.shape); c = np.array([xe[i] + .1, ye[j] + .1])
    S = S[np.hypot(S[:, 0] - c[0], S[:, 1] - c[1]) < 0.8]          # points near the densest stem spot
    r = cyl_fit(S, ground + 1.3) if len(S) >= 8 else None
    if r is not None and (species or '').strip().lower().split(' ')[0] in TRIM_STUBS_GENERA:
        for _ in range(2):
            r0 = ransac(S[:, :2])
            if r0 is None: break
            cx, cy, rr, _ = r0; S2 = S[np.hypot(S[:, 0] - cx, S[:, 1] - cy) - rr < 0.03]
            if len(S2) < 8: break
            S = S2; r = cyl_fit(S, ground + 1.3) or r
    return None if r is None else dict(dbh=2 * r['r'] * 100, taper=r['taper'], ninl=r['ninl'], n=len(S))


def plausible(d, min_inliers=20, taper_floor=-0.029, k=2.5):
    """Adds column keep. (1) per species, DBH must fit a robust log(DBH) ~ log(height) + log(crown_d) relation
    fitted on the drone measurements themselves (no ground reference), within k robust SDs;
    (2) at least min_inliers stem points on the cylinder and taper not stuck at its bound."""
    out = []
    for _, g in d.groupby('species'):
        g = g.copy(); ok = np.ones(len(g), bool)
        if len(g) >= 10:
            X = np.c_[np.ones(len(g)), np.log(g.height), np.log(g.crown_d)]; y = np.log(g.dbh.values)
            for _ in range(10):
                b, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None); res = y - X @ b
                s = 1.4826 * np.median(np.abs(res[ok] - np.median(res[ok]))); new = np.abs(res) < k * s
                if (new == ok).all(): break
                ok = new
        g['keep'] = ok & (g.ninl >= min_inliers) & (g.taper > taper_floor); out.append(g)
    return pd.concat(out)


def fill_missing(d, value='dbh', keep='keep'):
    """Gives every tree a DBH. Trees without a plausible measurement get an estimate from height and crown width,
    fitted per species group on the measured trees of the same plot; column dbh_method says which.
    Rejected drone fits were 11-32 cm off; the estimate instead was 2.8 cm (spruce plot) and 3.9 cm (pine/beech
    plot) off, giving 2.6-2.7 cm over all trees."""
    d = d.copy(); ok = d[keep].fillna(False).astype(bool) & d[value].notna()
    d['dbh_final'] = np.where(ok, d[value], np.nan); d['dbh_method'] = np.where(ok, 'measured', 'estimated from height and crown')
    for _, g in d.groupby('species'):
        kg = d.loc[g.index[ok[g.index]]]
        if len(kg) < 5: kg = d[ok]
        if len(kg) < 3: continue
        X = np.c_[np.ones(len(kg)), np.log(kg.height), np.log(kg.crown_d)]
        b, *_ = np.linalg.lstsq(X, np.log(kg[value]), rcond=None)
        m = g.index[~ok[g.index]]
        d.loc[m, 'dbh_final'] = np.exp(np.c_[np.ones(len(m)), np.log(d.loc[m, 'height']), np.log(d.loc[m, 'crown_d'])] @ b)
    return d
