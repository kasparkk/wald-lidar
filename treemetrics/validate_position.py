"""Stem position from the ground scan (lowest stem circle) against the drone scan (tilted cylinder through
the stem points 0.5-4 m) for the same trees, e.g. BioDiv-3DTrees SEW04_GD.
  python validate_position.py TLS_FOLDER ULS_FOLDER
"""
import sys, glob, os, numpy as np, pandas as pd, laspy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import treevol as TV, stem, drone_dbh as DD
from scipy.optimize import least_squares
TLS_DIR, ULS_DIR = sys.argv[1], sys.argv[2]   # folders with <tree>_TLS.laz and <tree>_ULS.laz
def read(f):
    l = laspy.read(f); return np.c_[l.x, l.y, l.z]
def tls_pos(P):
    g = TV.ground_height(P); R = stem.stem_profile(P, g)
    if len(R): return R[0, :2], g, 'circle'
    return np.median(P[P[:, 2] < g + 1.0, :2], axis=0), g, 'lowest metre'
def uls_pos(P, g):
    S = P[(P[:, 2] > g + .5) & (P[:, 2] < g + 4)]
    if len(S) < 8: return None, 'none'
    H, xe, ye = np.histogram2d(S[:, 0], S[:, 1], bins=[np.arange(S[:, 0].min(), S[:, 0].max() + .2, .2), np.arange(S[:, 1].min(), S[:, 1].max() + .2, .2)])
    i, j = np.unravel_index(H.argmax(), H.shape); c = np.array([xe[i] + .1, ye[j] + .1])
    S = S[np.hypot(S[:, 0] - c[0], S[:, 1] - c[1]) < .8]
    if len(S) < 8: return c, 'densest cell'
    r0 = DD.ransac(S[:, :2])
    if r0 is None: return np.median(S[:, :2], 0), 'median'
    cx, cy, rr, _ = r0; rr = float(np.clip(rr, .021, .69)); z13 = g + 1.3
    def res(p):
        dz = S[:, 2] - z13; return np.hypot(S[:, 0] - (p[0] + p[2] * dz), S[:, 1] - (p[1] + p[3] * dz)) - (p[4] + p[5] * dz)
    p = least_squares(res, [cx, cy, 0, 0, rr, -.005], loss='soft_l1', f_scale=.01, bounds=([-np.inf, -np.inf, -.3, -.3, .02, -.03], [np.inf, np.inf, .3, .3, .7, .005])).x
    return np.array([p[0], p[1]]), 'cylinder'
rows = []
for f in sorted(glob.glob(os.path.join(TLS_DIR, '*_TLS.laz'))):
    n = os.path.basename(f).replace('_TLS.laz', ''); P = read(f)
    xy, g, how = tls_pos(P)
    r = dict(tree=n, tx=xy[0], ty=xy[1], tls_how=how, g=g, h=P[:, 2].max() - g, tlsbbox=np.ptp(P[:, :2], 0).max())
    fu = os.path.join(ULS_DIR, n + '_ULS.laz')
    if os.path.exists(fu):
        Q = read(fu); r['nuls'] = len(Q)
        if len(Q) >= 8:
            # ground for ULS: lowest ULS points, else TLS ground
            gu = np.percentile(Q[:, 2], 0.5) if len(Q) > 200 else g
            u, uh = uls_pos(Q, g); r.update(uls_how=uh, gu=gu)
            if u is not None: r.update(ux=u[0], uy=u[1])
            r['top_dx'] = Q[Q[:, 2].argmax(), 0] - P[P[:, 2].argmax(), 0]; r['top_dy'] = Q[Q[:, 2].argmax(), 1] - P[P[:, 2].argmax(), 1]
    rows.append(r)
D = pd.DataFrame(rows); D['dist_m'] = np.hypot(D.ux - D.tx, D.uy - D.ty); D.to_csv('positions.csv', index=False)
print(D.dist_m.describe(percentiles=[.5, .9, .95]).round(3))
