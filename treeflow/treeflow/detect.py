"""Finding single trees.

Two routes, chosen by scan type:
- from above (airborne, drone): canopy height model -> tree tops -> crowns grown down from each top
- from below (ground, handheld): stems cut at breast height -> each stem seeds one tree
Both end in the same place: every vegetation point gets a tree id.
"""
import os
import numpy as np
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max


def canopy_height_model(c, hag, ground, res=0.5):
    x0, y0 = c.x.min(), c.y.min()
    nx = int(np.ceil((c.x.max() - x0) / res)) + 1
    ny = int(np.ceil((c.y.max() - y0) / res)) + 1
    i = np.clip(((c.y - y0) / res).astype(int), 0, ny - 1)
    j = np.clip(((c.x - x0) / res).astype(int), 0, nx - 1)
    chm = np.zeros((ny, nx))
    np.maximum.at(chm, (i, j), np.where(ground, 0, np.clip(hag, 0, None)))
    chm = ndi.grey_closing(chm, size=3)      # fills single empty cells between returns
    return dict(z=chm, x0=float(x0), y0=float(y0), res=res, nx=nx, ny=ny), (i, j)


def tree_tops(chm, min_height=5.0, sigma=0.8, cfg_a=0.8, cfg_b=0.07):
    """Local maxima with a window that grows with tree height (tall trees have wide crowns)."""
    res = chm['res']
    zs = ndi.gaussian_filter(chm['z'], sigma)
    cand = peak_local_max(zs, min_distance=1, threshold_abs=min_height)
    if not len(cand):
        return zs, np.empty((0, 2), int)
    h = zs[cand[:, 0], cand[:, 1]]
    order = np.argsort(-h)
    cand, h = cand[order], h[order]
    radius = np.clip(cfg_a + cfg_b * h, 1.0, 3.5) / res          # metres -> cells
    kept = []
    taken = np.zeros(len(cand), bool)
    for k in range(len(cand)):
        if taken[k]:
            continue
        kept.append(k)
        d = np.hypot(*(cand - cand[k]).T)
        taken |= d <= radius[k]
    return zs, cand[kept]


def stem_seeds(c, hag, ground, z_from=1.0, z_to=1.6, eps=0.1, min_pts=15):
    """Cut a slice at breast height, cluster it and fit a circle to each cluster."""
    from sklearn.cluster import DBSCAN
    m = (~ground) & (hag >= z_from) & (hag <= z_to)
    if m.sum() < min_pts:
        return []
    xy = np.c_[c.x[m], c.y[m]]
    lab = DBSCAN(eps=eps, min_samples=min_pts // 3).fit_predict(xy)
    seeds = []
    for k in set(lab) - {-1}:
        p = xy[lab == k]
        if len(p) < min_pts:
            continue
        fit = fit_circle(p)
        if fit is None:
            continue
        cx, cy, r, rmse = fit
        if 0.03 <= r <= 0.8 and rmse < max(0.02, 0.25 * r):
            seeds.append(dict(x=cx, y=cy, r=r, rmse=rmse, n=len(p)))
    return seeds


try:  # optional: same circle fit in Rust (treeflow/rust_kreis), several times faster
    if os.environ.get('TREEFLOW_OHNE_RUST'):
        raise ImportError
    from treeflow_rust import fit_circle as _fit_circle_rust
except ImportError:
    _fit_circle_rust = None


def fit_circle(p, iters=200, tol=0.02, rng=np.random.default_rng(0)):
    """RANSAC circle fit, then least squares on the inliers. Returns (cx, cy, r, rmse) or None.
    Uses the Rust version if it is installed (pip install treeflow/rust_kreis), else Python."""
    if _fit_circle_rust is not None:
        return _fit_circle_rust(np.ascontiguousarray(p[:, :2], dtype=np.float64), iters, tol)
    return fit_circle_python(p, iters, tol, rng)


def fit_circle_python(p, iters=200, tol=0.02, rng=np.random.default_rng(0)):
    """The original Python version (reference for the Rust version)."""
    if len(p) < 6:
        return None
    best = None
    for _ in range(iters):
        a, b, cc = p[rng.choice(len(p), 3, replace=False)]
        d = 2 * (a[0] * (b[1] - cc[1]) + b[0] * (cc[1] - a[1]) + cc[0] * (a[1] - b[1]))
        if abs(d) < 1e-9:
            continue
        ux = ((a @ a) * (b[1] - cc[1]) + (b @ b) * (cc[1] - a[1]) + (cc @ cc) * (a[1] - b[1])) / d
        uy = ((a @ a) * (cc[0] - b[0]) + (b @ b) * (a[0] - cc[0]) + (cc @ cc) * (b[0] - a[0])) / d
        r = np.hypot(a[0] - ux, a[1] - uy)
        if r > 1.0:
            continue
        inl = np.abs(np.hypot(p[:, 0] - ux, p[:, 1] - uy) - r) < tol
        if best is None or inl.sum() > best.sum():
            best = inl
    if best is None or best.sum() < 6:
        return None
    q = p[best]
    A = np.c_[2 * q[:, 0], 2 * q[:, 1], np.ones(len(q))]
    bb = (q ** 2).sum(1)
    sol, *_ = np.linalg.lstsq(A, bb, rcond=None)
    cx, cy = sol[0], sol[1]
    r = np.sqrt(max(sol[2] + cx ** 2 + cy ** 2, 0))
    rmse = float(np.sqrt(np.mean((np.hypot(q[:, 0] - cx, q[:, 1] - cy) - r) ** 2)))
    return float(cx), float(cy), float(r), rmse


def segment(c, hag, ground, kind, cfg):
    """Returns per-point tree id (0 = no tree), the list of seeds and the CHM used."""
    res = cfg.get('chm_res', 0.5 if kind == 'airborne' else 0.25)
    chm, (pi, pj) = canopy_height_model(c, hag, ground, res)
    det = cfg.get('detector', 'treeflow')
    if det.startswith('lidr'):
        from . import lidr
        tree_id, seeds = lidr.segment(c, hag, ground, res, cfg)
        # crown map: each cell belongs to the tree of its highest point
        lab = np.zeros(chm['z'].shape, int)
        m = tree_id > 0
        order = np.argsort(hag[m])                       # later (higher) points overwrite lower ones
        lab[pi[m][order], pj[m][order]] = tree_id[m][order]
        names = {'lidr-dalponte': 'lidR: local maxima + Dalponte 2016 crowns',
                 'lidr-silva': 'lidR: local maxima + Silva 2016 crowns',
                 'lidr-li2012': 'lidR: Li et al. 2012 point-cloud segmentation'}
        return tree_id, seeds, chm, lab, names.get(det, det)
    zs, tops = tree_tops(chm, cfg.get('min_tree_height', 5.0), cfg.get('chm_smooth', 0.8), cfg.get('window_a', 0.8), cfg.get('window_b', 0.07))
    method = 'tree tops in canopy height model'
    stems = []
    if kind != 'airborne' or cfg.get('use_stems', False):
        stems = stem_seeds(c, hag, ground)
        if len(stems) >= max(3, 0.3 * len(tops)):
            method = 'stems at breast height'
    markers = np.zeros(chm['z'].shape, int)
    seeds = []
    if method.startswith('stems'):
        for k, s in enumerate(stems, 1):
            i = int((s['y'] - chm['y0']) / res); j = int((s['x'] - chm['x0']) / res)
            markers[np.clip(i, 0, chm['ny'] - 1), np.clip(j, 0, chm['nx'] - 1)] = k
            seeds.append(dict(id=k, x=s['x'], y=s['y'], stem=s))
    else:
        for k, (i, j) in enumerate(tops, 1):
            markers[i, j] = k
            seeds.append(dict(id=k, x=chm['x0'] + (j + 0.5) * res, y=chm['y0'] + (i + 0.5) * res, stem=None))
        # attach stems found by the slice to the nearest top (gives DBH where stems are visible)
        for s in stems:
            if not seeds:
                break
            d = [np.hypot(s['x'] - t['x'], s['y'] - t['y']) for t in seeds]
            k = int(np.argmin(d))
            if d[k] < 2.0 and (seeds[k]['stem'] is None or s['n'] > seeds[k]['stem']['n']):
                seeds[k]['stem'] = s
    mask = zs > cfg.get('crown_min_height', 2.0)
    lab = watershed(-zs, markers, mask=mask)
    tree_id = np.where(ground | (hag < cfg.get('min_point_height', 1.0)), 0, lab[pi, pj])
    return tree_id, seeds, chm, lab, method
