"""Turn trees that are already segmented (by 3Dtrees.earth, SegmentAnyTree, lidR, FSCT, treeflow …)
into compact 3D tree models for a web viewer.

Input: one point cloud with a per-point tree id, or a folder of single-tree files.
Output: a bundle folder
    index.json                      sites, tiles, attribute list
    <site>/tile_<i>_<j>.json        tree models + attributes for one tile
    <site>/pts_<i>_<j>.json         a few hundred points per tree, loaded only when a tree is opened
    <site>/terrain.json             coarse ground model

Each tree model is small (≈ 100 numbers): crown as rings of radii at several heights, stem as a list
of measured (or estimated) radii, plus its values. That keeps whole regions loadable in a browser.
"""
import json, os, glob
import numpy as np
import laspy
from . import ground as groundmod
from .detect import fit_circle

TREE_FIELDS = ['treeID', 'treeid', 'tree_id', 'TreeID', 'PredInstance', 'instance', 'instance_id',
               'pred_instance', 'segment_id', 'tid']
L, K = 6, 12                     # crown rings and directions per ring


def _find_field(las, name=None):
    dims = list(las.point_format.dimension_names)
    if name:
        if name not in dims:
            raise ValueError(f'field {name!r} not in file; fields: {dims}')
        return name
    for f in TREE_FIELDS:
        if f in dims:
            return f
    raise ValueError(f'no tree id field found; fields: {dims}. Pass --tree-field.')


def load_segmented(path, tree_field=None):
    """Returns x, y, z, cls, tree_id arrays. A folder of single-tree files gets ids by file order."""
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, '*.la[sz]')))
        xs, ys, zs, cs, ts, names = [], [], [], [], [], {}
        for k, f in enumerate(files, 1):
            las = laspy.read(f)
            n = len(las.x)
            xs.append(np.asarray(las.x)); ys.append(np.asarray(las.y)); zs.append(np.asarray(las.z))
            cs.append(np.asarray(las.classification)); ts.append(np.full(n, k))
            names[k] = os.path.splitext(os.path.basename(f))[0]
        return (np.concatenate(xs), np.concatenate(ys), np.concatenate(zs), np.concatenate(cs),
                np.concatenate(ts), names, 'one file per tree')
    las = laspy.read(path)
    f = _find_field(las, tree_field)
    t = np.asarray(las[f]).astype(np.int64)
    t[t < 0] = 0
    return (np.asarray(las.x), np.asarray(las.y), np.asarray(las.z), np.asarray(las.classification),
            t, {}, f'field "{f}"')


def load_rgb(path):
    """Camera colour per point (0-255), or None if the file has none."""
    if os.path.isdir(path):
        return None
    las = laspy.read(path)
    if not {'red', 'green', 'blue'} <= set(las.point_format.dimension_names):
        return None
    rgb = np.c_[np.asarray(las.red), np.asarray(las.green), np.asarray(las.blue)].astype(np.float32)
    if rgb.max() > 255:
        rgb /= 256.0                                   # 16-bit colour -> 8-bit
    if rgb.max() == 0:
        return None
    return rgb.astype(np.uint8)


def crown_rings(px, py, h, cx, cy, base, top):
    """Radius of the crown in K directions at L heights (90th percentile of the points)."""
    dh = (top - base) / L
    if dh <= 0.1 or len(h) < 8:
        return None
    li = np.clip(((h - base) / dh).astype(int), 0, L - 1)
    ki = (((np.arctan2(py - cy, px - cx) + np.pi) / (2 * np.pi)) * K).astype(int) % K
    d = np.hypot(px - cx, py - cy)
    inside = (h >= base) & (h <= top)
    R = np.full((L, K), np.nan)
    for l in range(L):
        for k in range(K):
            m = inside & (li == l) & (ki == k)
            if m.sum():
                R[l, k] = np.percentile(d[m], 90)
    # fill gaps around each ring, then from neighbouring rings
    for l in range(L):
        row = R[l]
        if np.isnan(row).all():
            continue
        idx = np.nonzero(~np.isnan(row))[0]
        for k in np.nonzero(np.isnan(row))[0]:
            dist = np.minimum(np.abs(idx - k), K - np.abs(idx - k))
            row[k] = row[idx[np.argsort(dist)[:2]]].mean()
    for _ in range(2):
        for l in range(L):
            if np.isnan(R[l]).all():
                nb = [R[j] for j in (l - 1, l + 1) if 0 <= j < L and not np.isnan(R[j]).all()]
                if nb:
                    R[l] = np.mean(nb, axis=0)
    if np.isnan(R).any():
        R = np.where(np.isnan(R), np.nanmedian(R) if not np.isnan(R).all() else 1.0, R)
    # light smoothing
    S = R.copy()
    for l in range(L):
        for k in range(K):
            nb = [R[l, (k + 1) % K], R[l, (k - 1) % K], R[max(l - 1, 0), k], R[min(l + 1, L - 1), k]]
            S[l, k] = 0.5 * R[l, k] + 0.125 * sum(nb)
    return np.clip(S, 0.2, 25.0)


def stem_profile(px, py, h, cx, cy, crown_base, dense, tree_height=None):
    """Measured stem: circle fits every metre below the crown (dense scans only).
    A fit only counts if it looks like a real stem cross-section: most points lie on the circle,
    they surround it (not just one side), and the size is plausible for the tree's height."""
    out = []
    if not dense:
        return out
    r_max = 0.6
    if tree_height:
        r_max = min(r_max, (2.0 * tree_height + 10) / 200)       # e.g. 10 m tree -> max 30 cm diameter
    r_prev = None
    for z in np.arange(0.5, max(1.0, min(crown_base, 15.0)), 1.0):
        m = (h >= z) & (h < z + 0.3) & (np.hypot(px - cx, py - cy) < 0.8)
        if m.sum() < 12:
            continue
        P = np.c_[px[m], py[m]]
        fit = fit_circle(P)
        if not fit:
            continue
        fx, fy, r, rmse = fit
        d = np.abs(np.hypot(P[:, 0] - fx, P[:, 1] - fy) - r)
        inl = d < max(0.02, 0.15 * r)
        if inl.mean() < 0.6 or not (0.02 < r < r_max) or rmse > 0.25 * r:
            continue
        ang = np.sort(np.arctan2(P[inl, 1] - fy, P[inl, 0] - fx))
        gap = np.max(np.diff(np.r_[ang, ang[0] + 2 * np.pi])) if len(ang) > 2 else 2 * np.pi
        if gap > np.deg2rad(200):            # points only on one side: not trustworthy
            continue
        if r_prev is not None and r > 1.3 * r_prev:
            continue
        out.append([round(float(z + 0.15), 2), round(r * 100, 1), round(fx, 2), round(fy, 2)])
        r_prev = r
    return out


def build_models(x, y, z, cls, tid, cfg=None, rgb=None):
    cfg = cfg or {}
    # height above ground: ground class if present, else treeflow's ground filter on the whole cloud
    from .io import Cloud
    c = Cloud(x, y, z, cls.astype(np.uint8), np.ones(len(x), np.uint8), np.ones(len(x), np.uint8),
              np.zeros(len(x), np.float32))
    hag, g, dtm, how = groundmod.normalise(c, {'force_ground_filter': cfg.get('force_ground_filter', False)})
    area = max((x.max() - x.min()) * (y.max() - y.min()), 1)
    dense = bool(len(x) / area > 150)
    order = np.argsort(tid, kind='stable'); ts = tid[order]
    ids = np.unique(tid[tid > 0])
    trees, rng = [], np.random.default_rng(0)
    for k in ids:
        lo, hi = np.searchsorted(ts, k), np.searchsorted(ts, k + 1)
        idx = order[lo:hi]
        idx = idx[~g[idx] & (hag[idx] > 0.5)]
        if len(idx) < cfg.get('min_points_per_tree', 15):
            continue
        h = hag[idx]; px, py = x[idx], y[idx]
        top = float(np.percentile(h, 99.5))
        hi_pts = h >= np.percentile(h, 90)
        cx, cy = float(px[hi_pts].mean()), float(py[hi_pts].mean())
        # stem axis from the low points if there are enough (dense scans)
        low = h < 3.0
        if dense and low.sum() >= 20:
            cx, cy = float(np.median(px[low])), float(np.median(py[low]))
        away = np.hypot(px - cx, py - cy) > 0.8
        cb = float(np.percentile(h[away & (h > 1.5)], 5)) if (away & (h > 1.5)).sum() >= 10 else top * 0.5
        rings = crown_rings(px, py, h, cx, cy, cb, top)
        stem = stem_profile(px, py, h, cx, cy, cb, dense, top)
        dbh = None
        if stem:
            near = min(stem, key=lambda s: abs(s[0] - 1.3))
            if abs(near[0] - 1.3) < 1.0:
                dbh = round(2 * near[1], 1)
        cw = float(2 * np.mean(rings.max(0))) if rings is not None else float(2 * np.percentile(np.hypot(px - cx, py - cy), 90))
        g0 = float(groundmod.dtm_at(dtm, np.array([cx]), np.array([cy]))[0])
        sel = idx if len(idx) <= 250 else rng.choice(idx, 250, replace=False)
        colour, pcol = None, None
        if rgb is not None:
            upper = idx[h >= cb + 0.6 * (top - cb)]           # sunlit upper crown, as seen from above
            if len(upper) >= 5:
                colour = '#%02x%02x%02x' % tuple(int(v) for v in np.median(rgb[upper], axis=0))
            c = rgb[sel].astype(np.int64)
            pcol = ((c[:, 0] << 16) | (c[:, 1] << 8) | c[:, 2]).tolist()
        trees.append(dict(
            id=int(k), x=round(cx, 2), y=round(cy, 2), g=round(g0, 2), h=round(top, 2), cb=round(cb, 2),
            cw=round(cw, 2), n=int(len(idx)), dbh=dbh, dbh_src='measured (stem fit)' if dbh else None,
            rings=None if rings is None else np.round(rings * 10).astype(int).ravel().tolist(),   # decimetres
            stem=[[s[0], s[1]] for s in stem],
            pts=np.round(np.c_[x[sel] - cx, z[sel] - (g0 + 0), y[sel] - cy] * 10).astype(int).ravel().tolist(),
            **({'rgb': colour} if colour else {}), **({'pcol': pcol} if pcol else {})))
    info = dict(ground_source=how, dense=dense, density=round(len(x) / area, 1))
    if rgb is not None and g.any():
        # mean camera colour of the ground points in 5 m cells, for colouring the terrain
        st = 5.0; gx0, gy0 = float(x.min()), float(y.min())
        nx = int((x.max() - gx0) // st) + 1; ny = int((y.max() - gy0) // st) + 1
        ci = ((y[g] - gy0) // st).astype(int) * nx + ((x[g] - gx0) // st).astype(int)
        cnt = np.bincount(ci, minlength=nx * ny)
        mean = np.stack([np.bincount(ci, rgb[g][:, k].astype(float), minlength=nx * ny) for k in range(3)], 1) / np.maximum(cnt, 1)[:, None]
        info['_ground_rgb'] = dict(x0=gx0, y0=gy0, step=st, nx=nx, ny=ny, rgb=mean, has=cnt > 0)
    return trees, dtm, info


def write_bundle(trees, dtm, info, out_dir, site, title, tile=100.0, labels=None, crs=None):
    os.makedirs(os.path.join(out_dir, site), exist_ok=True)
    if labels:
        for t in trees:
            t.update(labels.get(t['id'], {}))
    x0 = min(t['x'] for t in trees); y0 = min(t['y'] for t in trees)
    x0 = np.floor(x0 / tile) * tile; y0 = np.floor(y0 / tile) * tile
    tiles = {}
    for t in trees:
        key = (int((t['x'] - x0) // tile), int((t['y'] - y0) // tile))
        tiles.setdefault(key, []).append(t)
    tile_list = []
    for (i, j), ts in tiles.items():
        models = [{k: v for k, v in t.items() if k not in ('pts', 'pcol')} for t in ts]
        json.dump(dict(trees=models), open(os.path.join(out_dir, site, f'tile_{i}_{j}.json'), 'w'), separators=(',', ':'))
        json.dump({str(t['id']): ({'p': t['pts'], 'c': t['pcol']} if t.get('pcol') else t['pts']) for t in ts}, open(os.path.join(out_dir, site, f'pts_{i}_{j}.json'), 'w'),
                  separators=(',', ':'))
        tile_list.append(dict(i=i, j=j, n=len(ts), hmax=round(max(t['h'] for t in ts), 1)))
    # coarse terrain for the ground surface (5 m)
    step = 5.0
    xs = np.arange(x0, x0 + (max(k[0] for k in tiles) + 1) * tile + 0.1, step)
    ys = np.arange(y0, y0 + (max(k[1] for k in tiles) + 1) * tile + 0.1, step)
    X, Y = np.meshgrid(xs, ys)
    Z = groundmod.dtm_at(dtm, X.ravel(), Y.ravel())
    terr = dict(x0=float(x0), y0=float(y0), step=step, nx=len(xs), ny=len(ys), z=np.round(Z, 2).tolist())
    G = info.pop('_ground_rgb', None)
    if G is not None:
        gi = np.clip(((Y.ravel() - G['y0']) // G['step']).astype(int), 0, G['ny'] - 1)
        gj = np.clip(((X.ravel() - G['x0']) // G['step']).astype(int), 0, G['nx'] - 1)
        k = gi * G['nx'] + gj
        col = np.where(G['has'][k][:, None], G['rgb'][k], np.array([120, 120, 100])).astype(int)
        terr['c'] = ((col[:, 0] << 16) | (col[:, 1] << 8) | col[:, 2]).tolist()
    json.dump(terr, open(os.path.join(out_dir, site, 'terrain.json'), 'w'), separators=(',', ':'))
    extra = sorted({k for t in trees for k in t} - {'id', 'x', 'y', 'g', 'h', 'cb', 'cw', 'n', 'dbh', 'dbh_src', 'rings', 'stem', 'pts', 'pcol', 'rgb'})
    site_info = dict(id=site, title=title, x0=float(x0), y0=float(y0), tile=tile, tiles=tile_list, trees=len(trees),
                     rings=[L, K], crs=crs, extra_fields=extra, **info)
    ip = os.path.join(out_dir, 'index.json')
    index = json.load(open(ip)) if os.path.exists(ip) else dict(sites=[])
    index['sites'] = [s for s in index['sites'] if s['id'] != site] + [site_info]
    json.dump(index, open(ip, "w"), indent=1, default=float)
    return site_info


def load_labels(path):
    """CSV with a column 'id' (tree id) plus any values to show, e.g. species, source."""
    import csv
    out = {}
    for r in csv.DictReader(open(path)):
        k = int(float(r.pop('id')))
        out[k] = {a: (float(v) if v.replace('.', '', 1).replace('-', '', 1).isdigit() else v) for a, v in r.items() if v != ''}
    return out


def visualise(path, out_dir, site=None, title=None, tree_field=None, labels=None, tile=100.0, cfg=None):
    x, y, z, cls, tid, names, how = load_segmented(path, tree_field)
    rgb = load_rgb(path)
    trees, dtm, info = build_models(x, y, z, cls, tid, cfg, rgb)
    info['has_colour'] = rgb is not None
    info['trees_from'] = how
    lab = load_labels(labels) if labels else {}
    if names:
        for k, nm in names.items():
            lab.setdefault(k, {})['name'] = nm
    site = site or os.path.splitext(os.path.basename(path.rstrip('/')))[0]
    try:
        crs = laspy.read(path).header.parse_crs().to_string() if not os.path.isdir(path) else None
    except Exception:
        crs = None
    return write_bundle(trees, dtm, info, out_dir, site, title or site, tile, lab, crs)
