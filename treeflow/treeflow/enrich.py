"""Optional: add aerial-photo values per crown (needs a 4-band RGB + near-infrared image in the
same coordinate system as the point cloud, e.g. LGB DOP20 RGBI)."""
import numpy as np


def ortho_features(rows, lab, chm, image_path):
    import rasterio
    from rasterio.windows import from_bounds
    res = chm['res']
    x0, y0 = chm['x0'], chm['y0']
    x1, y1 = x0 + chm['nx'] * res, y0 + chm['ny'] * res
    with rasterio.open(image_path) as d:
        if d.count < 4:
            return rows, 'image has no near-infrared band; skipped'
        w = from_bounds(x0, y0, x1, y1, d.transform)
        A = d.read([1, 2, 3, 4], window=w, boundless=True, fill_value=0).astype(float)
        tr = d.window_transform(w)
    R, G, B, N = A
    rr, cc = np.mgrid[0:N.shape[0], 0:N.shape[1]]
    px = tr.c + (cc + 0.5) * tr.a
    py = tr.f + (rr + 0.5) * tr.e
    li = np.clip(((py - y0) / res).astype(int), 0, chm['ny'] - 1)
    lj = np.clip(((px - x0) / res).astype(int), 0, chm['nx'] - 1)
    L = lab[li, lj]
    ndvi = (N - R) / (N + R + 1e-6)
    flat_l, flat_n, flat_v = L.ravel(), N.ravel(), ndvi.ravel()
    order = np.argsort(flat_l, kind='stable')
    starts = np.searchsorted(flat_l[order], np.arange(lab.max() + 2))
    for r in rows:
        k = r['tree']
        idx = order[starts[k]:starts[k + 1]]
        if len(idx) < 10:
            r['ndvi'] = None; r['nir'] = None
            continue
        nv = flat_n[idx]; sunlit = nv >= np.median(nv)      # ignore the shaded half of the crown
        r['nir'] = round(float(nv[sunlit].mean()), 1)
        r['ndvi'] = round(float(flat_v[idx][sunlit].mean()), 3)
    vals = np.array([r['ndvi'] for r in rows if r.get('ndvi') is not None])
    if len(vals) >= 5:
        med = np.median(vals); mad = 1.4826 * np.median(np.abs(vals - med))
        for r in rows:
            v = r.get('ndvi')
            if v is None:
                r['type_photo'] = 'unknown'
            elif v > med + 2.5 * max(mad, 0.02):
                r['type_photo'] = 'differs (broadleaf?)'
            else:
                r['type_photo'] = 'like the stand'
    return rows, 'ok'
