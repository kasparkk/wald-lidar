"""Writing results: table, GeoJSON, summary, a one-page report and data for the Click-a-Tree viewer."""
import csv, json, html
import numpy as np


def write_csv(rows, path):
    if not rows:
        open(path, 'w').write('')
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)


def write_geojson(rows, crs, path):
    to_ll = None
    if crs:
        try:
            from pyproj import Transformer
            to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
        except Exception:
            to_ll = None
    feats = []
    for r in rows:
        if to_ll:
            lon, lat = to_ll.transform(r['x'], r['y'])
            geom = dict(type='Point', coordinates=[round(lon, 7), round(lat, 7)])
        else:
            geom = dict(type='Point', coordinates=[r['x'], r['y']])
        feats.append(dict(type='Feature', geometry=geom, properties=r))
    json.dump(dict(type='FeatureCollection', features=feats,
                   crs_note='WGS84 lon/lat' if to_ll else f'file coordinates ({crs or "unknown CRS"})'),
              open(path, 'w'))


def summary(rows, info, area_ha):
    h = np.array([r['height_m'] for r in rows]) if rows else np.array([])
    d = np.array([r['dbh_cm'] for r in rows if r['dbh_cm'] and r.get('dbh_source') == 'stem circle fit'])
    return dict(
        trees=len(rows), area_ha=round(area_ha, 3),
        trees_per_ha=round(len(rows) / area_ha, 1) if area_ha else None,
        height_median_m=round(float(np.median(h)), 1) if len(h) else None,
        height_p95_m=round(float(np.percentile(h, 95)), 1) if len(h) else None,
        stems_measured=int(len(d)),
        dbh_median_cm=round(float(np.median(d)), 1) if len(d) else None,
        **{k: (round(info[k], 1) if k == 'density' else info[k]) for k in ('kind', 'density', 'method', 'ground_source') if k in info})


def _hist_svg(values, lo=0, hi=40, step=2, w=520, hgt=160):
    bins = np.arange(lo, hi + step, step)
    cnt, _ = np.histogram(values, bins)
    m = max(cnt.max(), 1) if len(cnt) else 1
    bw = (w - 40) / len(cnt)
    bars = ''.join(f'<rect x="{30 + i * bw:.1f}" y="{hgt - 20 - c / m * (hgt - 40):.1f}" width="{bw - 2:.1f}" '
                   f'height="{c / m * (hgt - 40):.1f}" fill="var(--accent)"><title>{bins[i]}–{bins[i + 1]} m: {c}</title></rect>'
                   for i, c in enumerate(cnt))
    ticks = ''.join(f'<text x="{30 + i * bw:.1f}" y="{hgt - 4}" font-size="11" fill="var(--muted)">{bins[i]}</text>'
                    for i in range(0, len(cnt), 5))
    return (f'<svg viewBox="0 0 {w} {hgt}" width="100%" role="img" aria-label="Height histogram">'
            f'{bars}{ticks}<text x="{w - 40}" y="{hgt - 4}" font-size="11" fill="var(--muted)">m</text></svg>')


def write_report(rows, summ, info, log, validation, path):
    e = html.escape
    kv = ''.join(f'<dt>{e(str(k))}</dt><dd>{e(str(v))}</dd>' for k, v in summ.items())
    steps = ''.join(f'<li>{e(s)}</li>' for s in log)
    val = ''
    if validation:
        val = '<h2>Check against reference trees</h2><dl>' + ''.join(
            f'<dt>{e(k)}</dt><dd>{e(str(v))}</dd>' for k, v in validation.items()) + '</dl>'
    warn = ''.join(f'<li>{e(w)}</li>' for w in info.get('warnings', []))
    doc = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>treeflow report</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>:root{{--bg:#f4f6f2;--fg:#17201b;--muted:#4f5c55;--accent:#3d7a55;--line:#cfd8d1}}
@media (prefers-color-scheme:dark){{:root{{--bg:#111714;--fg:#e4ebe6;--muted:#a3b2aa;--accent:#7fbf8f;--line:#2b3631}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif;max-width:760px;margin:0 auto;padding:24px 16px}}
h1{{margin:0 0 4px}} h2{{margin-top:28px;font-size:18px}} dl{{display:grid;grid-template-columns:auto 1fr;gap:4px 16px}}
dt{{color:var(--muted)}} dd{{margin:0;font-family:ui-monospace,monospace}} .src{{color:var(--muted);font-size:13px}}</style>
<h1>treeflow report</h1><p class="src">{e(info.get('source',''))}</p>
<h2>Result</h2><dl>{kv}</dl>
<h2>Tree heights</h2>{_hist_svg([r['height_m'] for r in rows])}
{val}
<h2>What was done</h2><ol>{steps}</ol>
<h2>Limits</h2><ul>{warn}</ul></html>"""
    open(path, 'w').write(doc)


def viewer_json(c, hag, tree_id, ground, dtm, rows, bbox, path, max_pts_tree=900, seed=1):
    """Data for the Click-a-Tree page: local metres, height above the lowest ground in the patch."""
    from .ground import dtm_at
    x0, y0, x1, y1 = bbox
    S = float(max(x1 - x0, y1 - y0))
    rng = np.random.default_rng(seed)
    core = (c.x >= x0) & (c.x < x1) & (c.y >= y0) & (c.y < y1)
    gmin = float(np.min(dtm_at(dtm, c.x[core & ground], c.y[core & ground]))) if (core & ground).any() else float(c.z[core].min())
    out = dict(S=S, trees=[], ground=[], other=[])
    order = np.argsort(tree_id, kind='stable'); ts = tree_id[order]
    for r in rows:
        k = r['tree']
        lo, hi = np.searchsorted(ts, k), np.searchsorted(ts, k + 1)
        sel = order[lo:hi]; sel = sel[core[sel]]
        if len(sel) > max_pts_tree:
            sel = rng.choice(sel, max_pts_tree, replace=False)
        g = float(dtm_at(dtm, np.array([r['x']]), np.array([r['y']]))[0]) - gmin
        out['trees'].append(dict(x=round(r['x'] - x0, 2), y=round(r['y'] - y0, 2), g=round(g, 2), h=r['height_m'],
                                 cw=r['crown_diameter_m'], cb=r['crown_base_m'] or 0, n=r['points'],
                                 below=int(r['points'] * r['later_returns_pct'] / 100), dbh=r['dbh_cm'],
                                 dbh_src=r['dbh_source'], ndvi=r.get('ndvi'), type=r.get('type_photo'),
                                 p=np.round(np.c_[c.x[sel] - x0, c.z[sel] - gmin, c.y[sel] - y0], 2).ravel().tolist()))
    gs = np.nonzero(core & ground)[0]
    if len(gs) > 9000: gs = rng.choice(gs, 9000, replace=False)
    out['ground'] = np.round(np.c_[c.x[gs] - x0, c.z[gs] - gmin, c.y[gs] - y0], 2).ravel().tolist()
    os_ = np.nonzero(core & ~ground & (tree_id == 0) & (hag > 1.0))[0]
    if len(os_) > 4000: os_ = rng.choice(os_, 4000, replace=False)
    out['other'] = np.round(np.c_[c.x[os_] - x0, c.z[os_] - gmin, c.y[os_] - y0], 2).ravel().tolist()
    n = int(S / 2) + 1
    tz = [round(float(dtm_at(dtm, np.array([x0 + xx * 2]), np.array([y0 + yy * 2]))[0]) - gmin, 2)
          for yy in range(n) for xx in range(n)]
    out['terrain'] = dict(step=2, n=n, z=tz)
    out['meta'] = dict(x0=x0, y0=y0)
    json.dump(out, open(path, 'w'), separators=(',', ':'))
