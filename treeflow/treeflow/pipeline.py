"""The automatic chain: file in, trees out. Each step logs what it did so the report explains itself."""
import json, os, time
import numpy as np
from . import io, ground, detect, metrics, enrich, export, validate

DEFAULTS = dict(
    buffer_m=10.0,            # extra margin read around a bbox so edge trees keep their whole crown
    min_tree_height=5.0,      # tops lower than this are ignored
    crown_min_height=2.0,
    min_point_height=1.0,     # points below this (above ground) are not assigned to trees
    min_points_per_tree=15,
    estimate_dbh=False,       # True: fill missing DBH from height (labelled as estimate)
    viewer=True,              # write viewer.json for Click-a-Tree (only for areas <= 150 m)
    detector='treeflow',      # or lidr-dalponte, lidr-silva, lidr-li2012 (needs R + lidR)
)


def run(path, out_dir, cfg=None, bbox=None, ortho=None, reference=None):
    cfg = {**DEFAULTS, **(cfg or {})}
    os.makedirs(out_dir, exist_ok=True)
    log, warnings, t0 = [], [], time.time()
    step = lambda msg: log.append(f'{msg}  [{time.time() - t0:.1f} s]')

    c = io.read(path, bbox=bbox, buffer=cfg['buffer_m'] if bbox else 0)
    info = io.describe(c)
    info['source'] = str(path)
    step(f"Read {info['points']:,} points, {info['density']:.0f} per m² → treated as {info['kind']} scan")
    if bbox is None:
        bbox = info['extent']
        warnings.append('Whole file processed: trees cut by the file edge are measured only partly.')

    hag, g, dtm, how = ground.normalise(c, cfg)
    info['ground_source'] = how
    step(f'Ground from {how}; terrain model at {dtm["res"]} m; heights converted to height above ground')

    if cfg.get('detector', 'treeflow').startswith('lidr'):
        from . import lidr
        if not lidr.available():
            warnings.append('lidR not installed: fell back to treeflow\'s own detector.')
            cfg['detector'] = 'treeflow'
    tree_id, seeds, chm, lab, method = detect.segment(c, hag, g, info['kind'], cfg)
    info['method'] = method
    step(f'Trees found by {method}: {len(seeds)} candidates; crowns split by watershed on a {chm["res"]} m canopy model')

    rows = metrics.tree_table(c, hag, tree_id, seeds, chm, lab, cfg, core_bbox=bbox)
    n_dbh = sum(1 for r in rows if r['dbh_cm'])
    step(f'Measured {len(rows)} trees inside the area; stem diameter measured for {n_dbh}')
    if cfg['estimate_dbh']:
        metrics.add_dbh_estimate(rows)
        step('Filled missing DBH with a height-based estimate (marked in dbh_source)')

    if ortho:
        rows, msg = enrich.ortho_features(rows, lab, chm, ortho)
        step(f'Aerial photo values per crown: {msg}')

    validation = None
    if reference:
        validation = validate.compare(rows, validate.load_reference(reference))
        step('Compared with reference trees')

    if info['kind'] == 'airborne':
        warnings += ['From above, trees under the canopy are often missed and neighbouring crowns can merge.',
                     'Stem diameter is rarely measurable from airborne data.']
    warnings.append('Species is not derived from geometry; aerial photo values only flag crowns that differ from the stand.')
    if not reference:
        warnings.append('Not checked against reference trees. Run with --reference to get detection and error rates.')
    info['warnings'] = warnings

    area_ha = (bbox[2] - bbox[0]) * (bbox[3] - bbox[1]) / 1e4
    summ = export.summary(rows, info, area_ha)
    export.write_csv(rows, os.path.join(out_dir, 'trees.csv'))
    export.write_geojson(rows, c.crs, os.path.join(out_dir, 'trees.geojson'))
    json.dump(dict(summary=summ, info=info, validation=validation, log=log, config=cfg),
              open(os.path.join(out_dir, 'summary.json'), 'w'), indent=1, default=str)
    if cfg.get('write_segmented'):
        keep = np.isin(tree_id, [r['tree'] for r in rows])
        tid_out = np.where(keep, tree_id, 0)
        core = (c.x >= bbox[0]) & (c.x < bbox[2]) & (c.y >= bbox[1]) & (c.y < bbox[3])
        import laspy
        h = laspy.LasHeader(point_format=1, version='1.4'); h.scales = [0.01] * 3
        h.offsets = [c.x.min(), c.y.min(), c.z.min()]
        h.add_extra_dim(laspy.ExtraBytesParams(name='treeID', type=np.int32))
        las = laspy.LasData(h)
        las.x, las.y, las.z = c.x[core], c.y[core], c.z[core]
        las.classification = np.where(g[core], 2, 1).astype(np.uint8)
        las.treeID = tid_out[core].astype(np.int32)
        las.write(os.path.join(out_dir, 'segmented.laz'))
        step('Wrote segmented.laz (every point with its tree id, field "treeID")')
    if cfg['viewer'] and max(bbox[2] - bbox[0], bbox[3] - bbox[1]) <= 150:
        export.viewer_json(c, hag, tree_id, g, dtm, rows, bbox, os.path.join(out_dir, 'viewer.json'))
        step('Wrote viewer.json for Click-a-Tree')
    step('Wrote trees.csv, trees.geojson, summary.json, report.html')
    export.write_report(rows, summ, info, log, validation, os.path.join(out_dir, 'report.html'))
    return rows, summ, validation


def run_folder(folder, out_dir, **kw):
    """Process every .las/.laz in a folder; one sub-folder of results per file plus an overview CSV."""
    files = sorted(f for f in os.listdir(folder) if f.lower().endswith(('.las', '.laz')))
    overview = []
    for f in files:
        _, s, v = run(os.path.join(folder, f), os.path.join(out_dir, os.path.splitext(f)[0]), **kw)
        overview.append(dict(file=f, **s, **({k: v[k] for k in v} if v else {})))
    export.write_csv(overview, os.path.join(out_dir, 'overview.csv'))
    return overview


def run_compare(path, out_dir, detectors, cfg=None, **kw):
    """Same input through several detectors; one table to see where they agree."""
    table = []
    for d in detectors:
        _, s, v = run(path, os.path.join(out_dir, d), cfg={**(cfg or {}), 'detector': d}, **kw)
        table.append(dict(detector=d, **{k: s[k] for k in ('trees', 'trees_per_ha', 'height_median_m', 'method')},
                          **(v or {})))
    export.write_csv(table, os.path.join(out_dir, 'compare.csv'))
    return table
