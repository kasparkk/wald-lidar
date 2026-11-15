"""Tree model, wood volume, DBH and height for every segmented tree in a folder (one LAS/LAZ file per tree,
ground or backpack scan). Needs RayCloudTools (rayimport, rayextract) built on Linux/WSL.

  python3 run_trees.py TREE_FOLDER --out results [--bin PATH/raycloudtools/build/raycloudtools]
                       [--species species.csv]       # columns: tree,species (e.g. "Pinus sylvestris")

Writes results/trees.csv and one model file per tree (results/models/<tree>.txt: x, y, z, radius, parent;
parent -1 = root; units metres). Columns:
  height_m            top minus ground at the stem base
  dbh_cm              mean of circle and girth at 1.3 m above the uphill ground; if not cleanly captured there,
                      projected from slices 0.6-2.5 m (see dbh_method); empty if not measurable
  dbh_method          direct | slices | one-sided | estimated from height and crown (every tree gets a DBH)
  dbh_corrected_cm    + bark offset where known (pine +0.9 cm)
  volume_l            RayExtract + measured stem, overlap-free, uncorrected
  volume_corrected_l  + twig term for conifers, x 0.94 for broadleaf (needs species; empty otherwise)
  volume_rayextract_l plain RayExtract, for comparison
  stem_measured_share how far up the stem could be measured (share of tree height); below 0.5 the volume
                      averages a cone and RayExtract's taper above the last circle (see stem.stem_volume_model)
  x, y, ground_z      stem position (centre of the lowest stem circle, else of the lowest metre of points) and
                      ground height, in the coordinates of the input files
  lon, lat            the same position in WGS84 degrees, if the files carry a coordinate system or --crs is
                      given (e.g. --crs 25832 for UTM 32N); empty otherwise. Plots in local coordinates need
                      --offset X,Y (plot origin in that CRS) first.
"""
import argparse, csv, glob, os, subprocess, tempfile, time
import numpy as np, laspy
import treevol as TV, stem, corrections as C


def tool(binroot, name):
    if binroot:
        for c in glob.glob(os.path.join(binroot, name, name + '*')):
            if os.access(c, os.X_OK): return c
    return name   # on PATH


def write_inputs(P, z0, ply, ground):
    with open(ply, 'wb') as f:
        f.write(f'ply\nformat binary_little_endian 1.0\nelement vertex {len(P)}\nproperty float x\nproperty float y\nproperty float z\nend_header\n'.encode())
        f.write(P.astype('<f4').tobytes())
    # flat 25 cm ground mesh under the tree: RayExtract grows trees from it (single-tree files have no ground)
    xs = np.arange(P[:, 0].min(), P[:, 0].max() + 0.25, 0.25); ys = np.arange(P[:, 1].min(), P[:, 1].max() + 0.25, 0.25)
    gx, gy = np.meshgrid(xs, ys); nx, ny = len(xs), len(ys)
    V = np.c_[gx.ravel(), gy.ravel(), np.full(gx.size, z0)]
    F = [(j * nx + i, j * nx + i + 1, j * nx + i + nx + 1) for j in range(ny - 1) for i in range(nx - 1)]
    F += [(j * nx + i, j * nx + i + nx + 1, j * nx + i + nx) for j in range(ny - 1) for i in range(nx - 1)]
    with open(ground, 'wb') as f:
        f.write(f'ply\nformat binary_little_endian 1.0\nelement vertex {len(V)}\nproperty float x\nproperty float y\nproperty float z\nelement face {len(F)}\nproperty list uchar int vertex_indices\nend_header\n'.encode())
        f.write(V.astype('<f4').tobytes())
        rec = np.zeros(len(F), dtype=[('n', 'u1'), ('i', '<i4', 3)]); rec['n'] = 3; rec['i'] = F; f.write(rec.tobytes())


def read_trees(path):
    lines = [l for l in open(path) if l.strip() and not l.startswith('#')]
    n = len([h for h in lines[0].strip().split(',') if h.strip()])
    return [np.array([float(x) for x in l.replace(' ', '').strip().strip(',').split(',')]).reshape(-1, n) for l in lines[1:]]


def rayextract(P, binroot, crop=0.25):
    """Largest RayExtract tree for one point cloud, rows x, y, z, radius, parent. No scanner positions are
    needed (rays are faked from above). crop 0.25 m instead of the default 1 m keeps more crown wood."""
    z0 = float(np.percentile(P[:, 2], 0.1))
    with tempfile.TemporaryDirectory() as d:
        ply, ground = os.path.join(d, 'tree.ply'), os.path.join(d, 'ground.ply')
        write_inputs(P, z0, ply, ground)
        subprocess.run([tool(binroot, 'rayimport'), ply, 'ray', '0,0,-10', '--max_intensity', '0'], cwd=d, capture_output=True)
        subprocess.run([tool(binroot, 'rayextract'), 'trees', os.path.join(d, 'tree_raycloud.ply'), ground, '--crop_length', str(crop)], cwd=d, capture_output=True)
        tf = os.path.join(d, 'tree_raycloud_trees.txt')
        T = read_trees(tf) if os.path.exists(tf) else []
    return max(T, key=len)[:, :5] if T else None


def measure(P, species='', binroot=''):
    g = TV.ground_height(P)
    r = dict(height_m=round(float(P[:, 2].max() - g), 2))
    from scipy.spatial import ConvexHull
    r['crown_d_m'] = round(float(2 * np.sqrt(ConvexHull(P[:, :2]).volume / np.pi)), 2)
    d, how = TV.dbh_with_fallback(P, g)
    r['dbh_method'] = how
    r['dbh_cm'] = round(float(d), 1) if d else ''
    dc = C.correct_dbh(d, species); r['dbh_corrected_cm'] = round(float(dc), 1) if dc else ''
    R = stem.stem_profile(P, g)
    if len(R): xy = R[0, :2]
    else: xy = np.median(P[P[:, 2] < g + 1.0, :2], axis=0)
    r.update(x=round(float(xy[0]), 2), y=round(float(xy[1]), 2), ground_z=round(float(g), 2))
    S0 = rayextract(P, binroot); S = None
    if S0 is not None:
        v, S, frac = stem.stem_volume_model(S0, R, float(P[:, 2].max() - g), g)
        r['stem_measured_share'] = round(float(frac), 2)
        tot0, ov0 = TV.volume_no_overlap(S0)
        r.update(volume_l=round(float(v) * 1000, 1), volume_rayextract_l=round(float(tot0 - ov0) * 1000, 1),
                 volume_corrected_l=round(float(C.correct_volume(v, S, species)) * 1000, 1) if species else '')
    return r, S


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input'); ap.add_argument('--out', required=True); ap.add_argument('--bin', default='')
    ap.add_argument('--species', help='CSV with columns tree,species')
    ap.add_argument('--crs', help='EPSG code of the input coordinates, if the files do not carry one')
    ap.add_argument('--offset', default='0,0', help='X,Y added to local plot coordinates before conversion')
    a = ap.parse_args(); os.makedirs(os.path.join(a.out, 'models'), exist_ok=True)
    sp = {}
    if a.species:
        with open(a.species) as fh: sp = {row['tree']: row['species'] for row in csv.DictReader(fh)}
    rows = []; crs = None; off = np.array([float(v) for v in a.offset.split(',')])
    for f in sorted(glob.glob(os.path.join(a.input, '*.la[sz]'))):
        name = os.path.splitext(os.path.basename(f))[0]; t = time.time()
        l = laspy.read(f); P = np.c_[l.x, l.y, l.z]
        if crs is None:
            try: crs = l.header.parse_crs()
            except Exception: crs = None
        r, S = measure(P, sp.get(name, ''), a.bin)
        if S is not None: np.savetxt(os.path.join(a.out, 'models', name + '.txt'), S, fmt='%.4f', header='x y z radius parent')
        rows.append(dict(tree=name, species=sp.get(name, ''), seconds=round(time.time() - t, 1), **r)); print(rows[-1], flush=True)
    # every tree gets a DBH: where none could be measured, estimate it from height and crown width (fitted per
    # species group on the measured trees of this run, see drone_dbh.fill_missing)
    import pandas as pd, drone_dbh
    T = pd.DataFrame(rows)
    if len(T) and T.dbh_cm.replace('', np.nan).notna().sum() >= 3:
        T['dbh'] = pd.to_numeric(T.dbh_cm.replace('', np.nan)); T['keep'] = T.dbh.notna()
        T['species'] = T.species.fillna(''); T['height'] = T.height_m; T['crown_d'] = T.crown_d_m
        F = drone_dbh.fill_missing(T)
        for i, row in F.iterrows():
            if not row.keep and np.isfinite(row.dbh_final):
                rows[i]['dbh_cm'] = round(float(row.dbh_final), 1); rows[i]['dbh_method'] = 'estimated from height and crown'
                dc = C.correct_dbh(row.dbh_final, rows[i]['species']); rows[i]['dbh_corrected_cm'] = round(float(dc), 1)
    # geographic position: input CRS (from --crs or the file header) -> WGS84
    src = f'EPSG:{a.crs}' if a.crs else crs
    if src is not None and rows:
        from pyproj import Transformer
        tr = Transformer.from_crs(src, 'EPSG:4326', always_xy=True)
        for row in rows:
            lon, lat = tr.transform(row['x'] + off[0], row['y'] + off[1])
            row['lon'], row['lat'] = round(lon, 7), round(lat, 7)
    cols = ['tree', 'species', 'x', 'y', 'lon', 'lat', 'ground_z', 'height_m', 'dbh_cm', 'dbh_method', 'dbh_corrected_cm', 'crown_d_m', 'volume_l', 'volume_corrected_l', 'volume_rayextract_l', 'stem_measured_share', 'seconds']
    with open(os.path.join(a.out, 'trees.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore', restval=''); w.writeheader(); w.writerows(rows)


if __name__ == '__main__':
    main()
