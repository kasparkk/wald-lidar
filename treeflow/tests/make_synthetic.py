"""Synthetic test stands with known trees, so the pipeline can be scored before real reference data exists.
Writes a dense ground-scan-like file (no ground class) and an airborne-like file (with ground class),
plus reference.csv with the true positions, heights and DBH."""
import sys, csv
import numpy as np
import laspy

rng = np.random.default_rng(42)
S = 50.0
out = sys.argv[1] if len(sys.argv) > 1 else '.'


def ground_z(x, y):
    return 40 + 0.04 * x + 0.02 * y + 0.4 * np.sin(x / 7) * np.cos(y / 9)


trees = []
while len(trees) < 45:
    x, y = rng.uniform(3, S - 3, 2)
    if all(np.hypot(x - t['x'], y - t['y']) > 4.5 for t in trees):
        h = rng.uniform(14, 30)
        trees.append(dict(x=x, y=y, h=h, dbh=max(8, 1.4 * h - 4 + rng.normal(0, 3)),
                          r=np.clip(0.5 + 0.09 * h + rng.normal(0, 0.4), 1.5, 4.5), cb=h * rng.uniform(0.45, 0.65)))
# a few understorey trees
for _ in range(8):
    x, y = rng.uniform(3, S - 3, 2)
    h = rng.uniform(4, 8)
    trees.append(dict(x=x, y=y, h=h, dbh=max(4, 1.2 * h), r=1.0, cb=h * 0.3, under=True))


def tree_points(t, n_stem, n_crown):
    zs = rng.uniform(0, t['cb'], n_stem)
    a = rng.uniform(0, 2 * np.pi, n_stem)
    rr = t['dbh'] / 200 * (1 - 0.3 * zs / t['h']) + rng.normal(0, 0.004, n_stem)
    stem = np.c_[t['x'] + rr * np.cos(a), t['y'] + rr * np.sin(a), zs]
    u = rng.uniform(-1, 1, n_crown); th = rng.uniform(0, 2 * np.pi, n_crown); k = rng.uniform(0.6, 1, n_crown) ** 0.5
    ry = (t['h'] - t['cb']) / 2
    crown = np.c_[t['x'] + np.sqrt(1 - u * u) * np.cos(th) * t['r'] * k,
                  t['y'] + np.sqrt(1 - u * u) * np.sin(th) * t['r'] * k,
                  t['cb'] + ry + u * ry * k]
    return stem, crown


def write(path, P, cls, ret):
    h = laspy.LasHeader(point_format=1, version='1.4'); h.scales = [0.001] * 3; h.offsets = P.min(0)
    las = laspy.LasData(h); las.x, las.y, las.z = P.T
    las.classification = cls; las.return_number = ret; las.number_of_returns = np.maximum(ret, 1)
    las.write(path)


# dense ground scan: everything visible, ground unclassified
pts, ret = [], []
gx, gy = rng.uniform(0, S, (2, 250000))
pts.append(np.c_[gx, gy, ground_z(gx, gy) + rng.normal(0, 0.01, len(gx))])
for t in trees:
    st, cr = tree_points(t, 3000, 6000)
    for P in (st, cr):
        P[:, 2] += ground_z(t['x'], t['y'])
        pts.append(P)
P = np.vstack(pts)
write(f'{out}/synthetic_ground_scan.laz', P, np.zeros(len(P), np.uint8), np.ones(len(P), np.uint8))

# airborne-like: ~25 pulses/m² from above; first return = highest surface in a 0.2 m cell,
# plus some later returns through gaps; ground classified
cell = 0.2
key = lambda X: (np.floor(X[:, 0] / cell).astype(int) * 10000 + np.floor(X[:, 1] / cell).astype(int))
k = key(P); order = np.lexsort((-P[:, 2], k)); first = np.r_[True, k[order][1:] != k[order][:-1]]
top = P[order[first]]
top = top[rng.random(len(top)) < 0.9]                          # ~22 pts/m²
later = P[rng.random(len(P)) < 0.004]
A = np.vstack([top, later])
gz = ground_z(A[:, 0], A[:, 1])
cls = np.where(A[:, 2] - gz < 0.15, 2, 1).astype(np.uint8)
r = np.r_[np.ones(len(top)), np.full(len(later), 2)].astype(np.uint8)
write(f'{out}/synthetic_airborne.laz', A, cls, r)

with open(f'{out}/reference.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['x', 'y', 'height_m', 'dbh_cm', 'understorey'])
    for t in trees:
        w.writerow([round(t['x'], 2), round(t['y'], 2), round(t['h'], 2), round(t['dbh'], 1), int(t.get('under', False))])
print(len(trees), 'trees;', len(P), 'ground-scan points;', len(A), 'airborne points')
