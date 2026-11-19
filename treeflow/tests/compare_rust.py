"""Compare the Python and Rust circle fit on the stems of a real scan: speed and difference per stem."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from sklearn.cluster import DBSCAN
from treeflow import io, ground, detect
from treeflow.pipeline import DEFAULTS

if detect._fit_circle_rust is None:
    sys.exit('treeflow_rust is not installed (pip install ./rust_kreis)')
c = io.read(sys.argv[1])
hag, g, _, _ = ground.normalise(c, dict(DEFAULTS))
m = (~g) & (hag >= 1.0) & (hag <= 1.6)
xy = np.c_[c.x[m], c.y[m]]
lab = DBSCAN(eps=0.1, min_samples=5).fit_predict(xy)
stems = [xy[lab == k] for k in set(lab) - {-1} if (lab == k).sum() >= 15]
res = {}
for name, f in [('Python', detect.fit_circle_python), ('Rust', detect._fit_circle_rust)]:
    t = time.perf_counter()
    res[name] = [f(p) for p in stems]
    print(f'{name}: {(time.perf_counter() - t) * 1000:.0f} ms for {len(stems)} stems')
both = [(a, b) for a, b in zip(res['Python'], res['Rust']) if a and b]
only_one = sum((a is None) != (b is None) for a, b in zip(res['Python'], res['Rust']))
dc = [np.hypot(a[0] - b[0], a[1] - b[1]) * 1000 for a, b in both]
dr = [abs(a[2] - b[2]) * 2000 for a, b in both]  # diameter difference in mm
print(f'circle found by both: {len(both)}, by only one: {only_one}')
if both:
    print(f'centre difference mm: median {np.median(dc):.1f}, max {max(dc):.1f}; '
          f'diameter difference mm: median {np.median(dr):.1f}, max {max(dr):.1f}')
