"""Use lidR (R package, the standard tool for airborne tree segmentation) instead of treeflow's own
detector. treeflow keeps doing everything around it: reading, ground model, measurements, photo values,
checking and output. Needs R with lidR installed (Rscript on the PATH)."""
import os, shutil, subprocess, tempfile
import numpy as np
import laspy

SCRIPT = os.path.join(os.path.dirname(__file__), 'lidr_segment.R')


def available():
    if not shutil.which('Rscript'):
        return False
    r = subprocess.run(['Rscript', '-e', 'cat(requireNamespace("lidR", quietly=TRUE))'],
                       capture_output=True, text=True)
    return r.stdout.strip().endswith('TRUE')


def segment(c, hag, ground, chm_res, cfg):
    algo = cfg.get('detector', 'lidr-dalponte').split('-', 1)[1]
    tmp = tempfile.mkdtemp(prefix='treeflow_lidr_')
    try:
        h = laspy.LasHeader(point_format=1, version='1.4')
        h.scales = [0.01, 0.01, 0.01]; h.offsets = [c.x.min(), c.y.min(), c.z.min()]
        h.add_extra_dim(laspy.ExtraBytesParams(name='pid', type=np.uint32))
        las = laspy.LasData(h)
        las.x, las.y, las.z = c.x, c.y, c.z
        las.classification = np.where(ground, 2, 1).astype(np.uint8)
        las.return_number = np.clip(c.ret, 1, 7); las.number_of_returns = np.clip(np.maximum(c.nret, c.ret), 1, 7)
        las.pid = np.arange(len(c), dtype=np.uint32)
        inp = os.path.join(tmp, 'in.las'); las.write(inp)
        args = ['Rscript', SCRIPT, inp, tmp, algo, str(cfg.get('lidr_ws_a', 3.0)), str(cfg.get('lidr_ws_b', 0.07)),
                str(cfg.get('min_tree_height', 5.0)), str(chm_res)]
        r = subprocess.run(args, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError('lidR failed: ' + (r.stderr or r.stdout)[-800:])
        pts = np.loadtxt(os.path.join(tmp, 'lidr_points.csv'), delimiter=',', skiprows=1, dtype=np.int64, ndmin=2)
        tree_id = np.zeros(len(c), dtype=np.int64)
        tree_id[pts[:, 0]] = pts[:, 1]
        tops_path = os.path.join(tmp, 'lidr_tops.csv')
        tops = np.loadtxt(tops_path, delimiter=',', skiprows=1, ndmin=2) if os.path.exists(tops_path) else None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    tree_id[ground | (hag < cfg.get('min_point_height', 1.0))] = 0
    seeds = []
    ids = np.unique(tree_id[tree_id > 0])
    if tops is not None and len(tops):
        top_xy = {int(t[0]): (t[1], t[2]) for t in tops}
    else:
        top_xy = {}
    for k in ids:
        if k in top_xy:
            x, y = top_xy[k]
        else:                       # li2012 has no tops: use the highest point of the tree
            m = np.nonzero(tree_id == k)[0]; j = m[np.argmax(hag[m])]; x, y = c.x[j], c.y[j]
        seeds.append(dict(id=int(k), x=float(x), y=float(y), stem=None))
    return tree_id, seeds
