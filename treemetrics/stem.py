"""Measured stem for RayExtract trees.

Thin horizontal slices every 50 cm from 0.8 m above the (uphill) ground, one RANSAC circle each, a smooth taper
curve through them; this curve replaces RayExtract's stem radii, with a cone from the last circle to the top.
Tested on 65 felled trees (Demol et al. 2021): 64 of 65 within +-25 % of weighed wood volume."""
from collections import defaultdict
import numpy as np
import treevol as TV


def smooth_stem(R):
    """Stems taper gradually: replace the single circle radii by a smooth taper curve through all of them
    (log-radius quadratic in height, outliers dropped), never growing upwards; the stem axis is smoothed too."""
    if len(R) < 3: return R
    z = R[:, 2]; y = np.log(R[:, 3]); deg = 2 if len(R) >= 6 else 1
    keep = np.ones(len(R), bool)
    for _ in range(3):
        c = np.polyfit(z[keep], y[keep], deg); res = y - np.polyval(c, z)
        sd = max(res[keep].std(), 0.03); new = np.abs(res) < 2.5 * sd
        if (new == keep).all(): break
        keep = new
    r = np.exp(np.polyval(c, z))
    R = R.copy()
    for j in (0, 1):                                           # gently smoothed stem axis
        cj = np.polyfit(z, R[:, j], min(3, len(R) - 1)); R[:, j] = np.polyval(cj, z)
    # a leaning stem cut horizontally gives an ellipse: correct the radius by the stem's tilt
    d = np.gradient(R[:, :3], z, axis=0); cos_t = 1 / np.sqrt(1 + d[:, 0] ** 2 + d[:, 1] ** 2)
    R[:, 3] = np.minimum.accumulate(r * cos_t)                # continuous, never thicker further up
    return R


def stem_profile(P, g, step=0.5):
    """Horizontal RANSAC circles every 50 cm from 0.8 m above ground g, while they stay plausible; returns the smoothed taper curve (x, y, z, r)."""
    rows, last, fails = [], None, 0
    for z in np.arange(g + 0.8, g + (P[:, 2].max() - g) * 0.95, step):
        s = P[(P[:, 2] > z - 0.05) & (P[:, 2] < z + 0.05)]
        if last is not None: s = s[np.hypot(s[:, 0] - last[0], s[:, 1] - last[1]) < max(0.4, 2.5 * last[2])]
        ok = len(s) >= 30
        if ok:
            cx, cy, r, inl, gap, _ = TV._ransac(s[:, :2])
            ok = not (inl < 0.4 or gap > 200 or r > 0.8 or r < 0.02 or (last is not None and r > last[2] * 1.25 + 0.01))
        if not ok:
            fails += 1
            if fails >= 2 and last is not None: break
            continue
        fails = 0; rows.append((cx, cy, z, r)); last = (cx, cy, r)
    R = np.array(rows) if rows else np.zeros((0, 4))
    return smooth_stem(R) if len(R) >= 3 else np.zeros((0, 4))


def replace_stem(S, R, mode='cone'):
    """Like stem.replace_stem, but with a choice of what happens above the last measured circle:
    cone = straight cone from the last circle to the top; rx = RayExtract's own radii, scaled to match at the last circle.
    Main stem = from the base always follow the child with the most branch length above it; never thicker upwards;
    side branches never thicker than the stem."""
    S = S.copy(); par = S[:, 4].astype(int); n = len(S); kids = defaultdict(list)
    for i, p in enumerate(par):
        if p >= 0: kids[p].append(i)
    seg = np.where(par >= 0, np.linalg.norm(S[:, :3] - S[np.maximum(par, 0), :3], axis=1), 0)
    depth = np.zeros(n, int); order = [i for i in range(n) if par[i] < 0]
    for u in order:
        for v in kids[u]: depth[v] = depth[u] + 1; order.append(v)
    L = seg.copy()
    for i in np.argsort(-depth):
        if par[i] >= 0: L[par[i]] += L[i]
    roots = [i for i in range(n) if par[i] < 0]; chain = [max(roots, key=lambda r: L[r])]
    while kids[chain[-1]]: chain.append(max(kids[chain[-1]], key=lambda c: L[c]))
    chain = np.array(chain)
    if len(R) < 3: return S
    z = S[chain, 2]; inside = z <= R[-1, 2]; above = chain[~inside]
    r_rx = S[chain, 3].copy()
    S[chain[inside], 3] = np.interp(z[inside], R[:, 2], R[:, 3])
    if len(above):
        H = S[chain, 2].max(); zt = R[-1, 2]; rt = R[-1, 3]
        cone = np.maximum(rt * np.clip((H - S[above, 2]) / max(H - zt, 0.5), 0, 1), 0.003)
        if mode == 'cone': S[above, 3] = cone
        elif mode == 'rx':
            j = np.argmin(np.abs(z - zt)); k = rt / max(np.median(r_rx[max(j - 2, 0): j + 3]), 1e-3)
            S[above, 3] = np.maximum(r_rx[~inside] * k, 0.003)
    S[chain, 3] = np.minimum.accumulate(S[chain, 3])
    onstem = set(chain.tolist())
    for c in chain:
        for k in kids[c]:
            if k not in onstem and S[k, 3] > S[c, 3]:
                f = S[c, 3] / S[k, 3]; st = [k]
                while st:
                    u = st.pop(); S[u, 3] *= f; st += kids[u]
    return S


def stem_volume_model(S, R, P_height, g):
    """Overlap-free volume [m^3] and the model to draw. When the stem could be measured to at least half the tree
    height, the cone version is used; below that, the mean of the cone and the RayExtract-taper version.
    Tested on the felled trees by hiding the upper stem: with the cone alone the corrected volume drops by 5, 11
    and 18 % when the stem is visible only to 50, 30 and 15 % of the height; with this rule by 4, 5 and 3 %.
    Fully measured stems are unchanged (6.2 vs 6.1 % mean deviation)."""
    A = replace_stem(S, R, 'cone'); va = np.subtract(*TV.volume_no_overlap(A))
    frac = (R[-1, 2] - g) / P_height if len(R) >= 3 else 0.0
    if frac >= 0.5: return va, A, frac
    vb = np.subtract(*TV.volume_no_overlap(replace_stem(S, R, 'rx')))
    return (va + vb) / 2, A, frac
