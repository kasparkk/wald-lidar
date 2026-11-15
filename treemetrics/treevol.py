"""Shared helpers: wood volume without overlaps at forks, ground height and DBH at 1.3 m."""
import numpy as np
from collections import defaultdict
from scipy.spatial import ConvexHull


def _seg_vol(S, i, frustum, side=()):
    p = int(S[i, 4]); L = np.linalg.norm(S[i, :3] - S[p, :3])
    if frustum:
        # a side branch starts with its own radius, not with the radius of the branch it leaves
        a = S[i, 3] if i in side else S[p, 3]; b = S[i, 3]
        return np.pi * L / 3 * (a * a + a * b + b * b)
    return np.pi * S[i, 3] ** 2 * L


def _dist_to_line(x, o, d):
    v = x - o; return np.linalg.norm(v - np.dot(v, d) * d)


def volume_no_overlap(S, frustum=False):
    """S rows: x, y, z, radius, parent (-1 = root). Each row i is the piece from its parent to i.
    Total volume of all pieces minus the parts of side branches that lie inside the branch they leave:
    at every fork the thickest child continues the host; every other child is followed outward and the
    pieces (or piece fractions) closer to the host axis than the host radius are subtracted."""
    par = S[:, 4].astype(int); n = len(S); kids = defaultdict(list)
    for i, p in enumerate(par):
        if p >= 0: kids[p].append(i)
    side = {c for p, ks in kids.items() if len(ks) > 1 for c in ks if c != max(ks, key=lambda q: S[q, 3])}
    total = sum(_seg_vol(S, i, frustum, side) for i in range(n) if par[i] >= 0)
    overlap = 0.0
    for p, ks in list(kids.items()):
        if len(ks) < 2: continue
        host = max(ks, key=lambda c: S[c, 3]); rh = S[host, 3]
        d = S[host, :3] - S[p, :3]
        if np.linalg.norm(d) < 1e-6 and par[p] >= 0: d = S[p, :3] - S[par[p], :3]
        if np.linalg.norm(d) < 1e-6: continue
        d = d / np.linalg.norm(d); o = S[p, :3]
        for c in ks:
            if c == host: continue
            prev, cur, steps = p, c, 0
            while steps < 400:
                de = _dist_to_line(S[cur, :3], o, d); ds = _dist_to_line(S[prev, :3], o, d)
                v = _seg_vol(S, cur, frustum, side)
                if de <= rh: overlap += v
                else:
                    overlap += v * float(np.clip((rh - ds) / max(de - ds, 1e-9), 0, 1)); break
                if not kids.get(cur): break
                prev, cur = cur, max(kids[cur], key=lambda q: S[q, 3]); steps += 1
    return total, overlap


def ground_height(P):
    """Ground at the stem base for a single-tree cloud without ground points: lowest stem points per
    30-degree sector around the stem; on slopes breast height is measured from the uphill side, so the
    ground is the upper quartile of the sector minima (equals the lowest point on flat ground)."""
    zmin = np.percentile(P[:, 2], 0.1)
    s = P[(P[:, 2] > zmin + 0.4) & (P[:, 2] < zmin + 0.7)]
    if len(s) < 30: return zmin
    c = np.median(s[:, :2], axis=0)
    r = np.percentile(np.hypot(*(s[:, :2] - c).T), 90)
    b = P[(P[:, 2] < zmin + 0.6) & (np.hypot(*(P[:, :2] - c).T) < r + 0.15)]
    ang = np.floor((np.arctan2(b[:, 1] - c[1], b[:, 0] - c[0]) + np.pi) / (np.pi / 6)).astype(int) % 12
    mins = np.array([b[ang == k, 2].min() for k in range(12) if (ang == k).sum() >= 5])
    return float(np.percentile(mins, 75)) if len(mins) >= 4 else zmin


def _ransac(xy, tol=0.015, it=300):
    x, y = xy[:, 0], xy[:, 1]; A = np.c_[2 * x, 2 * y, np.ones(len(x))]; bb = x * x + y * y
    rng = np.random.default_rng(0); best = None
    for _ in range(it):
        i = rng.choice(len(x), 3, replace=False)
        try: c = np.linalg.solve(A[i], bb[i])
        except np.linalg.LinAlgError: continue
        r = np.sqrt(max(c[2] + c[0] ** 2 + c[1] ** 2, 0)); inl = np.abs(np.hypot(x - c[0], y - c[1]) - r) < tol
        if best is None or inl.sum() > best.sum(): best = inl
    c = np.linalg.lstsq(A[best], bb[best], rcond=None)[0]
    r = np.sqrt(max(c[2] + c[0] ** 2 + c[1] ** 2, 0))
    ang = np.sort(np.arctan2(y[best] - c[1], x[best] - c[0]))
    gap = np.degrees(np.max(np.diff(np.r_[ang, ang[0] + 2 * np.pi]))) if best.sum() > 3 else 360
    return c[0], c[1], r, best.mean(), gap, best


def dbh(P, ground=None, h=1.3, half=0.05):
    """DBH at h above the ground (uphill side). Returns (circle_dbh, girth_dbh, ok):
    circle = RANSAC circle; girth = convex-hull perimeter / pi of the circle's inliers, like a tape
    around the bark (only where at least 300 degrees of the stem are covered, else the circle value)."""
    g = ground_height(P) if ground is None else ground
    s = P[(P[:, 2] > g + h - half) & (P[:, 2] < g + h + half)]
    if len(s) < 30: return np.nan, np.nan, False
    cx, cy, r, inl, gap, best = _ransac(s[:, :2])
    ok = inl >= 0.5 and gap <= 200 and r <= 0.8
    girth = 2 * r
    pts = s[best, :2]
    if gap <= 60 and len(pts) >= 10:
        girth = ConvexHull(pts).area / np.pi          # in 2D, .area is the perimeter
    return 2 * r, girth, ok


# ---- fallback when the stem is not cleanly captured at 1.3 m (branches in the way, one-sided scan)
# Tested on the 62 clean Belgian trees: ignoring the 1.2-1.4 m slices, alt_dbh is as accurate as the direct
# value (0.58 vs 0.60 cm off tape). With half of each stem removed, half_dbh is 1.3 cm off and 0.8 cm too small
# on average (HALF_BIAS_CM is added in dbh_with_fallback).
HALF_BIAS_CM = 0.8


def alt_dbh(P, g, skip=(1.2, 1.4)):
    """DBH from slices above/below 1.3 m: circle+girth at 0.9-2.5 m, robust line through the plausible ones, value at 1.3."""
    H, D = [], []
    for h in np.arange(0.9, 2.51, 0.1):
        if skip and skip[0] <= h <= skip[1] + 1e-9: continue
        c, gi, ok = dbh(P, g, h=h)
        if ok and np.isfinite(c): H.append(h); D.append((c + gi) / 2 * 100)
    H, D = np.array(H), np.array(D)
    if len(H) < 3 or np.min(np.abs(H - 1.3)) > 0.6: return None, len(H)
    keep = np.ones(len(H), bool)
    for _ in range(3):
        b = np.polyfit(H[keep], D[keep], 1); res = D - np.polyval(b, H)
        s = max(1.4826 * np.median(np.abs(res[keep])), 0.3); new = np.abs(res) < 3 * s
        if new.sum() < 3 or (new == keep).all(): break
        keep = new
    if b[0] > 1.0: return None, len(H)        # thicker upwards by >1 cm per m: fork/branch, not a stem
    return float(np.polyval(b, 1.3)), int(keep.sum())


def half_dbh(P, g, lo=0.6, hi=2.5, min_inl=0.85, max_gap=250):
    """For stems scanned from one side only: circle fits (no girth) with clean inliers but up to 250 deg gap."""
    H, D = [], []
    for h in np.arange(lo, hi + 1e-9, 0.1):
        s = P[(P[:, 2] > g + h - .05) & (P[:, 2] < g + h + .05)]
        if len(s) < 30: continue
        cx, cy, r, inl, gap, best = _ransac(s[:, :2])
        if inl >= min_inl and gap <= max_gap and 0.02 < r <= 0.8: H.append(h); D.append(2 * r * 100)
    H, D = np.array(H), np.array(D)
    if len(H) < 3 or np.min(np.abs(H - 1.3)) > 0.6: return None, len(H)
    keep = np.ones(len(H), bool)
    for _ in range(3):
        b = np.polyfit(H[keep], D[keep], 1); res = D - np.polyval(b, H)
        s = max(1.4826 * np.median(np.abs(res[keep])), 0.3); new = np.abs(res) < 3 * s
        if new.sum() < 3 or (new == keep).all(): break
        keep = new
    if b[0] > 1.0: return None, len(H)
    return float(np.polyval(b, 1.3)), int(keep.sum())


def dbh_with_fallback(P, g):
    """DBH [cm] and how it was obtained: 'direct' (circle+girth at 1.3 m), 'slices' (slices 0.9-2.5 m projected
    to 1.3 m), 'one-sided' (half-scanned stem, circle only, bias-corrected) or (None, 'not measurable')."""
    c, gi, ok = dbh(P, g)
    if ok and np.isfinite(c): return (c + gi) / 2 * 100, 'direct'
    d, _ = alt_dbh(P, g, skip=None)
    if d is not None and 3 < d < 150: return d, 'slices'
    d, _ = half_dbh(P, g)
    if d is not None and 3 < d < 150: return d + HALF_BIAS_CM, 'one-sided'
    return None, 'not measurable'
