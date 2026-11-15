"""Ground: classify it if the file has no ground class, build a terrain model (DTM),
and turn every point's elevation into height above ground."""
import numpy as np
from scipy import ndimage as ndi
from scipy.interpolate import griddata


def _grid_index(x, y, x0, y0, res, nx, ny):
    j = np.clip(((x - x0) / res).astype(int), 0, nx - 1)
    i = np.clip(((y - y0) / res).astype(int), 0, ny - 1)
    return i, j


def classify_ground(c, res=1.0, window_m=12.0, tol=0.4):
    """Simple robust ground filter for files without class 2:
    lowest point per cell -> drop cells that sit clearly above a smoothed local minimum
    (tree trunks, low branches) -> points close to that surface are ground."""
    x0, y0 = c.x.min(), c.y.min()
    nx = int(np.ceil((c.x.max() - x0) / res)) + 1
    ny = int(np.ceil((c.y.max() - y0) / res)) + 1
    i, j = _grid_index(c.x, c.y, x0, y0, res, nx, ny)
    zmin = np.full((ny, nx), np.inf)
    np.minimum.at(zmin, (i, j), c.z)
    has = np.isfinite(zmin)
    filled = np.where(has, zmin, np.nanmax(np.where(has, zmin, np.nan)))
    w = max(3, int(window_m / res) | 1)
    opened = ndi.grey_opening(filled, size=(w, w))            # removes objects narrower than the window
    smooth = ndi.uniform_filter(opened, size=3)
    keep = has & (zmin - smooth < tol)
    yy, xx = np.nonzero(keep)
    surf = griddata((xx, yy), zmin[keep], np.meshgrid(np.arange(nx), np.arange(ny)), method='linear')
    nn = griddata((xx, yy), zmin[keep], np.meshgrid(np.arange(nx), np.arange(ny)), method='nearest')
    surf = np.where(np.isnan(surf), nn, surf)
    g = (c.z - surf[i, j]) < tol
    return g


def build_dtm(c, ground_mask, res=1.0):
    x0, y0 = c.x.min(), c.y.min()
    nx = int(np.ceil((c.x.max() - x0) / res)) + 1
    ny = int(np.ceil((c.y.max() - y0) / res)) + 1
    gx, gy, gz = c.x[ground_mask], c.y[ground_mask], c.z[ground_mask]
    # thin ground points to one per cell (lowest) to keep interpolation fast on dense clouds
    i, j = _grid_index(gx, gy, x0, y0, res / 2, nx * 2, ny * 2)
    key = i * (nx * 2) + j
    order = np.lexsort((gz, key))
    first = np.r_[True, key[order][1:] != key[order][:-1]]
    sel = order[first]
    gx, gy, gz = gx[sel], gy[sel], gz[sel]
    X, Y = np.meshgrid(x0 + (np.arange(nx) + 0.5) * res, y0 + (np.arange(ny) + 0.5) * res)
    dtm = griddata((gx, gy), gz, (X, Y), method='linear')
    nn = griddata((gx, gy), gz, (X, Y), method='nearest')
    dtm = np.where(np.isnan(dtm), nn, dtm)
    return dict(z=dtm, x0=float(x0), y0=float(y0), res=float(res), nx=nx, ny=ny)


def dtm_at(dtm, x, y):
    i, j = _grid_index(np.asarray(x), np.asarray(y), dtm['x0'], dtm['y0'], dtm['res'], dtm['nx'], dtm['ny'])
    return dtm['z'][i, j]


def normalise(c, cfg):
    """Returns (height_above_ground, ground_mask, dtm)."""
    if (c.cls == 2).sum() > 0.01 * len(c) and not cfg.get('force_ground_filter', False):
        g = c.cls == 2
        how = 'file classification'
    else:
        g = classify_ground(c, res=cfg.get('ground_res', 1.0))
        how = 'treeflow ground filter'
    dtm = build_dtm(c, g, res=cfg.get('dtm_res', 1.0))
    hag = c.z - dtm_at(dtm, c.x, c.y)
    return hag, g, dtm, how
