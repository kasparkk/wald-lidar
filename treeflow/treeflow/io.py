"""Reading point clouds and describing what kind of scan they are."""
from dataclasses import dataclass, field
import numpy as np
import laspy


@dataclass
class Cloud:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    cls: np.ndarray            # LAS classification (2 = ground), zeros if absent
    ret: np.ndarray            # return number, ones if absent
    nret: np.ndarray           # number of returns, ones if absent
    intensity: np.ndarray
    crs: str | None = None
    meta: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.x)

    def subset(self, m):
        return Cloud(self.x[m], self.y[m], self.z[m], self.cls[m], self.ret[m],
                     self.nret[m], self.intensity[m], self.crs, dict(self.meta))


def read(path, bbox=None, buffer=0.0):
    """Read LAS/LAZ. bbox = (xmin, ymin, xmax, ymax) in file coordinates; buffer adds a margin
    so trees on the edge are segmented with their full crown."""
    las = laspy.read(path)
    x, y, z = np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)
    names = set(las.point_format.dimension_names)
    n = len(x)
    get = lambda name, default: np.asarray(las[name]) if name in names else np.full(n, default)
    c = Cloud(x, y, z,
              get('classification', 0).astype(np.uint8),
              get('return_number', 1).astype(np.uint8),
              get('number_of_returns', 1).astype(np.uint8),
              get('intensity', 0).astype(np.float32))
    try:
        crs = las.header.parse_crs()
        c.crs = crs.to_string() if crs else None
    except Exception:
        c.crs = None
    c.meta['source'] = str(path)
    if bbox is not None:
        x0, y0, x1, y1 = bbox
        m = (x >= x0 - buffer) & (x < x1 + buffer) & (y >= y0 - buffer) & (y < y1 + buffer)
        c = c.subset(m)
        c.meta['bbox'] = list(bbox)
        c.meta['buffer'] = buffer
    return c


def describe(c: Cloud):
    """Point density and a guess of the scanner type. Decides which detection method to use."""
    x0, x1, y0, y1 = c.x.min(), c.x.max(), c.y.min(), c.y.max()
    area = max((x1 - x0) * (y1 - y0), 1e-6)
    density = len(c) / area
    if density < 60:
        kind = 'airborne'          # ALS: crowns from above, few stems
    elif density < 1500:
        kind = 'drone'             # ULS: crowns detailed, some stems
    else:
        kind = 'ground'            # TLS / handheld: stems well captured
    has_ground = bool((c.cls == 2).sum() > 0.01 * len(c))
    multi = bool((c.ret > 1).any())
    return dict(points=int(len(c)), area_m2=float(area), density=float(density), kind=kind,
                has_ground_class=has_ground, multi_return=multi,
                extent=[float(x0), float(y0), float(x1), float(y1)])
