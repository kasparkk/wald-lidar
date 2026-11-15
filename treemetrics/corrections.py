"""Corrections fitted against 65 felled and weighed trees (Demol, Gielen & Verbeeck 2021, Belgium:
beech, Scots pine, ash, larch; leaf-off terrestrial scans). Only tested on this one data set, so treat the
constants as a starting point and refit them when you have your own felled trees (see validate_felled.py).

Volume (leave-one-out on the 65 trees, mean absolute deviation per tree):
  RayExtract + measured stem, uncorrected            10.6 %
  + twig term for conifers, broadleaf x 0.94          6.3 %   (65/65 within +-25 %, 52/65 within +-10 %;
                                                              6.1 % before the short-stem rule in stem.py)
  classic form factor x DBH^2 x height, for compare   9.8 % with tape DBH, 11.3 % with point DBH

DBH (points vs tape, 62 trees passing the circle/girth check): 0.60 cm, with the pine bark offset 0.44 cm.
"""
import numpy as np

# Missing fine twigs (pipe model). Every branch tip in the model is cut at some radius r; the twig crown
# above it has volume ~ r^2 * length with length ~ r, so ~ r^3. Fitted on the 35 conifers (no intercept):
#   added volume [m^3] = TWIG_COEF * sum(r_tip^3) [r in m]
# r^3 clearly beat r^2, tip count and length of thin branches. For broadleaf the fitted term came out
# negative (the model is already slightly high there), so broadleaf gets a plain factor instead.
TWIG_COEF = 2656.0
BROADLEAF_FACTOR = 0.94
MAX_TWIG_SHARE = 0.25
CONIFER_GENERA = ('pinus', 'larix', 'picea', 'abies', 'pseudotsuga', 'taxus', 'thuja', 'juniperus', 'tsuga', 'cedrus')

# Point DBH vs tape: beech, ash, larch within 0.25 cm; pine about 0.9 cm too small (median), probably
# because the tape spans the furrowed bark while the points follow the grooves. Unverified elsewhere.
DBH_OFFSET_CM = {'pinus sylvestris': 0.9}


def is_conifer(species):
    s = (species or '').strip().lower()
    return s in ('conifer', 'coniferous', 'nadelholz') or s.split(' ')[0] in CONIFER_GENERA


TIP_R_MAX = 0.05   # m. Thicker "tips" are model artefacts (a stem piece ending mid-tree), not cut twigs.
# The Belgian models have no tip above 2.2 cm, so the cap does not change the fit; on Schorfheide pines three
# trees had one 12-16 cm "tip" each that would otherwise have added 460-625 %.


def tip_r3(S, r_max=TIP_R_MAX):
    """Sum of r^3 over the branch tips of a model S (rows x, y, z, radius, parent; parent -1 = root),
    ignoring tips thicker than r_max."""
    par = S[:, 4].astype(int)
    has_kid = np.zeros(len(S), bool); has_kid[par[par >= 0]] = True
    r = S[~has_kid & (par >= 0), 3]
    return float((r[r <= r_max] ** 3).sum())


def correct_volume(volume_m3, S, species):
    """Corrected wood volume [m^3] of one tree: conifers + twig term, broadleaf x 0.94.
    volume_m3: overlap-free volume of the RayExtract + measured stem model S (treevol.volume_no_overlap)."""
    if is_conifer(species):
        # capped at the largest addition seen in the calibration data (+25 %): beyond that the tip radii lie
        # outside the fitted range (Belgian tips <= 2.2 cm) and the r^3 term would extrapolate
        return volume_m3 + min(TWIG_COEF * tip_r3(S), MAX_TWIG_SHARE * volume_m3)
    return volume_m3 * BROADLEAF_FACTOR


def dbh_checked(circle_cm, girth_cm, ok, max_diff_cm=None):
    """DBH as the mean of circle and girth, or None when the measurement is not plausible.
    ok comes from treevol.dbh (enough inliers, stem covered, radius sane). This check caught a 13 cm error
    (a fork at 1.3 m) on the felled trees. A check against a DBH-height relation does NOT work for ground
    scans: it flagged 13 good trees and missed that one."""
    if not ok or not np.isfinite(circle_cm):
        return None
    if max_diff_cm is not None and abs(circle_cm - girth_cm) > max_diff_cm:
        return None
    return (circle_cm + girth_cm) / 2


def correct_dbh(dbh_cm, species):
    if dbh_cm is None:
        return None
    return dbh_cm + DBH_OFFSET_CM.get((species or '').strip().lower(), 0.0)
