"""Conifer vs broadleaf from leaf-off laser scans of segmented trees: echo intensity plus crown shape.

Results (balanced accuracy, conifer vs broadleaf, leaf-off drone scans of two BioDiv-3DTrees plots):
  trained on                     SEW04_GD (pine/beech)   AEW33_GD (spruce/ash/maple)   tested on the other plot
  intensity only                 97 % (within)           91 % (within)                 92 % and 97 %
  crown shape only               95 % (within)           66 % (within)                 47 % and 37 %  (chance)
  all features                   95 % (within)           87 % (within)                 83 % and 70 %
Use intensity. Crown shape is species-specific: spruce crowns are narrow and pointed, pine crowns broad and flat,
so a shape model trained on pine/beech does not recognise spruce at all; it only helps when the training data
contain the same species. Aircraft scans (LGB): 94-97 % on two of four tiles, 67-73 % on the other two;
transfer between flight campaigns about 70 %. Pine and spruce echo 2-4x stronger than leafless broadleaf.
Leaf-off scans only.

Intensity is not comparable between flight strips or campaigns (raw ground intensity 181-1525 between LGB
flights). Normalise per strip to the ground echo (normalise_intensity) and calibrate each campaign on its own
labelled trees. Crown shape needs no calibration but well-segmented trees.

  F = tree_features(P, intensity_norm, ground_z)    # one tree
  model = train(features_df, labels)                 # labels: True = broadleaf; intensity features by default
  pred = model.predict(features_df[DEFAULT])
"""
import numpy as np, pandas as pd
from scipy.spatial import ConvexHull
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_predict, LeaveOneOut

SHAPE = ['crown_d', 'p25', 'p50', 'p90', 'sd']            # crown shape without height
INTENSITY = ['i_mean', 'i_top', 'i_p90', 'i_sd']
FEATURES = SHAPE + ['density'] + INTENSITY
DEFAULT = INTENSITY      # the only features that transferred between plots and species


def normalise_intensity(intensity, classification, return_number, point_source_id):
    """Scale each flight strip so its first-return ground echoes have median 1 (LAS fields)."""
    it = np.asarray(intensity, float); out = np.full(len(it), np.nan)
    g = (np.asarray(classification) == 2) & (np.asarray(return_number) == 1)
    ref = np.median(it[g])
    for k in np.unique(point_source_id):
        s = point_source_id == k; m = s & g
        out[s] = it[s] / (np.median(it[m]) if m.sum() > 1000 else ref)
    return out


def tree_features(P, intensity, ground):
    """Features of one segmented tree: P N x 3, intensity normalised as above (or raw within one flight)."""
    if len(P) < 10: return None
    h = P[:, 2].max() - ground; H = (P[:, 2] - ground) / h
    A = ConvexHull(P[:, :2]).volume; it = np.asarray(intensity, float); top = H > 0.8
    return dict(height=h, crown_d=2 * np.sqrt(A / np.pi), density=len(P) / A,
                p25=np.percentile(H, 25), p50=np.percentile(H, 50), p90=np.percentile(H, 90), sd=H.std(),
                i_mean=it.mean(), i_top=it[top].mean(), i_p90=np.percentile(it, 90), i_sd=it.std())


def model():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, class_weight='balanced'))


def train(F, y, cols=DEFAULT):
    return model().fit(F[cols], y)


def cross_validate(F, y, cols=DEFAULT):
    p = cross_val_predict(model(), F[cols], y, cv=LeaveOneOut()); return float(np.mean(p == np.asarray(y)))
