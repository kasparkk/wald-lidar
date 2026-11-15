"""Tree species from DetailView features + echo + crown shape, adapted to one plot with a few known trees.

DetailView alone was poor on leaf-off broadleaf and on drone scans (BioDiv plots: 32-80 % balanced accuracy for the
plot's own species). A small logistic head on its features plus echo strength and crown shape, trained on a few
trees of the same plot whose species is known (e.g. from the field inventory), reached:
  known trees per species      3      10     20
  ground scan, pine/beech     93 %   94 %   97 %
  ground scan, spruce/broadl. 92 %   98 %   96 %
  drone, pine/beech           93 %   96 %   96 %
  drone, spruce/broadl.       85 %   88 %   90 %
(20 random draws each; test = all other trees of the plot.) It does NOT transfer to an unseen plot better than
DetailView alone: for conifer/broadleaf on a new plot use DetailView (ground scans) or echo strength (drone scans).

  import species_local as SL
  table = SL.join(dv_csv, feature_csv)            # outputs of dv_extract.py and species.tree_features per tree
  pred, prob = SL.fit_predict(table, known)       # known: dict file -> species group for the labelled trees
"""
import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

ECHO = ['i_mean', 'i_top', 'i_p90', 'i_sd']; SHAPE = ['crown_d', 'p25', 'p50', 'p90', 'sd', 'height']
Z = [f'z{i}' for i in range(256)]

def join(dv_csv, feature_csv):
    E = pd.read_csv(dv_csv).rename(columns={f'p{i}': f'q{i}' for i in range(33)})
    return pd.read_csv(feature_csv).merge(E, on='file')

def fit_predict(T, known, n_pca=10):
    """T: one plot, one scanner. known: {file: species group}. Returns predicted group and its probability per tree."""
    tr = T.file.isin(known).values; y = T.file[tr].map(known).values
    Q = [f'q{i}' for i in range(33)]
    pca = make_pipeline(StandardScaler(), PCA(min(n_pca, tr.sum() - 1), random_state=0)).fit(T.loc[tr, Z])
    A = np.c_[pca.transform(T[Z]), T[Q].values, T[ECHO + SHAPE].values]
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, class_weight='balanced', C=0.5)).fit(A[tr], y)
    P = m.predict_proba(A)
    return pd.Series(m.classes_[P.argmax(1)], index=T.file), pd.Series(P.max(1), index=T.file)
