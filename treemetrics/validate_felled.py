"""Test against real wood: runs run_trees.measure on the 65 felled and weighed Belgian trees
(Demol, Gielen & Verbeeck 2021, Zenodo 4557401, CC BY 4.0) and refits the corrections leave-one-out,
so every tree is predicted with constants estimated without it.

  python3 validate_felled.py DEMOL_FOLDER [--bin PATH/raycloudtools/build/raycloudtools]

DEMOL_FOLDER must contain Destructive_and_qsm_data_DEMOL.csv and pointclouds_clean/ (unpacked from
pointclouds_clean.7z). Takes about 15-20 minutes. Writes felled_results.csv. Use the same script with your
own felled trees to refit corrections.TWIG_COEF, BROADLEAF_FACTOR and DBH_OFFSET_CM.
"""
import argparse, os, sys
import numpy as np, pandas as pd
import treevol as TV, corrections as C, run_trees as RT

SPECIES = {'FSYL': 'Fagus sylvatica', 'PSYL': 'Pinus sylvestris', 'FEXC': 'Fraxinus excelsior', 'LXDC': 'Larix decidua'}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder'); ap.add_argument('--bin', default='')
    a = ap.parse_args()
    d = pd.read_csv(os.path.join(a.folder, 'Destructive_and_qsm_data_DEMOL.csv'))
    d['key'] = d.tree_name.str.replace('-', '').str.replace(r'([A-Z]+)0*(\d+)', r'\1\2', regex=True)
    rows = []
    for _, t in d.iterrows():
        f = os.path.join(a.folder, 'pointclouds_clean', t.key + '.txt')
        if not os.path.exists(f): continue
        P = np.loadtxt(f, usecols=(0, 1, 2)); sp = SPECIES.get(t.key[:4], '')
        r, S = RT.measure(P, sp, a.bin)
        if S is None: continue
        rows.append(dict(tree=t.key, species=sp, ref_l=t.Volume_total_tree_harvested, tape_dbh_cm=t.DBH,
                         s3=C.tip_r3(S), conifer=C.is_conifer(sp), **r))
        print(rows[-1]['tree'], r, flush=True)
    m = pd.DataFrame(rows).replace('', np.nan)
    for c in ['volume_l', 'volume_rayextract_l', 'dbh_cm', 'height_m']: m[c] = m[c].astype(float)
    n = len(m); idx = np.arange(n); loo = np.empty(n)
    for i in range(n):            # refit both corrections without tree i
        tr = idx != i; con = m.conifer.values
        s = tr & con; X = m.s3.values[s] * 1000; y = (m.ref_l - m.volume_l).values[s]; coef = X @ y / (X @ X)
        s = tr & ~con; fac = np.exp(np.median(np.log(m.ref_l.values[s] / m.volume_l.values[s])))
        loo[i] = m.volume_l.iloc[i] + coef * m.s3.iloc[i] * 1000 if con[i] else m.volume_l.iloc[i] * fac
    m['volume_loo_l'] = loo; m.to_csv('felled_results.csv', index=False)
    print(f'\nWood volume against weighed wood ({n} trees)')
    for c, lbl in [('volume_rayextract_l', 'RayExtract'), ('volume_l', 'RayExtract + measured stem'), ('volume_loo_l', '+ corrections (leave-one-out)')]:
        q = m[c] / m.ref_l - 1
        print(f'  {lbl:32s} within 25 %: {(q.abs() <= .25).sum()}/{n}  within 10 %: {(q.abs() <= .1).sum()}/{n}  mean {q.abs().mean()*100:4.1f} %  total {(m[c].sum()/m.ref_l.sum()-1)*100:+.0f} %')
    e = (m.dbh_cm - m.tape_dbh_cm).dropna()
    print(f'DBH vs tape: {len(e)} plausible trees, mean deviation {e.abs().mean():.2f} cm, bias {e.mean():+.2f} cm')
    for s, g in m.groupby('species'):
        eg = (g.dbh_cm - g.tape_dbh_cm).dropna(); print(f'   {s:20s} median offset {eg.median():+.2f} cm')


if __name__ == '__main__':
    main()
