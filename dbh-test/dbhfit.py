import numpy as np, pandas as pd
d = pd.read_csv(__import__('sys').argv[1] if len(__import__('sys').argv) > 1 else 'labels.csv')
d = d[(d.dTLS > 5) & (d.hTLS > 3) & (d.cpaHTLS > 0.5)].copy()
d['cd'] = 2 * np.sqrt(d.cpaHTLS / np.pi)
d['region'] = d.ep.str[:3].map({'SEW': 'Schorfheide (BB)', 'AEW': 'Alb', 'HEW': 'Hainich'})
def loo(sub, X):
    A = np.c_[np.ones(len(sub)), np.log(sub[X].values)]; Y = np.log(sub.dTLS.values); p = np.empty(len(sub))
    for i in range(len(sub)):
        m = np.arange(len(sub)) != i; b, *_ = np.linalg.lstsq(A[m], Y[m], rcond=None)
        p[i] = np.exp(A[i] @ b) * np.mean(np.exp(Y[m] - A[m] @ b))
    return p
viewer = lambda s: np.maximum(5, 1.5 * s.hTLS - 6)                       # current viewer rule (pine)
jucker = lambda s: 0.557 * (s.hTLS * s.cd) ** 0.809                      # Jucker et al. 2017 global (gymnosperms) from memory
def sc(p, t): e = p - t; return f'{np.mean(np.abs(e)):5.1f} cm ({np.mean(np.abs(e) / t) * 100:3.0f} %)'
for (reg, sp), sub in d.groupby(['region', 'species']):
    if len(sub) < 30: continue
    print(f'{reg:18s} {sp:18s} n={len(sub):4d} DBH {sub.dTLS.min():.0f}-{sub.dTLS.max():.0f} cm')
    print('   viewer rule 1.5h-6 :', sc(viewer(sub), sub.dTLS.values))
    print('   height only (fit)  :', sc(loo(sub, ['hTLS']), sub.dTLS.values))
    print('   height + crown     :', sc(loo(sub, ['hTLS', 'cd']), sub.dTLS.values))
pine = d[(d.region == 'Schorfheide (BB)') & (d.species == 'Pinus_sylvestris')]
A = np.c_[np.ones(len(pine)), np.log(pine[['hTLS', 'cd']].values)]; b, *_ = np.linalg.lstsq(A, np.log(pine.dTLS.values), rcond=None)
cf = np.mean(np.exp(np.log(pine.dTLS.values) - A @ b))
print('Schorfheide pine: DBH = %.4f * h^%.3f * cd^%.3f' % (np.exp(b[0]) * cf, b[1], b[2]), ' n=', len(pine))
A1 = np.c_[np.ones(len(pine)), np.log(pine[['hTLS']].values)]; b1, *_ = np.linalg.lstsq(A1, np.log(pine.dTLS.values), rcond=None)
cf1 = np.mean(np.exp(np.log(pine.dTLS.values) - A1 @ b1)); print('Schorfheide pine height-only: DBH = %.4f * h^%.3f' % (np.exp(b1[0]) * cf1, b1[1]))
print('ranges h', pine.hTLS.min(), pine.hTLS.max(), 'cd', pine.cd.min().round(1), pine.cd.max().round(1))
