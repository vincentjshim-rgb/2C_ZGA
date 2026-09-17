"""Validate src/tage_py.py against tAge golden values (TACO head-to-head, tests/testthat/test-predict.R)."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys, numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp
B = f'{ANALYSIS}/gz'
ex = pd.read_csv(f'{tp.PKG}/extdata/Exprs_example.csv', index_col=0)
pp = tp.preprocess(ex)
print('example: %d genes x %d samples; after filter %d; entrez %d' % (*ex.shape, pp['n_genes_after_filter'], pp['n_entrez']))
gold = {'EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl': ([-0.593993, -0.668232, -0.530628], 1e-3),
        'EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl': ([-9.7581, -7.7424, -6.3348], 1e-2)}
ok = True
rows = []
for f, (g, tol) in gold.items():
    v = tp.predict(f'{B}/tools_tAge_models/{f}', pp['scaled_diff'], 'mouse').values[:3]
    d = np.max(np.abs(v - np.array(g)))
    passed = d <= tol
    ok &= passed
    print('%-55s got %s gold %s maxdiff %.2g %s' % (f, np.round(v, 6), g, d, 'PASS' if passed else 'FAIL'))
    rows.append(dict(model=f, got=';'.join('%.6f' % x for x in v), gold=';'.join(map(str, g)), max_abs_diff=d, tol=tol, passed=passed))
pd.DataFrame(rows).to_csv(f'{B}/results/00_tage_py_validation.tsv', sep='\t', index=False)
print('OVERALL', 'PASS' if ok else 'FAIL')
