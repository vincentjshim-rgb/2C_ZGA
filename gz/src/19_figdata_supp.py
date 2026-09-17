"""Source data for Supplementary Figures S1, S2 and S7. No new analysis.
S1: GSE66582 library clock values (results/gate2b_R2_library_tage.tsv, copied).
S2: post-hoc zygotic timing intervals (results/gate1alt_intervals.tsv, copied).
S7: per-library quality control of the 76 Gate 3 libraries, recomputed with the rule of src/11_gate3_analysis.py
    (p_pseudoaligned >= 30 and log10 detected genes >= dataset median - 3 MAD); the passing set is checked against the
    stored Gate 3 library table before writing.
Output: figures/source_data/figS7_library_qc.tsv (+ copies)"""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import os, re
import numpy as np, pandas as pd

B = f'{ANALYSIS}/gz'
OUT = f'{B}/figures/source_data'
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
L3 = pd.read_csv(f'{B}/results/gate3_library_tage.tsv', sep='\t')
L3 = L3[L3.variant == 'V0'].set_index('run')
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
T2G = dict(zip(t2g.tx, t2g.gene_id))

rows = []
for gse, S in SEL.groupby('gse', sort=False):
    det, meta = {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{gse}/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        g = a.est_counts.groupby(a.index.to_series().map(T2G).values).sum()
        info = open(f'{q}/run_info.json').read()
        det[r.run_accession] = int((g >= 1).sum())
        meta.append(dict(gse=gse, run=r.run_accession, title=r.title_clean,
                         n_processed=int(re.search(r'"n_processed": (\d+)', info).group(1)),
                         p_pseudoaligned=float(re.search(r'"p_pseudoaligned": ([\d.]+)', info).group(1))))
    M = pd.DataFrame(meta).set_index('run')
    d = np.log10(pd.Series(det))
    mad = 1.4826 * np.median(np.abs(d - d.median()))
    M['detected_genes'] = pd.Series(det)
    M['detected_log10'] = d
    M['detected_cut_log10'] = d.median() - 3 * mad
    M['qc_pass'] = (M.p_pseudoaligned >= 30) & (d >= d.median() - 3 * mad)
    stored = set(L3.index[L3.gse == gse])
    assert set(M.index[M.qc_pass]) == stored, gse
    M['stage'] = [L3.loc[i, 'stage'] if i in stored else '' for i in M.index]
    M['arm'] = [L3.loc[i, 'arm'] if i in stored else '' for i in M.index]
    rows.append(M.reset_index())
Q = pd.concat(rows)
assert len(Q) == 76 and int((~Q.qc_pass).sum()) == 3, (len(Q), int((~Q.qc_pass).sum()))
Q.to_csv(f'{OUT}/figS7_library_qc.tsv', sep='\t', index=False)
print(Q[~Q.qc_pass][['gse', 'run', 'title', 'p_pseudoaligned', 'detected_genes']].to_string(index=False))
print(Q.groupby('gse').agg(n=('run', 'size'), pass_=('qc_pass', 'sum'), pal_min=('p_pseudoaligned', 'min'),
                           det_min=('detected_genes', 'min')).to_string())

for f in ['gate2b_R2_library_tage', 'gate2b_R2_intervals', 'gate1alt_intervals']:
    pd.read_csv(f'{B}/results/{f}.tsv', sep='\t').to_csv(f'{OUT}/{f}.tsv', sep='\t', index=False)
print('written to', OUT)
