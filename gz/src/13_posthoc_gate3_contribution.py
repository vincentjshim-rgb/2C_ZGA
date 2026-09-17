"""POST HOC (plan: plan/POSTHOC_gate3_contribution_frozen.md). Decompose the within-2-cell tAge drop of
P1 GSE280522 (and P3 GSE300734 as a direction check) into per-gene contributions, and split the contributions
A485 removes. Descriptive only; no p-value, no gate. Loading and QC mirror src/11_gate3_analysis.py."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import re, json, os
import numpy as np, pandas as pd, joblib
import sys
B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import tage_py as tp

MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
t2g = dict(zip(t2g.tx, t2g.gene_id))
GT = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
GT = GT[GT.Entrez.notna()].drop_duplicates('Ensembl')
SYM = dict(zip(GT.Entrez.astype('int64').astype(str), GT['Gene.Symbol']))


def annotate(gse, t):
    if gse == 'GSE280522':
        m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
        return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]
    m = re.search(r'(WT|KO)_(E2C|L2C)_RNA', t)                      # GSE300734
    return m.group(2), {'WT': 'control', 'KO': 'Brg1_matKO'}[m.group(1)]


def load_dataset(gse):
    S = SEL[SEL.gse == gse]
    cols, meta = {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{gse}/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', open(f'{q}/run_info.json').read()).group(1))
        cols[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(t2g).values).sum()
        stage, arm = annotate(gse, r.title_clean)
        meta.append(dict(run=r.run_accession, stage=stage, arm=arm, p_pseudoaligned=p_al))
    C = pd.DataFrame(cols).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M['qc_pass'] = (M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)
    M = M[M.qc_pass]
    return C[M.index], M


def entrez_matrix(C):
    """Counts collapsed to Entrez exactly as tage_py.preprocess does, for CPM-based labels."""
    e = tp.map_mouse_ensembl(tp.filter_genes(C))
    return e / e.sum(axis=0) * 1e6


def labels(cpm, M):
    e = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
    l = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
    mat = ((cpm[e] >= 10).mean(axis=1) >= 0.8) & (cpm[l].median(axis=1) <= cpm[e].median(axis=1) / 8)
    zyg = ((cpm[e] < 1).mean(axis=1) >= 0.8) & ((cpm[l] >= 10).mean(axis=1) >= 0.5)
    out = pd.Series('neither', index=cpm.index)
    out[mat.values] = 'maternal'
    out[zyg.values] = 'zygotic'
    return out


def contributions(gse, perturbed_arm):
    C, M = load_dataset(gse)
    pp = tp.preprocess(C, reference_samples=list(M.index[(M.arm == 'control') & (M.stage == 'E2C')]))
    X = pp['scaled_diff']                                   # features (Entrez) x samples
    model, feats = tp.load_clock(MODEL)
    sel = model.named_steps['featureselection'].get_support()
    coef = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
    coef = coef[coef != 0]
    out = {}
    for arm in ['control', perturbed_arm]:
        e = M.index[(M.arm == arm) & (M.stage == 'E2C')]
        l = M.index[(M.arm == arm) & (M.stage == 'L2C')]
        d = X.reindex(coef.index)[l].mean(axis=1) - X.reindex(coef.index)[e].mean(axis=1)
        out[arm] = (coef * d.fillna(0))                     # absent genes: imputed constant -> exactly 0
    cpm = entrez_matrix(C)
    lab = labels(cpm, M)
    ce = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
    cl = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
    lfc = np.log2((cpm[cl].mean(axis=1) + 1) / (cpm[ce].mean(axis=1) + 1))
    T = pd.DataFrame({'coef': coef, 'c_control': out['control'], 'c_perturbed': out[perturbed_arm]})
    T['delta'] = T.c_perturbed - T.c_control
    T['symbol'] = [SYM.get(i, '') for i in T.index]
    T['label'] = lab.reindex(T.index).fillna('not_expressed')
    T['log2FC_control_L2C_vs_E2C'] = lfc.reindex(T.index)
    return T.sort_values('c_control')


def summarise(T, tag, D_expected):
    D = T.c_control.sum()
    assert abs(D - D_expected) < 1e-9, (tag, D, D_expected)
    rows = []
    neg = T[T.c_control < 0].c_control.sort_values()         # genes pushing the score down
    for k in [1, 5, 10, 25, 50, 100]:
        rows.append(dict(dataset=tag, quantity='share_of_drop_top_genes', key=f'top{k}',
                         value=neg.head(k).sum() / D))
    cum = neg.cumsum() / D
    for q in [0.5, 0.8]:
        rows.append(dict(dataset=tag, quantity='genes_reaching_share_of_drop', key=f'{int(q*100)}%',
                         value=float((cum < q).sum() + 1)))
    up = T.log2FC_control_L2C_vs_E2C > 0
    for nm, mask in [('rises_x_coef_neg', up & (T.coef < 0)), ('rises_x_coef_pos', up & (T.coef > 0)),
                     ('falls_x_coef_neg', ~up & (T.coef < 0)), ('falls_x_coef_pos', ~up & (T.coef > 0))]:
        rows.append(dict(dataset=tag, quantity='direction_split_share_of_D', key=nm,
                         value=T.c_control[mask].sum() / D))
        rows.append(dict(dataset=tag, quantity='direction_split_share_of_I', key=nm,
                         value=T.delta[mask].sum() / T.delta.sum()))
    for nm in ['maternal', 'zygotic', 'neither', 'not_expressed']:
        m = T.label == nm
        rows.append(dict(dataset=tag, quantity='programme_share_of_D', key=nm, value=T.c_control[m].sum() / D))
        rows.append(dict(dataset=tag, quantity='programme_share_of_I', key=nm, value=T.delta[m].sum() / T.delta.sum()))
        rows.append(dict(dataset=tag, quantity='programme_n_genes', key=nm, value=float(m.sum())))
    rows.append(dict(dataset=tag, quantity='total', key='D_control', value=D))
    rows.append(dict(dataset=tag, quantity='total', key='D_perturbed', value=T.c_perturbed.sum()))
    rows.append(dict(dataset=tag, quantity='total', key='I', value=T.delta.sum()))
    rows.append(dict(dataset=tag, quantity='total', key='n_nonzero_clock_genes_present',
                     value=float((T.c_control != 0).sum())))
    return pd.DataFrame(rows)


A = pd.read_csv(f'{B}/results/gate3_arm_drops.tsv', sep='\t')


def expected(gse, arm):
    r = A[(A.gse == gse) & (A.variant == 'V0') & (A.arm == arm)]
    return float(r['drop'].iloc[0])


T1 = contributions('GSE280522', 'A485')
T3 = contributions('GSE300734', 'Brg1_matKO')
S = pd.concat([summarise(T1, 'P1_GSE280522_A485', expected('GSE280522', 'control')),
               summarise(T3, 'P3_GSE300734_Brg1KO', expected('GSE300734', 'control'))])

T1.to_csv(f'{B}/results/posthoc_gate3_contribution_genes.tsv', sep='\t')
T3.to_csv(f'{B}/results/posthoc_gate3_contribution_genes_P3.tsv', sep='\t')
S.to_csv(f'{B}/results/posthoc_gate3_contribution_summary.tsv', sep='\t', index=False)

pd.set_option('display.width', 250)
print('=== P1 GSE280522: 25 genes contributing most to the control drop (POST HOC, descriptive) ===')
print(T1.head(25)[['symbol', 'coef', 'c_control', 'c_perturbed', 'delta', 'label',
                   'log2FC_control_L2C_vs_E2C']].round(4).to_string())
print('\n=== summary ===')
print(S.round(4).to_string(index=False))
