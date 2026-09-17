"""Source data for the publication figures (plan: results/FIGURE_PLAN_gz_v1_KR.md). No new analysis: values already
defined in Gates 3-4 and post-hoc analyses are re-extracted per library for display, and each re-extraction is checked
against the stored result before it is written. PCA of PSI is a display (no test).
Output: figures/source_data/*.tsv"""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import os, re, sys
import numpy as np, pandas as pd

B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import tage_py as tp

OUT = f'{B}/figures/source_data'
os.makedirs(OUT, exist_ok=True)
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
L3 = pd.read_csv(f'{B}/results/gate3_library_tage.tsv', sep='\t')
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
T2G = dict(zip(t2g.tx, t2g.gene_id))
GT = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
GT = GT[GT.Entrez.notna()].drop_duplicates('Ensembl')
ENS2SYM = dict(zip(GT.Ensembl, GT['Gene.Symbol']))
SYM2ENS = {s: e for e, s in ENS2SYM.items()}
GTF_NAME = {}                                   # display fallback for genes without an Entrez symbol
with open(f'{B}/data/annotation/Mus_musculus.GRCm39.112.gtf') as fh:
    for ln in fh:
        if '\tgene\t' in ln:
            m1, m2 = re.search(r'gene_id "([^"]+)"', ln), re.search(r'gene_name "([^"]+)"', ln)
            if m1 and m2:
                GTF_NAME[m1.group(1)] = m2.group(1)
display = lambda g: ENS2SYM.get(g) or GTF_NAME.get(g, g)


def p1_title(t):
    m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep(\d)', t)
    return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)], int(m.group(3))


# ---- P1 libraries (Gate 3 QC set) and gene-level CPM --------------------------------------------------
qc = set(L3[(L3.gse == 'GSE280522') & (L3.variant == 'V0')].run)
S1 = SEL[(SEL.gse == 'GSE280522') & SEL.run_accession.isin(qc)]
cols, meta = {}, []
for _, r in S1.iterrows():
    a = pd.read_csv(f'{B}/data/perturb/GSE280522/quant/{r.run_accession}/abundance.tsv', sep='\t', index_col=0)
    cols[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(T2G).values).sum()
    st, arm, rep = p1_title(r.title_clean)
    meta.append(dict(run=r.run_accession, stage=st, arm=arm, rep=rep))
C = pd.DataFrame(cols).fillna(0)
M = pd.DataFrame(meta).set_index('run').sort_values(['arm', 'stage', 'rep'])
C = C[M.index]
cpm = C / C.sum(axis=0) * 1e6
lg = np.log2(cpm + 1)
assert len(M) == 23, len(M)

# zygotic set: control-arm V2 rule, as src/14_posthoc_rescue.py
e = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
l = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
zyg = ((cpm[e] < 1).mean(axis=1) >= 0.8) & ((cpm[l] >= 10).mean(axis=1) >= 0.5)
Z = lg.loc[zyg.values]
assert len(Z) == 13, len(Z)
score = Z.mean(axis=0)
ref = pd.read_csv(f'{B}/results/posthoc_rescue_programme_scores.tsv', sep='\t')
for arm, key in [('control', 'control L2C'), ('A485', 'A485 L2C'), ('A485+Dux', 'A485+Dux L2C')]:
    got = score[M.index[(M.arm == arm) & (M.stage == 'L2C')]].mean()
    want = ref[(ref.set == 'zygotic') & (ref.arm == key)]['mean'].iloc[0]
    assert abs(got - want) < 1e-9, (arm, got, want)
lib = M.copy()
lib['zygotic_score'] = score
lib.reset_index().to_csv(f'{OUT}/fig4b_zygotic_score_per_library.tsv', sep='\t', index=False)

# heatmap genes: zygotic set + 2-cell/DUX-target markers (fixed regex of src/14), z-scored per gene across libraries
MARK = re.compile(r'^(Zscan4|Tcstv|Tdpoz|Dux)', re.I)
mk = [g for g in lg.index if MARK.match(str(ENS2SYM.get(g, ''))) and lg.loc[g].max() > 0]
rows = []
for g in list(Z.index) + [g for g in mk if g not in Z.index]:
    v = lg.loc[g]
    z = (v - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else v * 0
    rows.append(dict(gene_id=g, symbol=display(g), group='zygotic set' if g in Z.index else '2C/DUX marker',
                     **{f'{run}': z[run] for run in M.index}))
H = pd.DataFrame(rows)
H.to_csv(f'{OUT}/fig4a_heatmap_z.tsv', sep='\t', index=False)
M.reset_index().to_csv(f'{OUT}/fig4a_heatmap_columns.tsv', sep='\t', index=False)
print('heatmap genes:', len(H), '(zygotic', int((H.group == 'zygotic set').sum()), ')')

# key contributor genes, log2(CPM+1) per library
KEY = ['Klf9', 'Neto2', 'Pi4k2a', 'Gpatch4', 'Psmb5']
kr = []
for s in KEY:
    g = SYM2ENS[s]
    for run in M.index:
        kr.append(dict(symbol=s, gene_id=g, run=run, arm=M.loc[run, 'arm'], stage=M.loc[run, 'stage'], log2cpm=lg.loc[g, run]))
pd.DataFrame(kr).to_csv(f'{OUT}/fig3d_key_genes_per_library.tsv', sep='\t', index=False)
T = pd.read_csv(f'{B}/results/posthoc_gate3_contribution_genes.tsv', sep='\t', index_col=0)
for s in KEY:
    ent = GT.loc[GT['Gene.Symbol'] == s, 'Entrez'].astype('int64').astype(str).iloc[0]
    lfc_file = T.loc[int(ent), 'log2FC_control_L2C_vs_E2C'] if int(ent) in T.index else T.loc[ent, 'log2FC_control_L2C_vs_E2C']
    g = SYM2ENS[s]
    lfc_here = np.log2((cpm.loc[g, l].mean() + 1) / (cpm.loc[g, e].mean() + 1))
    print(f'  {s}: log2FC (Entrez, decomposition) {lfc_file:.3f} vs Ensembl gene here {lfc_here:.3f}')

# ---- library-level clock values for Figure 2b (Gate 3, V0) -------------------------------------------
L3[L3.variant == 'V0'].to_csv(f'{OUT}/fig2b_clock_per_library.tsv', sep='\t', index=False)

# ---- PSI PCA (display) for the three primary datasets --------------------------------------------------
def read_expr(path):
    with open(path) as fh:
        c = fh.readline().rstrip('\n').split('\t')
    df = pd.read_csv(path, sep='\t', skiprows=1, header=None, index_col=0)
    df.columns = c
    return df


G4 = pd.read_csv(f'{B}/results/gate4_progression.tsv', sep='\t')
pcs = []
for gse in ['GSE280522', 'GSE221985', 'GSE300734']:
    runs = list(G4[G4.gse == gse].run)
    ev = pd.read_csv(f'{B}/results/gate4_events_filtered_{gse}.tsv', sep='\t', index_col=0)
    P = read_expr(f'{B}/data/gate4/{gse}.psi').loc[ev.index, runs]
    P = P.apply(lambda r: r.fillna(r.mean()), axis=1)
    X = P.sub(P.mean(axis=1), axis=0).T.values
    U, Sv, Vt = np.linalg.svd(X, full_matrices=False)
    var = Sv ** 2 / (Sv ** 2).sum()
    sc = U[:, :2] * Sv[:2]
    # orient PC1 so that control late two-cell is positive
    g = G4[G4.gse == gse].set_index('run')
    ctl_l = [i for i, r in enumerate(runs) if g.loc[r, 'arm'] == 'control' and g.loc[r, 'stage'] == 'L2C']
    if sc[ctl_l, 0].mean() < 0:
        sc[:, 0] *= -1
    for i, r in enumerate(runs):
        pcs.append(dict(gse=gse, run=r, arm=g.loc[r, 'arm'], stage=g.loc[r, 'stage'], PC1=sc[i, 0], PC2=sc[i, 1],
                        var_PC1=var[0], var_PC2=var[1], n_events=len(P)))
PCA = pd.DataFrame(pcs)
PCA.to_csv(f'{OUT}/fig5a_psi_pca.tsv', sep='\t', index=False)
print(PCA.groupby('gse')[['var_PC1', 'var_PC2', 'n_events']].first().round(3).to_string())

# ---- copies of stored tables used as-is ----------------------------------------------------------------
for f in ['gate3_arm_drops', 'gate3_secondary_single_stage', 'gate4_progression', 'gate4_arm_progression',
          'posthoc_gate3_contribution_genes', 'posthoc_gate3_contribution_genes_P3', 'posthoc_rescue_contributions',
          'posthoc_gate4_crossfit_folds', 'gate1_embryo_tage', 'gate1_intervals', 'gate2b_R1_embryo_tage',
          'gate2a_intervals', 'gate2a_simulation']:
    pd.read_csv(f'{B}/results/{f}.tsv', sep='\t').to_csv(f'{OUT}/{f}.tsv', sep='\t', index=False)
print('source data written to', OUT)
