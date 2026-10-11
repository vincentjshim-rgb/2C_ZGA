"""EXPLORATION A (plan: plan/EXPLORE_A_netgross_frozen.md, frozen 2026-10-12 before any gross value existed).

For every comparison already scored in this study, read the clock twice from the same per-gene terms
t_i = c_i (xbar_late,i - xbar_early,i):   net D = sum t_i (the reported clock change),  gross G = sum |t_i|,
down / up (negative / positive terms) and the cancellation C = |D| / G.

Loading and preprocessing copy the gate scripts line for line (01, 06, 08, 11/13, 31) and every recomputed net D is
asserted against the stored result before anything new is written, so the gross numbers rest on exactly the matrices
the paper already reports. Bootstrap: libraries resampled within stage (or arm x stage), 2,000 draws, seed 20261012.

Run with the gz environment python.
"""
import glob
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import tage_py as tp  # noqa: E402

MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
NBOOT, SEED = 2000, 20261012
rng = np.random.default_rng(SEED)

model, feats = tp.load_clock(MODEL)
sel = model.named_steps['featureselection'].get_support()
COEF = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
COEF = COEF[COEF != 0]
assert len(COEF) == 1839, len(COEF)

# literature major-ZGA intervals fixed in the Gate 1 plan; human from Braude et al. (1988)
ZGA = {'Mouse': ('Early-2-cell', 'Late-2-cell'), 'Pig': ('2-cell', '4-cell'), 'Cow': ('8-cell', '16-cell'),
       'Rabbit': ('8-cell', '16-cell'), 'Human': ('4-cell', '8-cell')}


# ====================================================================== the two readings of one comparison
def terms_matrix(X):
    """Clock features x samples, absent genes as 0 (they are imputed to a constant and contribute nothing)."""
    return X.reindex(COEF.index).fillna(0.0)


def budget(T, early_cols, late_cols):
    t = COEF.values * (T[late_cols].mean(axis=1).values - T[early_cols].mean(axis=1).values)
    return dict(D=float(t.sum()), G=float(np.abs(t).sum()), down=float(t[t < 0].sum()), up=float(t[t > 0].sum()))


def budget_boot(T, groups, pairs, nboot=NBOOT):
    """Bootstrap the stage means (libraries resampled within group) and return D and G per pair per draw."""
    M = T.values
    idx = {g: [T.columns.get_loc(c) for c in cols] for g, cols in groups.items()}
    D = np.empty((nboot, len(pairs)))
    G = np.empty((nboot, len(pairs)))
    for b in range(nboot):
        means = {g: M[:, rng.choice(ix, len(ix))].mean(axis=1) for g, ix in idx.items()}
        for k, (a, c) in enumerate(pairs):
            t = COEF.values * (means[c] - means[a])
            D[b, k] = t.sum()
            G[b, k] = np.abs(t).sum()
    return D, G


def interval_table(T, stage, order, dataset, species, zga=None):
    stages = [s for s in order if (stage == s).any()]
    groups = {s: list(stage.index[stage == s]) for s in stages}
    pairs = list(zip(stages[:-1], stages[1:]))
    Db, Gb = budget_boot(T, groups, pairs)
    amax = np.argmax(Gb, axis=1)
    rows = []
    for k, (a, c) in enumerate(pairs):
        bgt = budget(T, groups[a], groups[c])
        rows.append(dict(dataset=dataset, species=species, interval=f'{a}->{c}', n_early=len(groups[a]),
                         n_late=len(groups[c]), **bgt, C=abs(bgt['D']) / bgt['G'],
                         D_lo=np.percentile(Db[:, k], 2.5), D_hi=np.percentile(Db[:, k], 97.5),
                         G_lo=np.percentile(Gb[:, k], 2.5), G_hi=np.percentile(Gb[:, k], 97.5),
                         P_max_G=float(np.mean(amax == k)), is_zga=bool(zga is not None and (a, c) == zga)))
    return pd.DataFrame(rows)


# ====================================================================== loaders, copied from the gate scripts
def load_gse225056(sp):                                   # src/01_gate1_cross_species.py
    ORDER = {'Mouse': ['Oocyte', 'Zygote', 'Early-2-cell', 'Late-2-cell', '4-cell', '8-cell', '16-cell'],
             'Cow': ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
             'Rabbit': ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
             'Pig': ['Oocyte', 'Zygote', '2-cell', '4-cell', 'Day2', 'Day3', 'Morula']}
    ORTHO = {'Cow': 'btaurus', 'Pig': 'sscrofa', 'Rabbit': 'ocuniculus'}
    C = pd.read_csv(f'{B}/data/GSE225056_{sp}_count_5prime.txt.gz', sep='\t', index_col=0)
    C.index = C.index.str.strip('"')
    C.columns = C.columns.str.strip('"')
    if sp == 'Cow':
        C = C[C.index.str.match(r'^(NM|XM)_')]
        acc = C.index.str.replace(r'\.\d+$', '', regex=True)
        nm = pd.read_csv(f'{B}/metadata/cow_refseq_mrna_to_ensembl.tsv', sep='\t').set_axis(['gid', 'acc'], axis=1).dropna()
        nm = nm.drop_duplicates('acc')
        xm = pd.read_csv(f'{B}/metadata/cow_xm_symbols.tsv', sep='\t')
        xm = xm[xm.symbol.notna() & ~xm.symbol.astype(str).str.startswith('LOC')]
        nameO = pd.read_csv(f'{B}/metadata/orthologs_btaurus_to_mouse.tsv', sep='\t').set_axis(['gid', 'name', 'mgid', 'type'], axis=1)
        nm2g = nameO.dropna(subset=['name']).drop_duplicates(['gid', 'name'])[['gid', 'name']]
        uniq = nm2g.groupby('name').gid.nunique()
        name2gid = nm2g[nm2g.name.isin(uniq.index[uniq == 1])].drop_duplicates('name').set_index('name').gid
        mp = dict(zip(nm.acc, nm.gid))
        mp.update({a.split('.')[0]: name2gid[s] for a, s in zip(xm.accession, xm.symbol) if s in name2gid.index})
        gids = pd.Series(acc.map(lambda a: mp.get(a)), index=C.index)
        ok = gids.notna()
        C = C.loc[ok].groupby(gids[ok].values).sum()
    else:
        C = C[C.index.str.startswith('ENS')]
        C.index = C.index.str.replace(r'\.\d+$', '', regex=True)
        C = C.groupby(level=0).sum()
    if sp in ORTHO:
        O = pd.read_csv(f'{B}/metadata/orthologs_{ORTHO[sp]}_to_mouse.tsv', sep='\t')
        O.columns = ['gid', 'name', 'mouse_gid', 'type']
        O = O[O.type == 'ortholog_one2one'].drop_duplicates('gid')
        mp = dict(zip(O.gid, O.mouse_gid))
        ids = C.index.to_series().map(mp)
        ok = ids.notna()
        C = C.loc[ok].groupby(ids[ok].values).sum()
    stage = pd.Series([re.sub(r'_\d+$', '', c) for c in C.columns], index=C.columns)
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    keep = det >= det.median() - 3 * mad
    C, stage = C.loc[:, keep], stage[keep]
    pp = tp.preprocess(C, reference_samples=list(stage.index[stage == 'Oocyte']))
    return pp['scaled_diff'], stage, ORDER[sp]


def load_gse45719():                                      # src/06_gate2b_R1_deng.py
    PREFIX = [('zy', 'Zygote'), ('early2cell', 'Early-2C'), ('mid2cell', 'Mid-2C'), ('late2cell', 'Late-2C'),
              ('4cell', '4C'), ('8cell', '8C'), ('16cell', '16C'), ('earlyblast', 'EarlyBl'), ('midblast', 'MidBl'),
              ('lateblast', 'LateBl')]
    ORDER = [s for _, s in PREFIX]
    EXCL = re.compile(r'split|pooled|smartseq2|C57twocell|fibroblast|liver|BXC', re.I)
    cells = {}
    for f in sorted(glob.glob(f'{B}/data/GSE45719/*_expression.txt.gz')):
        name = re.sub(r'^GSM\d+_', '', os.path.basename(f)).replace('_expression.txt.gz', '')
        if EXCL.search(name):
            continue
        st = next((s for p, s in PREFIX if name.startswith(p)), None)
        if st is None:
            continue
        if st == '8C' and not re.match(r'^8cell_(1|2|5|8)-', name):
            continue
        if st == '16C' and not re.match(r'^16cell_(1|4|5|6)-', name):
            continue
        embryo = name if st == 'Zygote' else name.rsplit('-', 1)[0]
        x = pd.read_csv(f, sep='\t', usecols=[0, 3])
        x.columns = ['symbol', 'reads']
        cells[name] = (st, embryo, x.groupby('symbol').reads.sum())
    meta = pd.DataFrame([(n, s, e) for n, (s, e, _) in cells.items()], columns=['cell', 'stage', 'embryo'])
    X = pd.DataFrame({n: v for n, (_, _, v) in cells.items()}).fillna(0)
    gm = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
    u = gm.drop_duplicates('Gene.Symbol', keep=False)
    sym2ens = dict(zip(u['Gene.Symbol'], u.Ensembl))
    ens = X.index.to_series().map(sym2ens)
    X = X.loc[ens.notna()].groupby(ens[ens.notna()].values).sum()
    P = X.T.groupby(meta.set_index('cell').embryo.reindex(X.columns).values).sum().T
    emb = meta.groupby('embryo').agg(stage=('stage', 'first'), n_cells=('cell', 'size')).loc[P.columns]
    det = np.log10((P >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    keep = det >= det.median() - 3 * mad
    P, emb = P.loc[:, keep], emb[keep]
    pp = tp.preprocess(P, reference_samples=list(emb.index[emb.stage == 'Zygote']))
    return pp['scaled_diff'], emb.stage, ORDER


def load_gse66582():                                      # src/08_gate2b_R2_wu.py
    RUNS = {'SRR2927026': 'MII', 'SRR2927027': 'MII', 'SRR1840514': 'Zygote', 'SRR1840515': 'Zygote',
            'SRR2927028': 'Early-2C', 'SRR2927029': 'Early-2C', 'SRR1840516': '2C', 'SRR1840517': '2C',
            'SRR1840518': '4C', 'SRR1840519': '4C', 'SRR1840520': '8C', 'SRR1840521': '8C',
            'SRR1840522': 'ICM', 'SRR1840523': 'ICM', 'SRR1840526': 'ICM'}
    ORDER = ['MII', 'Zygote', 'Early-2C', '2C', '4C', '8C']
    t2g = pd.read_csv(f'{ANALYSIS}/metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
    t2g = dict(zip(t2g.tx, t2g.gene_id))
    cols = {}
    for r in RUNS:
        a = pd.read_csv(f'{B}/data/GSE66582/quant/{r}/abundance.tsv', sep='\t', index_col=0)
        cols[r] = a.est_counts.groupby(a.index.to_series().map(t2g).values).sum()
    C = pd.DataFrame(cols).fillna(0)
    stage = pd.Series(RUNS)[C.columns]
    pp = tp.preprocess(C, reference_samples=list(stage.index[stage == 'MII']))
    return pp['scaled_diff'], stage, ORDER


SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
T2G = pd.read_csv(f'{ANALYSIS}/metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
T2G = dict(zip(T2G.tx, T2G.gene_id))


def annotate(gse, t):                                     # src/11_gate3_analysis.py
    if gse == 'GSE280522':
        m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
        return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]
    if gse == 'GSE221985':
        m = re.search(r'RNA_(E2C|L2C)_(Ctrl|KO)_rep', t)
        return m.group(1), {'Ctrl': 'control', 'KO': 'Tardbp_matKO'}[m.group(2)]
    m = re.search(r'(WT|KO)_(E2C|L2C)_RNA', t)
    return m.group(2), {'WT': 'control', 'KO': 'Brg1_matKO'}[m.group(1)]


def load_perturbation(gse):                               # src/13_posthoc_gate3_contribution.py
    S = SEL[SEL.gse == gse]
    cols, meta = {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{gse}/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', open(f'{q}/run_info.json').read()).group(1))
        cols[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(T2G).values).sum()
        st, arm = annotate(gse, r.title_clean)
        meta.append(dict(run=r.run_accession, stage=st, arm=arm, p_pseudoaligned=p_al))
    C = pd.DataFrame(cols).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M = M[(M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)]
    C = C[M.index]
    pp = tp.preprocess(C, reference_samples=list(M.index[(M.arm == 'control') & (M.stage == 'E2C')]))
    return pp['scaled_diff'], M


def load_gse36552():                                      # src/31_human_decomposition.py (inputs now in data/)
    D = f'{B}/data/GSE36552'
    ORDER = ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'Late blastocyst']
    S = pd.read_csv(f'{D}/samples.tsv', sep='\t').query("stage != 'ESC'")
    S['embryo'] = [f'{s}_{m.group(1)}' if (m := re.search(r'#\s*(\d+)', t)) else f'{s}_1' for t, s in zip(S.title, S.stage)]
    cols = {}
    for gsm, emb in zip(S.gsm, S.embryo):
        p = f'{D}/expr/{gsm}.txt.gz'
        if not os.path.exists(p):
            continue
        with gzip.open(p, 'rt') as fh:
            d = pd.read_csv(fh, sep='\t', usecols=[0, 1], names=['sym', 'n'], skiprows=1)
        cols.setdefault(emb, []).append(d.groupby('sym').n.sum())
    P = pd.DataFrame({e: pd.concat(v, axis=1).fillna(0).sum(axis=1) for e, v in cols.items()}).fillna(0)
    meta = S.drop_duplicates('embryo').set_index('embryo').stage.reindex(P.columns)
    H = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_human.csv').dropna(subset=['Entrez', 'Gene.Symbol'])
    sym2h = dict(zip(H.drop_duplicates('Gene.Symbol')['Gene.Symbol'].astype(str),
                     H.drop_duplicates('Gene.Symbol')['Entrez'].astype('int64')))
    O = pd.read_csv(f'{tp.PKG}/extdata/metadata/Table_of_orthologs.csv').dropna(subset=['Entrez.Human', 'Entrez.Mouse'])
    h2m = dict(zip(O['Entrez.Human'].astype('int64'), O['Entrez.Mouse'].astype('int64')))
    Mt = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv').dropna(subset=['Entrez', 'Ensembl'])
    m2e = dict(zip(Mt.drop_duplicates('Entrez')['Entrez'].astype('int64'), Mt.drop_duplicates('Entrez')['Ensembl'].astype(str)))
    step = P.index.to_series().map(sym2h).map(h2m).map(m2e)
    ok = step.notna()
    C = P.loc[ok].groupby(step[ok].values).sum()
    pp = tp.preprocess(C, reference_samples=list(meta.index[meta == 'Oocyte']))
    return pp['scaled_diff'], meta, ORDER


# ====================================================================== run, checking every net against the record
def check(label, got, want, tol=1e-6):
    assert abs(got - want) < tol, f'{label}: recomputed net {got:.6f} vs stored {want:.6f}'


stored1 = pd.read_csv(f'{B}/results/gate1_intervals.tsv', sep='\t').query("analysis == 'primary'")
stored3 = pd.read_csv(f'{B}/results/gate3_arm_drops.tsv', sep='\t').query("variant == 'V0'")
tables, arms = [], []

for sp in ['Mouse', 'Cow', 'Pig', 'Rabbit']:
    X, stage, order = load_gse225056(sp)
    T = terms_matrix(X)
    tab = interval_table(T, stage, order, 'GSE225056', sp, zga=ZGA[sp])
    for _, r in tab.iterrows():
        check(f'{sp} {r.interval}', r.D, float(stored1[(stored1.species == sp) & (stored1.interval == r.interval)].delta.iloc[0]))
    tables.append(tab)
    print(f'{sp}: nets match Gate 1; gross max at {tab.loc[tab.G.idxmax(), "interval"]} (P_max {tab.G.max() and tab.loc[tab.G.idxmax(), "P_max_G"]:.2f})')

X, stage, order = load_gse45719()
tab = interval_table(terms_matrix(X), stage, order, 'GSE45719', 'Mouse')
e = pd.read_csv(f'{B}/results/gate2b_R1_embryo_tage.tsv', sep='\t', index_col=0)
m = e.groupby('stage').tAge_primary.mean()
check('GSE45719 Early-2C->Mid-2C', tab.set_index('interval').loc['Early-2C->Mid-2C', 'D'], m['Mid-2C'] - m['Early-2C'])
tables.append(tab)

X, stage, order = load_gse66582()
tab = interval_table(terms_matrix(X), stage, order, 'GSE66582', 'Mouse')
s2 = pd.read_csv(f'{B}/results/gate2b_R2_intervals.tsv', sep='\t').query("analysis == 'primary'").set_index('interval').delta
for _, r in tab.iterrows():
    check(f'GSE66582 {r.interval}', r.D, float(s2[r.interval]))
tables.append(tab)

X, meta, order = load_gse36552()
tab = interval_table(terms_matrix(X), meta, order, 'GSE36552', 'Human', zga=ZGA['Human'])
check('GSE36552 8-cell->Morula', tab.set_index('interval').loc['8-cell->Morula', 'D'], -0.098, tol=5e-4)
tables.append(tab)

for gse, pert in [('GSE280522', ['A485', 'A485+Dux']), ('GSE221985', ['Tardbp_matKO']), ('GSE300734', ['Brg1_matKO'])]:
    X, M = load_perturbation(gse)
    T = terms_matrix(X)
    groups = {(a, s): list(M.index[(M.arm == a) & (M.stage == s)]) for a in M.arm.unique() for s in ['E2C', 'L2C']}
    pairs = [((a, 'E2C'), (a, 'L2C')) for a in ['control'] + pert]
    Db, Gb = budget_boot(T, groups, pairs)
    for k, a in enumerate(['control'] + pert):
        bgt = budget(T, groups[(a, 'E2C')], groups[(a, 'L2C')])
        check(f'{gse} {a}', bgt['D'], float(stored3[(stored3.gse == gse) & (stored3.arm == a)]['drop'].iloc[0]))
        row = dict(dataset=gse, arm=a, n_E2C=len(groups[(a, 'E2C')]), n_L2C=len(groups[(a, 'L2C')]), **bgt,
                   C=abs(bgt['D']) / bgt['G'], G_lo=np.percentile(Gb[:, k], 2.5), G_hi=np.percentile(Gb[:, k], 97.5))
        if a != 'control':
            dG = Gb[:, k] - Gb[:, 0]
            row.update(dG_vs_control=bgt['G'] - budget(T, groups[('control', 'E2C')], groups[('control', 'L2C')])['G'],
                       dG_lo=np.percentile(dG, 2.5), dG_hi=np.percentile(dG, 97.5))
        if a == 'A485+Dux':
            dR = Gb[:, k] - Gb[:, 1]
            row.update(dG_rescue=bgt['G'] - budget(T, groups[('A485', 'E2C')], groups[('A485', 'L2C')])['G'],
                       dG_rescue_lo=np.percentile(dR, 2.5), dG_rescue_hi=np.percentile(dR, 97.5))
        arms.append(row)
    print(f'{gse}: nets match Gate 3')

I = pd.concat(tables, ignore_index=True)
A = pd.DataFrame(arms)
I.to_csv(f'{B}/results/explore_A_netgross_intervals.tsv', sep='\t', index=False)
A.to_csv(f'{B}/results/explore_A_netgross_arms.tsv', sep='\t', index=False)

# ====================================================================== readings against the frozen rules
L = []
L.append('# Exploration A — net and gross (plan: plan/EXPLORE_A_netgross_frozen.md)\n')
L.append('All recomputed nets matched the stored Gate 1 / Gate 2b / Gate 3 / human values before any gross was written.\n')
L.append('## P1 — species: interval of maximal gross vs literature major-ZGA interval\n')
L.append('| species | gross-max interval | P_max | ZGA interval | match | net-min interval |')
L.append('|---|---|---|---|---|---|')
hits = 0
for sp in ['Mouse', 'Cow', 'Pig', 'Rabbit']:
    t = I[(I.dataset == 'GSE225056') & (I.species == sp)]
    gmax = t.loc[t.G.idxmax()]
    nmin = t.loc[t.D.idxmin()]
    ok = bool(gmax.is_zga and gmax.P_max_G >= 0.6)
    hits += ok
    L.append(f'| {sp} | {gmax.interval} | {gmax.P_max_G:.2f} | {"->".join(ZGA[sp])} | {"yes" if ok else "no"} | {nmin.interval} |')
L.append(f'\n**P1: {hits} of 4 — {"SUPPORTED" if hits >= 3 else "UNSUPPORTED"}** (net: coincidence in mouse alone).\n')
alt = pd.read_csv(f'{B}/results/gate1alt_intervals.tsv', sep='\t')
L.append('Secondary (reported only): interval of maximal *measured* zygotic gain (Gate 1-alt) vs gross-max interval\n')
L.append('| species | measured zygotic-gain max | gross max | same? |')
L.append('|---|---|---|---|')
for sp in ['Mouse', 'Cow', 'Pig', 'Rabbit']:
    a = alt[alt.species == sp]
    zmax = a.loc[a.d_zygotic.idxmax(), 'interval']
    t_ = I[(I.dataset == 'GSE225056') & (I.species == sp)]
    gmax = t_.loc[t_.G.idxmax(), 'interval']
    L.append(f'| {sp} | {zmax} | {gmax} | {"yes" if zmax == gmax else "no"} |')
L.append('')
L.append('## P2 — perturbation and rescue (gross per arm)\n')
L.append('| dataset | arm | n | D | G | down | up | C | ΔG vs control [95% CI] |')
L.append('|---|---|---|---|---|---|---|---|---|')
signs = []
for _, r in A.iterrows():
    dg = '' if pd.isna(r.get('dG_vs_control', np.nan)) else f'{r.dG_vs_control:+.3f} [{r.dG_lo:+.3f}, {r.dG_hi:+.3f}]'
    L.append(f'| {r.dataset} | {r.arm} | {int(r.n_E2C)}/{int(r.n_L2C)} | {r.D:+.3f} | {r.G:.3f} | {r.down:+.3f} | {r.up:+.3f} | {r.C:.2f} | {dg} |')
    if r.arm in ('A485', 'Tardbp_matKO', 'Brg1_matKO'):
        signs.append(r.dG_vs_control < 0)
resc = A[A.arm == 'A485+Dux'].iloc[0]
L.append(f'\nPerturbed below control in {sum(signs)} of 3 primaries — **P2 perturbation {"SUPPORTED" if all(signs) else "UNSUPPORTED"}**; '
         f'rescue G(A485+DUX) − G(A485) = {resc.dG_rescue:+.3f} [{resc.dG_rescue_lo:+.3f}, {resc.dG_rescue_hi:+.3f}] — '
         f'**{"excludes zero" if resc.dG_rescue_lo > 0 else "includes zero"}**.\n')
L.append('## P3 — human (GSE36552): gross-max interval vs 4-cell→8-cell\n')
t = I[I.dataset == 'GSE36552']
gmax, nmin = t.loc[t.G.idxmax()], t.loc[t.D.idxmin()]
L.append(f'gross max at **{gmax.interval}** (P_max {gmax.P_max_G:.2f}); net min at {nmin.interval}. '
         f'**P3 {"SUPPORTED" if gmax.is_zga and gmax.P_max_G >= 0.6 else "UNSUPPORTED"}**.\n')
L.append('## P4 — mouse replication series\n')
for ds in ['GSE45719', 'GSE66582']:
    t = I[I.dataset == ds]
    gmax = t.loc[t.G.idxmax()]
    L.append(f'- {ds}: gross max at **{gmax.interval}** (P_max {gmax.P_max_G:.2f}); net min at {t.loc[t.D.idxmin(), "interval"]}')
L.append('\n## Every interval\n')
def md_table(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(f'{v:.4f}' if isinstance(v, float) else str(v) for v in r.values) + ' |')
    return '\n'.join(out)


L.append(md_table(I))
open(f'{B}/results/explore_A_netgross.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L[:40]))
