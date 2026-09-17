"""Gate 1 (plan: plan/GATE1_PLAN_frozen.md). Cross-species timing of the transcriptomic-age drop vs major ZGA.
seed 20260914. Exploratory."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys, re, json
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
M = f'{B}/tools_tAge_models'
rng = np.random.default_rng(20260914)
PRIMARY = 'EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
SENS = ['EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl', 'EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl',
        'EN_Chronoage_Mouse_Multitissue_scaleddiff.pkl']
ORDER = {'Mouse': ['Oocyte', 'Zygote', 'Early-2-cell', 'Late-2-cell', '4-cell', '8-cell', '16-cell'],
         'Cow': ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
         'Rabbit': ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
         'Pig': ['Oocyte', 'Zygote', '2-cell', '4-cell', 'Day2', 'Day3', 'Morula'],
         'Rhesus': ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell']}
CELLNUM = {('Oocyte', 'Zygote'): 'fertilisation', ('Zygote', '2-cell'): '1->2', ('Zygote', 'Early-2-cell'): '1->2',
           ('Early-2-cell', 'Late-2-cell'): 'within-2-cell', ('Late-2-cell', '4-cell'): '2->4', ('2-cell', '4-cell'): '2->4',
           ('4-cell', '8-cell'): '4->8', ('4-cell', 'Day2'): '4->8', ('8-cell', '16-cell'): '8->16', ('Day2', 'Day3'): '8->16',
           ('16-cell', 'Morula'): '16->morula', ('Day3', 'Morula'): '16->morula'}
ZGA = {'Mouse': ('Early-2-cell', 'Late-2-cell'), 'Pig': ('2-cell', '4-cell'), 'Cow': ('8-cell', '16-cell'),
       'Rabbit': ('8-cell', '16-cell')}
ORTHO = {'Cow': 'btaurus', 'Pig': 'sscrofa', 'Rabbit': 'ocuniculus'}

gm = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')


def load(sp):
    C = pd.read_csv(f'{B}/data/GSE225056_{sp}_count_5prime.txt.gz', sep='\t', index_col=0)
    C.index = C.index.str.strip('"')
    C.columns = C.columns.str.strip('"')
    if sp == 'Cow':
        # addendum 1: RefSeq transcripts -> cow Ensembl gene (NM_ via Ensembl xref; XM_ via NCBI title symbol)
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
    raw_genes = C.shape[0]
    if sp in ORTHO:
        O = pd.read_csv(f'{B}/metadata/orthologs_{ORTHO[sp]}_to_mouse.tsv', sep='\t')
        O.columns = ['gid', 'name', 'mouse_gid', 'type']
        O = O[O.type == 'ortholog_one2one'].drop_duplicates('gid')
        mp = dict(zip(O.gid, O.mouse_gid))
        ids = C.index.to_series().map(mp)
        ok = ids.notna()
        C = C.loc[ok].groupby(ids[ok].values).sum()
    elif sp == 'Rhesus':
        # addendum 1: M. mulatta Ensembl one-to-one orthologs (package table is M. fascicularis)
        O = pd.read_csv(f'{B}/metadata/orthologs_mmulatta_to_mouse.tsv', sep='\t')
        O.columns = ['gid', 'name', 'mouse_gid', 'type']
        O = O[O.type == 'ortholog_one2one'].drop_duplicates('gid')
        mp = dict(zip(O.gid, O.mouse_gid))
        ids = C.index.to_series().map(mp)
        ok = ids.notna()
        C = C.loc[ok].groupby(ids[ok].values).sum()
    stage = pd.Series([re.sub(r'_\d+$', '', c) for c in C.columns], index=C.columns)
    return C, stage, raw_genes


def qc(C):
    det = np.log10((C >= 1).sum(axis=0))
    med = det.median()
    mad = 1.4826 * np.median(np.abs(det - med))
    return det >= med - 3 * mad, det


def score(C, stage, model):
    pp = tp.preprocess(C, reference_samples=list(stage.index[stage == 'Oocyte']))
    key = 'yugene_diff' if 'yugenediff' in model else 'scaled_diff'
    return tp.predict(f'{M}/{model}', pp[key], species=None), tp.feature_coverage(f'{M}/{model}', pp[key])


def intervals(v, stage, order, nboot=2000, zga=None):
    stages = [s for s in order if (stage == s).any()]
    groups = [v[stage == s].values for s in stages]
    means = np.array([g.mean() for g in groups])
    pairs = list(zip(stages[:-1], stages[1:]))
    d = np.diff(means)
    boot = np.empty((nboot, len(d)))
    for b in range(nboot):
        bm = np.array([rng.choice(g, len(g)).mean() for g in groups])
        boot[b] = np.diff(bm)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    amin = np.argmin(boot, axis=1)
    rows = []
    for k, (a, c) in enumerate(pairs):
        rows.append(dict(interval=f'{a}->{c}', cell_number=CELLNUM.get((a, c), '?'), delta=d[k], ci_lo=lo[k],
                         ci_hi=hi[k], ci_excludes_0=bool(lo[k] > 0 or hi[k] < 0),
                         p_boot_most_negative=float(np.mean(amin == k)),
                         is_zga=bool(zga is not None and (a, c) == zga)))
    return pd.DataFrame(rows), {s: float(m) for s, m in zip(stages, means)}


emb_rows, int_rows, summ = [], [], {}
for sp in ['Mouse', 'Cow', 'Pig', 'Rabbit', 'Rhesus']:
    C, stage, raw_genes = load(sp)
    keep, det = qc(C)
    C, stage = C.loc[:, keep], stage[keep]
    missing = [s for s in ORDER[sp] if not (stage == s).any()]
    if missing:
        print(sp, 'not scored: stages missing after QC', missing)
        continue
    tot = np.log10(C.sum(axis=0))
    res = {}
    for model in [PRIMARY] + SENS:
        if 'Mouse_Multitissue' in model and sp != 'Mouse':
            continue
        v, cov = score(C, stage, model)
        res[model] = v
        lab = 'primary' if model == PRIMARY else model.replace('.pkl', '')
        I, means = intervals(v, stage, ORDER[sp], zga=ZGA.get(sp))
        I.insert(0, 'analysis', lab)
        I.insert(0, 'species', sp)
        I['feature_coverage'] = cov
        I['n_genes_mapped'] = C.shape[0]
        I['n_genes_raw'] = raw_genes
        int_rows.append(I)
        if model == PRIMARY:
            X = np.column_stack([np.ones(len(v)), tot.values, det[keep].values])
            beta, *_ = np.linalg.lstsq(X, v.values, rcond=None)
            vr = pd.Series(v.values - X @ beta, index=v.index)
            Ir, _ = intervals(vr, stage, ORDER[sp], zga=ZGA.get(sp))
            Ir.insert(0, 'analysis', 'primary_residualised_depth_detected')
            Ir.insert(0, 'species', sp)
            int_rows.append(Ir)
            summ[sp] = dict(n_libraries=int(keep.sum()), n_excluded_qc=int((~keep).sum()), feature_coverage=cov,
                            n_genes_mapped=int(C.shape[0]), stage_n={k: int(x) for k, x in stage.value_counts().items()},
                            stage_means=means)
    E = pd.DataFrame({'species': sp, 'stage': stage, 'log10_total_counts': tot, 'log10_detected': det[keep]})
    for model, v in res.items():
        E[model.replace('.pkl', '')] = v
    emb_rows.append(E)

E = pd.concat(emb_rows)
I = pd.concat(int_rows)
E.to_csv(f'{B}/results/gate1_embryo_tage.tsv', sep='\t')
I.to_csv(f'{B}/results/gate1_intervals.tsv', sep='\t', index=False)
json.dump(summ, open(f'{B}/results/gate1_summary.json', 'w'), indent=1)

# ---------------- decision (primary analysis, 4 primary species)
def decide(P, need=3):
    informative = {sp: bool(g.ci_excludes_0.any()) for sp, g in P.groupby('species')}
    zga_hit, most_neg_cell = {}, {}
    for sp, g in P.groupby('species'):
        r = g.iloc[int(g.delta.values.argmin())]
        most_neg_cell[sp] = r.cell_number
        zga_hit[sp] = bool(r.is_zga and r.ci_hi < 0 and g[g.is_zga].p_boot_most_negative.iloc[0] >= 0.6)
    early = zga_hit.get('Mouse', False) or zga_hit.get('Pig', False)
    late = zga_hit.get('Cow', False) or zga_hit.get('Rabbit', False)
    cn = pd.Series(most_neg_cell).value_counts()
    if sum(informative.values()) < need:
        v = 'UNINFORMATIVE'
    elif sum(zga_hit.values()) >= need and early and late:
        v = 'ZGA-LINKED'
    elif cn.iloc[0] >= need:
        v = 'DIVISION-LINKED'
    else:
        v = 'MIXED'
    return v, informative, most_neg_cell, zga_hit


P = I[(I.analysis == 'primary') & I.species.isin(['Mouse', 'Cow', 'Pig', 'Rabbit'])]
verdict, informative, most_neg_cell, zga_hit = decide(P)
sens_rows = []
for lab in sorted(I.analysis.unique()):
    Q = I[(I.analysis == lab) & I.species.isin(['Mouse', 'Cow', 'Pig', 'Rabbit'])]
    if Q.species.nunique() >= 3:
        sens_rows.append(dict(analysis=lab, species=','.join(sorted(Q.species.unique())), verdict=decide(Q)[0]))
Qnc = P[P.species != 'Cow']
sens_rows.append(dict(analysis='primary_cow_excluded (3 of 3 required)', species=','.join(sorted(Qnc.species.unique())),
                      verdict=decide(Qnc, need=3)[0] if Qnc.species.nunique() == 3 else 'n/a'))
S = pd.DataFrame(sens_rows)
S.to_csv(f'{B}/results/gate1_sensitivity_verdicts.tsv', sep='	', index=False)

pd.set_option('display.width', 250)
cols = ['species', 'interval', 'cell_number', 'delta', 'ci_lo', 'ci_hi', 'p_boot_most_negative', 'is_zga']
print(P[cols].round(4).to_string(index=False))
print('\ninformative:', informative)
print('most negative interval (cell number):', most_neg_cell)
print('ZGA criterion met:', zga_hit)
print('VERDICT:', verdict)
print(S.to_string(index=False))
for sp, s in summ.items():
    print(sp, 'libs', s['n_libraries'], 'qc-excluded', s['n_excluded_qc'], 'genes mapped', s['n_genes_mapped'],
          'clock coverage %.3f' % s['feature_coverage'])
with open(f'{B}/results/gate1_decision.md', 'w', encoding='utf-8') as fh:
    fh.write('# Gate 1 decision (primary analysis, frozen rules)\n\n')
    fh.write(f'VERDICT: **{verdict}**\n\n')
    fh.write(f'- informative species: {informative}\n')
    fh.write(f'- most negative interval (cell-number label): {most_neg_cell}\n')
    fh.write(f'- ZGA-interval criterion met: {zga_hit}\n\n')
    fh.write('Intervals and sensitivity analyses: results/gate1_intervals.tsv\n')
