"""Gate 1-alt (plan: plan/GATE1ALT_PLAN_frozen.md). POST HOC / EXPLORATORY.
Data-driven ZGA timing (zygotic-gene share, clock genes excluded) vs the Gate 1 clock drop, four species.
seed 20260915."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
src = open(f'{B}/src/01_gate1_cross_species.py', encoding='utf-8').read().split('emb_rows, int_rows, summ')[0]
ns = {}
exec(src, ns)                       # reuse Gate 1 loading, QC, ORDER, CELLNUM (no scoring runs here)
load, qc, ORDER, CELLNUM, ZGA_LIT = ns['load'], ns['qc'], ns['ORDER'], ns['CELLNUM'], ns['ZGA']
rng = np.random.default_rng(20260915)
SPECIES = ['Mouse', 'Pig', 'Cow', 'Rabbit']

_, feats = tp.load_clock(f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl')
gm = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
gm = gm[gm.Entrez.notna()]
clock_ens = set(gm[gm.Entrez.astype('int64').astype(str).isin(set(feats))].Ensembl)

G = pd.read_csv(f'{B}/results/gate1_intervals.tsv', sep='\t')
G = G[G.analysis == 'primary']


def frac_ge(M, cols, thr):
    return (M[cols] >= thr).mean(axis=1)


def frac_lt(M, cols, thr):
    return (M[cols] < thr).mean(axis=1)


def delta_table(score, stage, order, nboot=2000):
    groups = [score[stage == s].values for s in order]
    means = np.array([g.mean() for g in groups])
    d = np.diff(means)
    boot = np.empty((nboot, len(d)))
    for b in range(nboot):
        boot[b] = np.diff([rng.choice(g, len(g)).mean() for g in groups])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    return means, d, lo, hi, boot


stage_rows, int_rows, set_rows = [], [], []
for sp in SPECIES:
    C, stage, _ = load(sp)
    keep, _ = qc(C)
    C, stage = C.loc[:, keep], stage[keep]
    C = C.loc[~C.index.isin(clock_ens)]
    cpm = C / C.sum(axis=0) * 1e6
    order = ORDER[sp]
    cols = {s: list(stage.index[stage == s]) for s in order}
    early = ['Oocyte', 'Zygote']
    later = [s for s in order if s not in early]
    zyg = (frac_lt(cpm, cols['Oocyte'], 1) >= 0.8) & (frac_lt(cpm, cols['Zygote'], 1) >= 0.8) & \
          pd.concat([frac_ge(cpm, cols[s], 10) >= 0.5 for s in later], axis=1).any(axis=1)
    oo_med = cpm[cols['Oocyte']].median(axis=1)
    mat = (frac_ge(cpm, cols['Oocyte'], 10) >= 0.8) & (frac_ge(cpm, cols['Zygote'], 10) >= 0.8) & \
          pd.concat([cpm[cols[s]].median(axis=1) <= oo_med / 8 for s in later], axis=1).any(axis=1)
    A = cpm.loc[zyg].sum(axis=0) / 1e6
    Bm = cpm.loc[mat].sum(axis=0) / 1e6
    set_rows += [dict(species=sp, set='zygotic', n_genes=int(zyg.sum()), n_genes_considered=int(len(cpm))),
                 dict(species=sp, set='maternal', n_genes=int(mat.sum()), n_genes_considered=int(len(cpm)))]
    mA, dA, loA, hiA, bA = delta_table(A, stage, order)
    mB, dB, loB, hiB, bB = delta_table(Bm, stage, order)
    pA = np.bincount(bA.argmax(axis=1), minlength=len(dA)) / len(bA)
    for k, s in enumerate(order):
        stage_rows.append(dict(species=sp, stage=s, n=len(cols[s]), zygotic_share=mA[k], maternal_share=mB[k]))
    g = G[G.species == sp].reset_index(drop=True)
    for k, (a, c) in enumerate(zip(order[:-1], order[1:])):
        assert g.interval[k] == f'{a}->{c}'
        int_rows.append(dict(species=sp, interval=f'{a}->{c}', cell_number=CELLNUM[(a, c)],
                             d_zygotic=dA[k], d_zygotic_lo=loA[k], d_zygotic_hi=hiA[k], p_boot_max_zygotic=pA[k],
                             d_maternal=dB[k], d_maternal_lo=loB[k], d_maternal_hi=hiB[k],
                             d_clock=g.delta[k], d_clock_lo=g.ci_lo[k], d_clock_hi=g.ci_hi[k],
                             literature_zga=bool((a, c) == ZGA_LIT[sp])))

S = pd.DataFrame(stage_rows)
I = pd.DataFrame(int_rows)
S.to_csv(f'{B}/results/gate1alt_stage_scores.tsv', sep='\t', index=False)
I.to_csv(f'{B}/results/gate1alt_intervals.tsv', sep='\t', index=False)
pd.DataFrame(set_rows).to_csv(f'{B}/results/gate1alt_gene_sets.tsv', sep='\t', index=False)

# ---------------- decision
summary, match = {}, 0
for sp in SPECIES:
    g = I[I.species == sp].reset_index(drop=True)
    kz, kc, km = int(g.d_zygotic.values.argmax()), int(g.d_clock.values.argmin()), int(g.d_maternal.values.argmin())
    summary[sp] = dict(data_zga=g.interval[kz], data_zga_cell=g.cell_number[kz], p_boot=round(float(g.p_boot_max_zygotic[kz]), 3),
                       clock_drop=g.interval[kc], maternal_clearance=g.interval[km],
                       literature_zga=g.interval[g.literature_zga].iloc[0], clock_matches_data_zga=bool(kz == kc),
                       clock_matches_maternal=bool(km == kc), data_matches_literature=bool(bool(g.literature_zga[kz])))
    match += int(kz == kc)
separable = len({v['data_zga_cell'] for v in summary.values()}) >= 2


def pooled_rho(Ix):
    return pd.Series(Ix.d_zygotic.values).rank().corr(pd.Series(-Ix.d_clock.values).rank())


rho = pooled_rho(I)
null = np.empty(10000)
for b in range(10000):
    parts = []
    for sp in SPECIES:
        g = I[I.species == sp]
        parts.append(pd.DataFrame({'d_zygotic': g.d_zygotic.values, 'd_clock': rng.permutation(g.d_clock.values)}))
    null[b] = pooled_rho(pd.concat(parts))
p_perm = (1 + np.sum(null >= rho)) / (1 + len(null))
cond1, cond2 = match >= 3, (rho > 0 and p_perm < 0.05)
if not separable:
    verdict = 'NOT SEPARABLE'
elif cond1 and cond2:
    verdict = 'ALIGNED'
elif cond1 or cond2:
    verdict = 'PARTIAL'
else:
    verdict = 'NOT ALIGNED'

pd.set_option('display.width', 250)
print(pd.DataFrame(set_rows).to_string(index=False))
print(I[['species', 'interval', 'd_zygotic', 'd_zygotic_lo', 'd_zygotic_hi', 'p_boot_max_zygotic', 'd_maternal', 'd_clock',
         'literature_zga']].round(4).to_string(index=False))
print(pd.DataFrame(summary).T.to_string())
print('clock matches data ZGA: %d/4 | separable %s | pooled Spearman %.3f, permutation p %.4f' % (match, separable, rho, p_perm))
print('VERDICT (post hoc):', verdict)
with open(f'{B}/results/gate1alt_decision.md', 'w', encoding='utf-8') as fh:
    fh.write('# Gate 1-alt decision (POST HOC / EXPLORATORY; frozen rules in plan/GATE1ALT_PLAN_frozen.md)\n\n')
    fh.write(f'VERDICT: **{verdict}**\n\n- clock drop = data-driven ZGA interval in {match}/4 species\n')
    fh.write(f'- ZGA timing separable from cell number: {separable}\n')
    fh.write(f'- pooled Spearman(dZygotic, -dClock) = {rho:.3f}, within-species permutation p = {p_perm:.4f}\n\n')
    fh.write(pd.DataFrame(summary).T.to_string() + '\n')
