"""Gate 2a (plan: plan/GATE2A_PLAN_frozen.md). Is the clock drop a compositional artefact of maternal clearance /
zygotic activation? Variants V0-V2 on real embryos, V3 composition-only simulation. seed 20260916."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
src = open(f'{B}/src/01_gate1_cross_species.py', encoding='utf-8').read().split('emb_rows, int_rows, summ')[0]
ns = {}
exec(src, ns)
load, qc, ORDER, CELLNUM = ns['load'], ns['qc'], ns['ORDER'], ns['CELLNUM']
rng = np.random.default_rng(20260916)
MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
SPECIES = ['Mouse', 'Cow', 'Pig', 'Rabbit']

_, feats = tp.load_clock(MODEL)
gm = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
gm = gm[gm.Entrez.notna()]
clock_ens = set(gm[gm.Entrez.astype('int64').astype(str).isin(set(feats))].Ensembl)

A = pd.read_csv(f'{B}/results/gate1alt_intervals.tsv', sep='\t')
DATA_ZGA = {sp: g.reset_index(drop=True).interval[int(g.d_zygotic.values.argmax())] for sp, g in A.groupby('species')}


def gene_sets(C, stage, order):
    cpm = C / C.sum(axis=0) * 1e6
    cols = {s: list(stage.index[stage == s]) for s in order}
    later = [s for s in order if s not in ('Oocyte', 'Zygote')]
    fge = lambda s, t: (cpm[cols[s]] >= t).mean(axis=1)
    flt = lambda s, t: (cpm[cols[s]] < t).mean(axis=1)
    oo_med = cpm[cols['Oocyte']].median(axis=1)
    mat = (fge('Oocyte', 10) >= 0.8) & (fge('Zygote', 10) >= 0.8) & \
          pd.concat([cpm[cols[s]].median(axis=1) <= oo_med / 8 for s in later], axis=1).any(axis=1)
    zyg = (flt('Oocyte', 1) >= 0.8) & (flt('Zygote', 1) >= 0.8) & \
          pd.concat([fge(s, 10) >= 0.5 for s in later], axis=1).any(axis=1)
    return set(cpm.index[mat]), set(cpm.index[zyg]), cpm


def clock(C, ref):
    pp = tp.preprocess(C, reference_samples=ref)
    return tp.predict(MODEL, pp['scaled_diff'], species=None)


def deltas(v, stage, order, nboot=2000):
    groups = [v[stage == s].values for s in order]
    d = np.diff([g.mean() for g in groups])
    boot = np.array([np.diff([rng.choice(g, len(g)).mean() for g in groups]) for _ in range(nboot)])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    return d, lo, hi


rows, sim_rows, set_rows = [], [], []
for sp in SPECIES:
    C, stage, _ = load(sp)
    keep, _ = qc(C)
    C, stage = C.loc[:, keep], stage[keep]
    order = ORDER[sp]
    mat, zyg, cpm = gene_sets(C, stage, order)
    D = mat | zyg
    ref = list(stage.index[stage == 'Oocyte'])
    set_rows.append(dict(species=sp, maternal=len(mat), zygotic=len(zyg), dynamic=len(D),
                         dynamic_clock_features=len(D & clock_ens), clock_features_in_matrix=len(set(C.index) & clock_ens)))
    variants = {'V0_original': C,
                'V1_drop_nonclock_dynamic': C.loc[~C.index.isin(D - clock_ens)],
                'V2_drop_all_dynamic': C.loc[~C.index.isin(D)]}
    for vname, Cv in variants.items():
        v = clock(Cv, ref)
        d, lo, hi = deltas(v, stage, order)
        for k, (a, c) in enumerate(zip(order[:-1], order[1:])):
            rows.append(dict(species=sp, variant=vname, interval=f'{a}->{c}', cell_number=CELLNUM[(a, c)], delta=d[k],
                             ci_lo=lo[k], ci_hi=hi[k], data_zga=bool(f'{a}->{c}' == DATA_ZGA[sp]), n_genes=Cv.shape[0]))
    # V3 composition-only simulation on real oocytes
    a, c = DATA_ZGA[sp].split('->')
    share = cpm.loc[list(mat)].sum(axis=0) / 1e6
    r = share[stage == c].mean() / share[stage == a].mean()
    oo = C[ref]
    sim = oo.copy()
    sim.loc[sim.index.isin(mat)] = sim.loc[sim.index.isin(mat)] * r
    sim.columns = ['sim_' + x for x in sim.columns]
    v = clock(pd.concat([oo, sim], axis=1), ref)
    vo, vs = v[ref].values, v[sim.columns].values
    boot = [rng.choice(vs, len(vs)).mean() - rng.choice(vo, len(vo)).mean() for _ in range(2000)]
    sim_rows.append(dict(species=sp, data_zga=DATA_ZGA[sp], maternal_share_before=share[stage == a].mean(),
                         maternal_share_after=share[stage == c].mean(), ratio=r, delta_sim=vs.mean() - vo.mean(),
                         ci_lo=np.percentile(boot, 2.5), ci_hi=np.percentile(boot, 97.5)))

R = pd.DataFrame(rows)
Sm = pd.DataFrame(sim_rows)
R.to_csv(f'{B}/results/gate2a_intervals.tsv', sep='\t', index=False)
Sm.to_csv(f'{B}/results/gate2a_simulation.tsv', sep='\t', index=False)
pd.DataFrame(set_rows).to_csv(f'{B}/results/gate2a_gene_sets.tsv', sep='\t', index=False)

# ---------------- decision (mouse, cow)
res = {}
for sp in ['Mouse', 'Cow']:
    z = R[(R.species == sp) & R.data_zga].set_index('variant')
    d0, d1, d2 = z.loc['V0_original', 'delta'], z.loc['V1_drop_nonclock_dynamic', 'delta'], z.loc['V2_drop_all_dynamic', 'delta']
    ds = Sm.set_index('species').loc[sp, 'delta_sim']
    res[sp] = dict(interval=DATA_ZGA[sp], d_V0=d0, d_V1=d1, d_V2=d2, V2_ci_hi=z.loc['V2_drop_all_dynamic', 'ci_hi'],
                   V2_ci_lo=z.loc['V2_drop_all_dynamic', 'ci_lo'], retention_V1=d1 / d0, retention_V2=d2 / d0,
                   d_sim=ds, sim_fraction_of_V0=ds / d0)
robust = all(r['d_V2'] < 0 and r['V2_ci_hi'] < 0 and r['retention_V2'] >= 0.5 and abs(r['d_sim']) < 0.5 * abs(r['d_V0'])
             for r in res.values())
comp_a = all((r['V2_ci_lo'] <= 0 <= r['V2_ci_hi']) or r['retention_V2'] < 0.25 for r in res.values())
comp_b = all(r['sim_fraction_of_V0'] >= 0.75 for r in res.values())
verdict = 'RESET-ROBUST' if robust else ('COMPOSITIONAL' if (comp_a or comp_b) else 'INTERMEDIATE')

pd.set_option('display.width', 250)
print(pd.DataFrame(set_rows).to_string(index=False))
print(R[R.data_zga][['species', 'variant', 'interval', 'delta', 'ci_lo', 'ci_hi', 'n_genes']].round(4).to_string(index=False))
print(Sm.round(4).to_string(index=False))
print(pd.DataFrame(res).T.round(3).to_string())
print('VERDICT:', verdict)
with open(f'{B}/results/gate2a_decision.md', 'w', encoding='utf-8') as fh:
    fh.write('# Gate 2a decision (plan/GATE2A_PLAN_frozen.md)\n\n')
    fh.write(f'VERDICT: **{verdict}**\n\n')
    fh.write(pd.DataFrame(res).T.round(3).to_string() + '\n\nSimulation:\n' + Sm.round(4).to_string(index=False) + '\n')
