"""Gate 2 figure. a: independent mouse replication (GSE45719, embryo pseudobulk), transcriptomic age by stage,
primary and dynamic-gene-removed (V2). b: composition test in GSE225056 (mouse, cow): clock change in the data-driven
ZGA interval under V0/V1/V2 and the composition-only simulation, 95% bootstrap CI.
Sources: results/gate2b_R1_embryo_tage.tsv, results/gate2a_intervals.tsv, results/gate2a_simulation.tsv. seed 20260917.
"""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, f'{ANALYSIS}/src')
import figstyle as fs

fs.use()
INK, INK2, MUTED, GRID, BAND = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#f0e6cf'
BLUE, ORANGE = '#2a78d6', '#eb6834'
B = f'{ANALYSIS}/gz'
rng = np.random.default_rng(20260917)
E = pd.read_csv(f'{B}/results/gate2b_R1_embryo_tage.tsv', sep='\t', index_col=0)
R = pd.read_csv(f'{B}/results/gate2a_intervals.tsv', sep='\t')
S = pd.read_csv(f'{B}/results/gate2a_simulation.tsv', sep='\t').set_index('species')
ORDER = ['Zygote', 'Early-2C', 'Mid-2C', 'Late-2C', '4C', '8C', '16C', 'EarlyBl', 'MidBl', 'LateBl']
LAB = ['Zy', 'E2C', 'M2C', 'L2C', '4C', '8C', '16C', 'EBl', 'MBl', 'LBl']

fig = plt.figure(figsize=(fs.W2, 80 * fs.MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.32, left=0.08, right=0.985, top=0.86, bottom=0.3)

ax = fig.add_subplot(gs[0])
ax.axvspan(0.8, 3.2, color=BAND, lw=0, zorder=0)
ax.text(2.0, 1.01, 'within 2-cell (major ZGA)', transform=ax.get_xaxis_transform(), ha='center', va='bottom',
        fontsize=5.8, color=INK2)
for col, colr, dx, lab in [('tAge_primary', INK, -0.12, 'all genes'), ('tAge_V2', BLUE, 0.12, 'maternal/zygotic genes removed')]:
    means, lo, hi = [], [], []
    for k, st in enumerate(ORDER):
        v = E.loc[E.stage == st, col].values
        ax.scatter(k + dx + (rng.random(len(v)) - 0.5) * 0.12, v, s=4, color=colr, alpha=0.35, edgecolor='none', zorder=2)
        bm = np.array([rng.choice(v, len(v)).mean() for _ in range(2000)])
        means.append(v.mean()); lo.append(np.percentile(bm, 2.5)); hi.append(np.percentile(bm, 97.5))
    x = np.arange(len(ORDER)) + dx
    ax.plot(x, means, color=colr, lw=0.9, zorder=3, label=lab)
    ax.vlines(x, lo, hi, color=colr, lw=0.9, zorder=3)
    ax.scatter(x, means, s=9, color=colr, zorder=4)
ax.axhline(0, color=GRID, lw=0.6)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(LAB, fontsize=6)
ax.tick_params(axis='x', length=0)
n = E.groupby('stage').size().reindex(ORDER)
for k, v in enumerate(n.values):
    ax.text(k, -0.2, 'n=%d' % v, transform=ax.get_xaxis_transform(), ha='center', fontsize=5, color=MUTED)
ax.set_ylabel('Transcriptomic age relative to zygote\n(multispecies clock)', color=INK2)
ax.set_title('Independent mouse series (GSE45719), one point per embryo', loc='left', fontsize=6.5, color=INK2, pad=12)
ax.legend(frameon=False, fontsize=5.8, loc='upper right')
ax.grid(axis='y', color=GRID, lw=0.4)
ax.set_axisbelow(True)
fs.panel(ax, 'a', dx=-0.12, dy=1.2)

ax = fig.add_subplot(gs[1])
labels = ['original', 'non-clock\ndynamic\nremoved', 'all\ndynamic\nremoved', 'maternal\nclearance\nsimulated']
for j, (sp, colr) in enumerate([('Mouse', BLUE), ('Cow', ORANGE)]):
    z = R[(R.species == sp) & R.data_zga].set_index('variant')
    vals = [z.loc[v, ['delta', 'ci_lo', 'ci_hi']].values for v in ['V0_original', 'V1_drop_nonclock_dynamic', 'V2_drop_all_dynamic']]
    vals.append(S.loc[sp, ['delta_sim', 'ci_lo', 'ci_hi']].values)
    x = np.arange(4) + (j - 0.5) * 0.34
    d = np.array([v[0] for v in vals]); lo = np.array([v[1] for v in vals]); hi = np.array([v[2] for v in vals])
    ax.bar(x, d, width=0.32, color=colr, alpha=[1, 1, 1, 0.45][0] if False else 1, zorder=2,
           label='%s (%s)' % (sp, z.index.size and R[(R.species == sp) & R.data_zga].interval.iloc[0].replace('->', ' to ')))
    ax.bars = None
    ax.vlines(x, lo, hi, color=INK, lw=0.8, zorder=3)
ax.axhline(0, color=INK2, lw=0.5)
ax.set_xticks(range(4))
ax.set_xticklabels(labels, fontsize=5.4)
ax.tick_params(axis='x', length=0)
ax.set_ylabel('Clock change in ZGA interval', color=INK2)
ax.set_title('Composition test (GSE225056)', loc='left', fontsize=6.5, color=INK2, pad=12)
ax.legend(frameon=False, fontsize=5.4, loc='upper center', bbox_to_anchor=(0.5, -0.36), ncol=1)
ax.grid(axis='y', color=GRID, lw=0.4)
ax.set_axisbelow(True)
fs.panel(ax, 'b', dx=-0.22, dy=1.2)
fs.save(fig, 'FigGZ_3_gate2_replication_composition', outdir=f'{B}/figures')
