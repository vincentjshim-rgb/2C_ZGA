"""Gate 1 descriptive figure: transcriptomic age across oocyte-to-morula stages in four species.

Each panel = one species (ordered by the stage of major ZGA). Grey dots = single embryos; black = stage mean with
95% bootstrap CI (embryos resampled within stage, 2000 replicates, seed 20260914). Shaded band = major ZGA
interval fixed from the literature before scoring (plan/GATE1_PLAN_frozen.md). y = tAge from
EN_Chronoage_Multispecies_Multitissue_scaleddiff, centred on the species' oocytes (normalised units).
Source data: results/gate1_embryo_tage.tsv. Descriptive only; the decision is in results/gate1_decision.md.
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
B = f'{ANALYSIS}/gz'
rng = np.random.default_rng(20260914)
COL = 'EN_Chronoage_Multispecies_Multitissue_scaleddiff'
E = pd.read_csv(f'{B}/results/gate1_embryo_tage.tsv', sep='\t', index_col=0)
I = pd.read_csv(f'{B}/results/gate1_intervals.tsv', sep='\t')
I = I[I.analysis == 'primary']

PANELS = [('Mouse', ['Oocyte', 'Zygote', 'Early-2-cell', 'Late-2-cell', '4-cell', '8-cell', '16-cell'],
           ['Oo', 'Zy', 'E2C', 'L2C', '4C', '8C', '16C'], ('Early-2-cell', 'Late-2-cell')),
          ('Pig', ['Oocyte', 'Zygote', '2-cell', '4-cell', 'Day2', 'Day3', 'Morula'],
           ['Oo', 'Zy', '2C', '4C', '5-8C', '8-16C', 'Mor'], ('2-cell', '4-cell')),
          ('Cow', ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
           ['Oo', 'Zy', '2C', '4C', '8C', '16C', 'Mor'], ('8-cell', '16-cell')),
          ('Rabbit', ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
           ['Oo', 'Zy', '2C', '4C', '8C', '16C', 'Mor'], ('8-cell', '16-cell'))]

fig, axes = plt.subplots(1, 4, figsize=(fs.W2, 58 * fs.MM), sharey=True,
                         gridspec_kw=dict(wspace=0.12, left=0.075, right=0.99, top=0.84, bottom=0.2))
rows = []
for j, (ax, (sp, order, labels, zga)) in enumerate(zip(axes, PANELS)):
    d = E[E.species == sp]
    x0, x1 = order.index(zga[0]), order.index(zga[1])
    ax.axvspan(x0, x1, color=BAND, lw=0, zorder=0)
    ax.text((x0 + x1) / 2, 1.0, 'major ZGA', transform=ax.get_xaxis_transform(), ha='center', va='bottom',
            fontsize=5.8, color=INK2)
    means, los, his = [], [], []
    for k, st in enumerate(order):
        v = d.loc[d.stage == st, COL].values
        ax.scatter(k + (rng.random(len(v)) - 0.5) * 0.35, v, s=3, color=MUTED, edgecolor='none', alpha=0.7, zorder=2)
        bm = np.array([rng.choice(v, len(v)).mean() for _ in range(2000)])
        means.append(v.mean()); los.append(np.percentile(bm, 2.5)); his.append(np.percentile(bm, 97.5))
        rows.append(dict(species=sp, stage=st, n=len(v), mean=v.mean(), ci_lo=los[-1], ci_hi=his[-1]))
    ax.plot(range(len(order)), means, color=INK, lw=0.9, zorder=3)
    ax.vlines(range(len(order)), los, his, color=INK, lw=1.0, zorder=3)
    ax.scatter(range(len(order)), means, s=10, color=INK, zorder=4)
    ax.axhline(0, color=GRID, lw=0.6, zorder=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontsize=5.8)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.6, len(order) - 0.4)
    g = I[I.species == sp]
    r = g.iloc[int(g.delta.values.argmin())]
    ax.set_title(sp, loc='left', fontsize=7, color=INK, pad=11)
    ax.text(0.02, 0.03, 'largest drop: %s' % r.interval.replace('->', ' to '), transform=ax.transAxes,
            fontsize=5.5, color=INK2, va='bottom')
    ax.grid(axis='y', color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    fs.panel(ax, 'abcd'[j], dx=-0.14 if j == 0 else -0.06, dy=1.2)
axes[0].set_ylabel('Transcriptomic age relative to oocyte\n(multispecies clock, normalised)', color=INK2)
fig.text(0.53, 0.035, 'Stage (ordered); panels ordered by stage of major ZGA', ha='center', fontsize=6.3, color=INK2)
fs.save(fig, 'FigGZ_1_gate1_trajectories', outdir=f'{B}/figures')
pd.DataFrame(rows).to_csv(f'{B}/results/figgz1_source_stage_means.tsv', sep='\t', index=False)
