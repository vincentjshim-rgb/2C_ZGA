"""Gate 1-alt descriptive figure (POST HOC). For each species, stage-to-stage change in zygotic-gene share (data-driven
ZGA, clock genes excluded) and in maternal-gene share, against the change in transcriptomic age (Gate 1 primary clock).

a-d  per species: bars = change per interval (top: zygotic share, up = activation; middle: maternal share,
     down = clearance; bottom: transcriptomic age, down = younger), 95% bootstrap CI.
     Shaded interval = data-driven ZGA (largest zygotic gain); dashed outline = literature ZGA.
e    pooled intervals: change in zygotic share (rank within species) vs change in transcriptomic age.
Source data: results/gate1alt_intervals.tsv. seed 20260915.
"""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, f'{ANALYSIS}/src')
import figstyle as fs

fs.use()
INK, INK2, MUTED, GRID, BAND = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#f0e6cf'
COLS = {'Mouse': '#2a78d6', 'Pig': '#eb6834', 'Cow': '#1b9e77', 'Rabbit': '#7b5ea7'}
B = f'{ANALYSIS}/gz'
I = pd.read_csv(f'{B}/results/gate1alt_intervals.tsv', sep='\t')
SHORT = {'Oocyte': 'Oo', 'Zygote': 'Zy', 'Early-2-cell': 'E2C', 'Late-2-cell': 'L2C', '2-cell': '2C', '4-cell': '4C',
         '8-cell': '8C', '16-cell': '16C', 'Morula': 'Mor', 'Day2': '5-8C', 'Day3': '8-16C'}
SPECIES = ['Mouse', 'Pig', 'Cow', 'Rabbit']

fig = plt.figure(figsize=(fs.W2, 128 * fs.MM))
outer = fig.add_gridspec(1, 2, width_ratios=[4.0, 1.25], wspace=0.28, left=0.07, right=0.99, top=0.93, bottom=0.12)
left = outer[0].subgridspec(3, 4, hspace=0.18, wspace=0.28)
ROWS = [('d_zygotic', 'Zygotic-gene share\n(change)'), ('d_maternal', 'Maternal-gene share\n(change)'),
        ('d_clock', 'Transcriptomic age\n(change)')]
for j, sp in enumerate(SPECIES):
    g = I[I.species == sp].reset_index(drop=True)
    labels = [f"{SHORT[a]}-{SHORT[b]}" for a, b in (s.split('->') for s in g.interval)]
    kz = int(g.d_zygotic.values.argmax())
    kl = int(np.where(g.literature_zga.values)[0][0])
    x = np.arange(len(g))
    for i, (col, ylab) in enumerate(ROWS):
        ax = fig.add_subplot(left[i, j])
        ax.axvspan(kz - 0.45, kz + 0.45, color=BAND, lw=0, zorder=0)
        ylo, yhi = (g[col + '_lo'], g[col + '_hi']) if col != 'd_clock' else (g.d_clock_lo, g.d_clock_hi)
        ax.bar(x, g[col], width=0.62, color=COLS[sp] if col == 'd_clock' else MUTED, zorder=2)
        ax.vlines(x, ylo, yhi, color=INK, lw=0.7, zorder=3)
        ax.axhline(0, color=INK2, lw=0.5, zorder=1)
        lo_all, hi_all = min(ylo.min(), 0), max(yhi.max(), 0)
        pad = 0.08 * (hi_all - lo_all)
        ax.set_ylim(lo_all - pad, hi_all + pad)
        ax.add_patch(Rectangle((kl - 0.47, lo_all - pad), 0.94, hi_all - lo_all + 2 * pad, fill=False, ec=INK2,
                               lw=0.6, ls=(0, (2, 2)), zorder=4))
        ax.set_xticks(x)
        ax.set_xticklabels(labels if i == 2 else [], rotation=60, ha='right', fontsize=5.2)
        ax.tick_params(axis='y', labelsize=5.2)
        ax.tick_params(axis='x', length=0)
        ax.grid(axis='y', color=GRID, lw=0.3)
        ax.set_axisbelow(True)
        if j == 0:
            ax.set_ylabel(ylab, fontsize=5.8, color=INK2)
        if i == 0:
            ax.set_title(sp, loc='left', fontsize=7, color=INK, pad=3)
            fs.panel(ax, 'abcd'[j], dx=-0.32 if j == 0 else -0.2, dy=1.28)

ax = fig.add_subplot(outer[1])
for sp in SPECIES:
    g = I[I.species == sp]
    rz = g.d_zygotic.rank() / len(g)
    ax.scatter(rz, g.d_clock, s=16, color=COLS[sp], edgecolor='white', lw=0.4, zorder=3, label=sp)
ax.axhline(0, color=INK2, lw=0.5)
ax.set_xlabel('Zygotic-gene gain\n(rank within species)', color=INK2)
ax.set_ylabel('Change in transcriptomic age', color=INK2)
rho = pd.Series(I.d_zygotic.values).rank().corr(pd.Series(-I.d_clock.values).rank())
ax.text(0.03, 0.03, 'pooled Spearman %.2f\npermutation p = 0.003\n(within-species, post hoc)' % rho,
        transform=ax.transAxes, fontsize=5.5, color=INK2, va='bottom')
ax.legend(frameon=False, fontsize=5.6, loc='upper right', handletextpad=0.2)
ax.grid(color=GRID, lw=0.3)
ax.set_axisbelow(True)
fs.panel(ax, 'e', dx=-0.3, dy=1.06)
fig.text(0.37, 0.015, 'Shaded = ZGA measured in these embryos; dashed = ZGA stage from the literature. Post hoc analysis.',
         ha='center', fontsize=5.8, color=INK2)
fs.save(fig, 'FigGZ_2_zga_vs_clock', outdir=f'{B}/figures')
