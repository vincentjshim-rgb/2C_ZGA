"""Diagnostic figure for exploration A: net and gross side by side for every interval and arm.

Not a publication figure — a reading aid for the direction decision. Output: results/explore_A_netgross.png
"""
import importlib.util
import os

import numpy as np
import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
spec = importlib.util.spec_from_file_location('figmod', f'{B}/src/18_figures_publication.py')
fm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fm)
INK, MUTED, GRID, BAND, AMBER, TEAL, BLUE = fm.INK, fm.MUTED, fm.GRID, fm.BAND, fm.AMBER, fm.TEAL, fm.BLUE

I = pd.read_csv(f'{B}/results/explore_A_netgross_intervals.tsv', sep='\t')
A = pd.read_csv(f'{B}/results/explore_A_netgross_arms.tsv', sep='\t')
SHORT = {'Oocyte': 'Oo', 'Zygote': 'Zy', 'Early-2-cell': 'E2C', 'Late-2-cell': 'L2C', '2-cell': '2C', '4-cell': '4C',
         '8-cell': '8C', '16-cell': '16C', 'Morula': 'Mor', 'Day2': 'D2', 'Day3': 'D3', 'Late blastocyst': 'LBl',
         'Early-2C': 'E2C', 'Mid-2C': 'M2C', 'Late-2C': 'L2C', 'EarlyBl': 'EBl', 'MidBl': 'MBl', 'LateBl': 'LBl',
         'MII': 'MII', '2C': '2C', '4C': '4C', '8C': '8C', '16C': '16C'}


def lab(iv):
    a, c = iv.split('->')
    return f'{SHORT.get(a, a)}→{SHORT.get(c, c)}'


panels = [('GSE225056', 'Mouse'), ('GSE225056', 'Cow'), ('GSE225056', 'Pig'), ('GSE225056', 'Rabbit'),
          ('GSE45719', 'Mouse'), ('GSE66582', 'Mouse'), ('GSE36552', 'Human')]
fig = fm.newfig(150)
for k, (ds, sp) in enumerate(panels):
    t = I[(I.dataset == ds) & (I.species == sp)].reset_index(drop=True)
    x = np.arange(len(t))
    col, row = k % 4, k // 4
    ax = fm.axmm(fig, 12 + col * 40, 8 + row * 48, 30, 30)
    for i, r in t.iterrows():
        if r.is_zga:
            ax.axvspan(i - 0.5, i + 0.5, color=BAND, lw=0, zorder=0)
    ax.bar(x - 0.18, t.G, width=0.36, color=MUTED, label='gross G', zorder=2)
    ax.bar(x + 0.18, t.D * 5, width=0.36, color=[BLUE if v < 0 else AMBER for v in t.D], label='net D × 5', zorder=2)
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([lab(v) for v in t.interval], rotation=90, fontsize=5)
    ax.set_title(f'{sp} · {ds}', fontsize=6, loc='left', pad=2)
    gmax = t.loc[t.G.idxmax()]
    ax.text(0.02, 0.97, f'gross max {lab(gmax.interval)} (P {gmax.P_max_G:.2f})\nnet min {lab(t.loc[t.D.idxmin(), "interval"])}',
            transform=ax.transAxes, fontsize=5, va='top', color=INK)
    if k == 0:
        ax.legend(fontsize=5, loc='upper right', frameon=False)
    fm.ygrid(ax)

# arms
ax = fm.axmm(fig, 12 + 3 * 40, 8 + 48, 30, 30)
t = A.reset_index(drop=True)
x = np.arange(len(t))
ax.bar(x - 0.18, t.G, width=0.36, color=MUTED, zorder=2)
ax.bar(x + 0.18, t.D * 5, width=0.36, color=[BLUE if v < 0 else AMBER for v in t.D], zorder=2)
ax.axhline(0, color=INK, lw=0.5)
ax.set_xticks(x)
ax.set_xticklabels([f'{r.dataset[-6:]}\n{r.arm}' for _, r in t.iterrows()], rotation=90, fontsize=4.5)
ax.set_title('Perturbation arms (E2C→L2C)', fontsize=6, loc='left', pad=2)
fm.ygrid(ax)

fig.text(12 / fm.WMM, 1 - 112 / fig._hmm,
         'Exploration A. Grey: gross G = Σ|cᵢΔᵢ|; coloured: net D (× 5 for visibility; blue = decrease). Shaded: literature major-ZGA '
         'interval of the frozen plan.\nP1 (gross max = literature ZGA stage): 2 of 4 — unsupported by the rule; gross max = measured '
         'zygotic-gain max in 3 of 4. P2: gross lower in all 3 perturbed arms; rescue +0.45 [0.13, 0.78].\nP3 human: gross max at '
         '4C→8C as predicted but P_max 0.48. P4: gross max within the two-cell stage in both mouse replication series.',
         fontsize=5.5, color=MUTED, va='top')
fig.savefig(f'{B}/results/explore_A_netgross.png', dpi=300)
print('wrote results/explore_A_netgross.png')
