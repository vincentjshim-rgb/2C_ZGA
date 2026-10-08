"""One figure for the feedback note: the running sums behind the argument, and the human-versus-mouse comparison.

Four panels from the same source data as Figures 3 and 5 (nothing recomputed):
  (a) mouse control arm, ranked contributions accumulated: −1.40 → −0.24
  (b) the same genes in control order under A485 and A485 + DUX
  (c) human 8-cell → morula: −1.24 → −0.10
  (d) per-gene contribution, human against mouse, for the genes non-zero in both

Output: results/for_feedback_VNG/send/story_figure.{pdf,png}. Run with the gz environment python.
"""
import importlib.util
import os

import numpy as np

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
OUT = f'{B}/results/for_feedback_VNG/send'
os.makedirs(OUT, exist_ok=True)

spec = importlib.util.spec_from_file_location('figmod', f'{B}/src/18_figures_publication.py')
fm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fm)
INK, INK2, MUTED, GRID, ROLE = fm.INK, fm.INK2, fm.MUTED, fm.GRID, fm.ROLE

CC = fm.rd('fig3a_4c_cumulative_contributions')
RC = fm.rd('posthoc_rescue_contributions', index_col=0)
HC = fm.rd('figHb_human_cumulative')
J = fm.rd('figHc_human_vs_mouse', index_col=0)

fig = fm.newfig(118)
W, H = 62, 36                      # panel size in mm
X = {'l': 15, 'r': 100}
Y = {'t': 9, 'b': 60}


def running_sum(ax, d, title, ylim):
    """A ranked running sum with the two numbers that matter written on it."""
    m, e = d.loc[d.cumulative.idxmin()], d.iloc[-1]
    ax.plot(d['rank'], d.cumulative, color=INK, lw=1.3, zorder=3)
    ax.axhline(0, color=GRID, lw=0.6)
    ax.annotate(f'{m.cumulative:.2f}\nsum of the downward\ncontributions', (m['rank'], m.cumulative),
                xytext=(-34, 42), textcoords='offset points', ha='left', va='bottom', fontsize=6, color=INK,
                arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED))
    ax.annotate(f'{e.cumulative:.2f}\nreported change', (e['rank'], e.cumulative), xytext=(-8, 10),
                textcoords='offset points', ha='right', va='bottom', fontsize=6, color=INK,
                arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED))
    ax.set_xlim(0, len(d) * 1.04); ax.set_ylim(*ylim)
    ax.set_xlabel('Clock genes, ranked by contribution')
    ax.set_ylabel('Running sum of\ncontributions')
    ax.set_title(title, fontsize=6, pad=3, loc='left', color=INK2)
    fm.ygrid(ax)
    return m, e


# ---- (a) mouse control
fm.letter(fig, 1, 2, 'a')
ax = fm.axmm(fig, X['l'], Y['t'], W, H)
d1 = CC[CC.series == 'P1_control'].reset_index(drop=True)
running_sum(ax, d1, 'Mouse, early → late two-cell · GSE280522 control arm', (-1.6, 0.3))

# ---- (b) the same genes in the three arms
fm.letter(fig, 86, 2, 'b')
ax = fm.axmm(fig, X['r'], Y['t'], W, H)
nneg = int((RC['control'] < 0).sum())
offs = {'control': -0.09, 'A485': 0.13, 'A485+Dux': 0.0}
for series, arm, lab in [('P1_control', 'control', 'Control'), ('P1_A485_ctrlorder', 'A485', 'A485'),
                         ('P1_A485+Dux_ctrlorder', 'A485+Dux', 'A485 + DUX')]:
    d = CC[CC.series == series]
    ax.plot(d['rank'], d.cumulative, color=ROLE[arm], lw=1.2, label=lab, zorder=3)
    v = float(d[d['rank'] == nneg].cumulative.iloc[0]); e = float(d.cumulative.iloc[-1])
    ax.text(nneg + 22, v + 0.05, f'{v:.2f}', ha='left', va='bottom', fontsize=6, color=ROLE[arm])
    ax.text(len(d) + 14, e + offs[arm], f'{e:.2f}', ha='left', va='center', fontsize=6, color=ROLE[arm])
ax.axvline(nneg, color=GRID, lw=0.6, ls=(0, (2, 2)), zorder=1)
ax.axhline(0, color=GRID, lw=0.6, zorder=1)
ax.text(nneg + 18, 0.22, f'the {nneg} genes negative\nin the control arm', fontsize=6, color=MUTED, va='top')
ax.set_xlim(0, 1839 + 150); ax.set_ylim(-1.6, 0.3)
ax.set_xlabel('Clock genes in control order')
ax.set_ylabel('Running sum of\ncontributions')
ax.set_title('The same genes under A485 and A485 + DUX', fontsize=6, pad=12, loc='left', color=INK2)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=6, handlelength=1.2, frameon=False,
          columnspacing=1.4, borderaxespad=0, handletextpad=0.5)
fm.ygrid(ax)

# ---- (c) human
fm.letter(fig, 1, 53, 'c')
ax = fm.axmm(fig, X['l'], Y['b'], W, H)
running_sum(ax, HC, 'Human, 8-cell → morula · GSE36552, 20 embryo pseudobulks', (-1.45, 0.35))

# ---- (d) human against mouse
fm.letter(fig, 86, 53, 'd')
ax = fm.axmm(fig, X['r'] + 14, Y['b'], H, H)
r = float(np.corrcoef(J.human_8cell_to_morula, J.mouse_E2C_to_L2C)[0, 1])
ax.scatter(J.mouse_E2C_to_L2C, J.human_8cell_to_morula, s=1.6, color=MUTED, alpha=0.5, lw=0, rasterized=True)
lim = 0.035
ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=0.5, ls='--')
ax.axhline(0, color=GRID, lw=0.5); ax.axvline(0, color=GRID, lw=0.5)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xticks([-0.03, 0, 0.03]); ax.set_yticks([-0.03, 0, 0.03])
ax.set_xlabel('Mouse, early → late two-cell')
ax.set_ylabel('Human, 8-cell → morula')
ax.set_title(f'Human against mouse, per gene · r = {r:.2f}, n = {len(J):,}', fontsize=6, pad=3,
             loc='left', color=INK2)

fig.text(X['l'] / fm.WMM, 1 - 107 / fig._hmm,
         'Contribution = clock coefficient × change in the preprocessed feature between the two stages; the terms sum '
         'to the reported clock difference (an identity, to 1e-9).\nMouse: 23 re-quantified libraries, GSE280522 '
         '(Xiao et al., 2025); the composition agrees with an independent study, GSE300734, at r = 0.74 (not shown).\n'
         'Human: GSE36552 (Yan et al., 2013), scored against the oocytes of the same dataset.\n'
         'Multi-species chronological tAge clock (Tyshkovskiy et al., 2026), applied as published.',
         fontsize=6, color=MUTED, va='top')
fm.save(fig, 'story_figure', OUT)
