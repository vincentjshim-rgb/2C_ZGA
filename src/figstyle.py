"""Publication figure style for a reproductive-ageing / bioinformatics journal.

Conventions followed (Aging Cell, npj Aging, Biology of Reproduction house style):
  - Arial, 7 pt base text, 8 pt bold panel letters
  - single column 85 mm, 1.5 column 114 mm, double column 174 mm
  - left and bottom spines only, 0.5 pt, ticks outward
  - no gridlines, no panel titles, no legend frames
  - every data point shown when n is small enough to show
  - vector PDF plus 600 dpi TIFF/PNG
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

MM = 1 / 25.4
W1, W15, W2 = 85 * MM, 114 * MM, 174 * MM       # journal column widths

# restrained, colour-blind-safe, print-safe
YOUNG = '#4E79A7'        # blue
MATURE = '#C0504D'       # brick red
GREY = '#7F7F7F'
LGREY = '#BFBFBF'
ACCENT = '#4E7A4E'       # muted green, for a third series only
SIRE_PAL = ['#4E79A7', '#6B93BC', '#8AAED0', '#A9C8E3',
            '#C0504D', '#CE706C', '#DC908C', '#E9B0AD']

_RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5,
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': False,
    'axes.labelpad': 2.0,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.major.size': 2.0,
    'ytick.major.size': 2.0,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.size': 1.2,
    'ytick.minor.size': 1.2,
    'lines.linewidth': 0.8,
    'lines.markersize': 3,
    'legend.frameon': False,
    'legend.handlelength': 1.2,
    'legend.handletextpad': 0.5,
    'legend.columnspacing': 1.0,
    'legend.borderpad': 0.0,
    'figure.dpi': 200,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.01,
    'pdf.fonttype': 42,          # editable text in Illustrator
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
}


def use():
    plt.rcParams.update(_RC)


def panel(ax, letter, dx=-0.19, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight='bold', va='top', ha='left')


def stars(p):
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 0.05: return '*'
    return 'n.s.'


def sig_bracket(ax, x1, x2, y, text, h=None, lw=0.5, fs=6.5, color='k'):
    """Standard significance bracket between two x positions."""
    if h is None:
        h = 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, c=color,
            solid_capstyle='butt', clip_on=False)
    ax.text((x1 + x2) / 2, y + h * 1.15, text, ha='center', va='bottom',
            fontsize=fs, color=color, clip_on=False)


def beeswarm(ax, x, values, width=0.28, rng=None, **kw):
    """Deterministic jitter that spreads points by local density."""
    v = np.asarray(values, float)
    if len(v) == 0:
        return
    order = np.argsort(v)
    nb = max(1, int(np.ceil(len(v) / 6)))
    edges = np.linspace(v.min(), v.max() + 1e-12, nb + 1)
    off = np.zeros(len(v))
    for i in range(nb):
        sel = np.where((v >= edges[i]) & (v < edges[i + 1]))[0]
        k = len(sel)
        if k == 0:
            continue
        sel = sel[np.argsort(v[sel])]
        pos = (np.arange(k) - (k - 1) / 2) / max(k, 2)
        off[sel] = pos * 2 * width
    ax.scatter(x + off, v, **kw)


def boxstats(ax, x, values, width=0.5, color='k', lw=0.6, median_lw=1.0):
    """Thin box + whisker drawn in the traditional (non-filled) style."""
    v = np.asarray(values, float)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    iqr = q3 - q1
    lo = v[v >= q1 - 1.5 * iqr].min()
    hi = v[v <= q3 + 1.5 * iqr].max()
    ax.plot([x, x], [lo, q1], c=color, lw=lw, zorder=3)
    ax.plot([x, x], [q3, hi], c=color, lw=lw, zorder=3)
    for yy in (lo, hi):
        ax.plot([x - width / 4, x + width / 4], [yy, yy], c=color, lw=lw, zorder=3)
    ax.add_patch(plt.Rectangle((x - width / 2, q1), width, q3 - q1,
                               fill=False, ec=color, lw=lw, zorder=3))
    ax.plot([x - width / 2, x + width / 2], [med, med], c=color, lw=median_lw, zorder=4)


def meanci(ax, x, m, lo, hi, color='k', ms=4, lw=0.8, marker='o', zorder=5):
    """Point estimate with 95% CI, forest-plot style."""
    ax.plot([lo, hi], [x, x], c=color, lw=lw, solid_capstyle='butt', zorder=zorder)
    ax.plot([lo, lo], [x - 0.12, x + 0.12], c=color, lw=lw, zorder=zorder)
    ax.plot([hi, hi], [x - 0.12, x + 0.12], c=color, lw=lw, zorder=zorder)
    ax.plot([m], [x], marker=marker, c=color, ms=ms, mec='none', zorder=zorder + 1)


def save(fig, stem, outdir='figures'):
    import os
    os.makedirs(outdir, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/{stem}.{ext}')
    plt.close(fig)
    print(f'  wrote {outdir}/{stem}.pdf and .png')
