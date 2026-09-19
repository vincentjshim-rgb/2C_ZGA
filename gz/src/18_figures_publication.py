"""Publication Figures 1-4 (plan: results/FIGURE_PLAN_gz_v1_KR.md; Figures 3 and 4 of the plan merged into Figure 3). Reads only figures/source_data/*.tsv (written by
src/17_figdata.py). Every figure is exactly 167 mm wide; axes are placed in millimetres (no tight bounding box).
FreeSans (Helvetica metrics), 5-6.5 pt. Bootstrap intervals re-drawn here use the seeds of the scripts that produced
them (Figure 1: 20260914 and 20260917), so the plotted intervals equal the stored ones.
Colour meaning is fixed across figures (validated for colour-vision deficiency):
  red / blue  = expression up / down (also heat-map z-scores and ZGA vs clearance in the schematic)
  dark / light grey = clock contributions that push the value down / up
  black = control arm, amber = ZGA-blocking arm, teal = A485 + DUX rescue arm."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
import matplotlib
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle, Wedge
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, f'{ANALYSIS}/src')
import figstyle as fs

B = f'{ANALYSIS}/gz'
SD = f'{B}/figures/source_data'
OUT = f'{B}/figures/pub'
for f in ['FreeSans.ttf', 'FreeSansBold.ttf', 'FreeSansOblique.ttf']:
    font_manager.fontManager.addfont(f'/usr/share/fonts/truetype/freefont/{f}')
fs.use()
plt.rcParams.update({'font.sans-serif': ['FreeSans', 'DejaVu Sans'], 'font.size': 6, 'axes.labelsize': 6,
                     'xtick.labelsize': 5.6, 'ytick.labelsize': 5.6, 'legend.fontsize': 5.6, 'axes.titlesize': 6,
                     'savefig.bbox': None, 'savefig.pad_inches': 0})
MM = fs.MM
WMM = 167
INK, INK2, MUTED, GRID, BAND = '#1a1a1a', '#4a4a4a', '#8c8c8c', '#e4e4e4', '#f0e6cf'
BLUE, RED, AMBER, TEAL = '#2f6db5', '#c8453a', '#d4940f', '#1f9e89'
LRED, LBLUE = '#efc4bf', '#c9d9ee'
DOWNC, UPC = '#555555', '#b8b8b8'
ROLE = {'control': INK, 'A485': AMBER, 'Tardbp_matKO': AMBER, 'Brg1_matKO': AMBER, 'A485+Dux': TEAL}
LABEL = {'control': 'Control', 'A485': 'A485', 'A485+Dux': 'A485 + DUX', 'Tardbp_matKO': 'Maternal Tardbp KO',
         'Brg1_matKO': 'Maternal Brg1 KO'}
DS = {'GSE280522': 'P1 · GSE280522', 'GSE221985': 'P2 · GSE221985', 'GSE300734': 'P3 · GSE300734'}
CMAP = LinearSegmentedColormap.from_list('bwr_study', [BLUE, '#f7f7f7', RED])


def rd(name, **kw):
    return pd.read_csv(f'{SD}/{name}.tsv', sep='\t', **kw)


def newfig(hmm):
    fig = plt.figure(figsize=(WMM * MM, hmm * MM))
    fig._hmm = hmm
    return fig


def axmm(fig, x, y, w, h, **kw):
    """Axes placed by its top-left corner (x, y) and size in millimetres."""
    return fig.add_axes([x / WMM, 1 - (y + h) / fig._hmm, w / WMM, h / fig._hmm], **kw)


def canvas(fig, x, y, w, h):
    """Drawing area in millimetre coordinates (origin top-left)."""
    ax = axmm(fig, x, y, w, h)
    ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.axis('off')
    return ax


def letter(fig, x, y, s):
    fig.text(x / WMM, 1 - y / fig._hmm, s, fontsize=8, fontweight='bold', va='top', ha='left')


def ygrid(ax):
    ax.grid(axis='y', color=GRID, lw=0.4)
    ax.set_axisbelow(True)


def short(sym):
    return 'ENSMUSG…' + sym[-6:] if str(sym).startswith('ENSMUSG') else sym


SHORT = {'control': 'Control', 'A485': 'A485', 'A485+Dux': 'A485\n+ DUX', 'Tardbp_matKO': 'mat.\nTardbp KO',
         'Brg1_matKO': 'mat.\nBrg1 KO'}


def arm_strip(ax, df, ycol, arms, ms=11, labels=True, d_text=False):
    """Per arm: early two-cell (open) and late two-cell (filled) library points side by side; short bar = mean.
    No connecting lines."""
    for i, arm in enumerate(arms):
        c = ROLE[arm]
        means = {}
        for st, dx, filled in [('E2C', -0.2, False), ('L2C', 0.2, True)]:
            v = df[(df.arm == arm) & (df.stage == st)][ycol].values
            xs = i + dx + (np.linspace(-0.055, 0.055, len(v)) if len(v) > 1 else 0)
            ax.scatter(xs, v, s=ms, facecolor=c if filled else 'white', edgecolor=c, lw=0.8, zorder=3)
            ax.plot([i + dx - 0.12, i + dx + 0.12], [v.mean()] * 2, color=c, lw=1.3, zorder=4)
            means[st] = v.mean()
        if d_text:
            ax.text(i, 1.01, f"D {means['L2C'] - means['E2C']:+.2f}", transform=ax.get_xaxis_transform(), ha='center',
                    va='bottom', fontsize=4.9, color=INK2)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels([SHORT[a] for a in arms] if labels else [], linespacing=1.0)
    ax.set_xlim(-0.6, len(arms) - 0.4)
    ax.tick_params(axis='x', length=0)


def stage_handles():
    return [plt.Line2D([], [], marker='o', ls='', mfc='white', mec=INK, ms=3.4, label='early 2-cell (E2C)'),
            plt.Line2D([], [], marker='o', ls='', mfc=INK, mec=INK, ms=3.4, label='late 2-cell (L2C)')]


def arm_handles(pairs):
    return [Rectangle((0, 0), 1, 1, color=ROLE[a], label=l) for a, l in pairs]


def tile(c, x, y, w, h, fc, ec='white', lw=0.6, text=None, tc=INK, fs_=4.8, bold=False):
    c.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc=fc, ec=ec, lw=lw))
    if text:
        c.text(x, y, text, ha='center', va='center', fontsize=fs_, color=tc, fontweight='bold' if bold else 'normal')


def embryo_icon(c, x, y, kind, r=3.0):
    """Schematic preimplantation embryo in millimetre canvas coordinates."""
    cyto, edge = '#f1eadc', INK2
    c.add_patch(Circle((x, y), r, fc='white', ec=edge, lw=0.55))
    def cell(cx, cy, rr, fc=cyto):
        c.add_patch(Circle((cx, cy), rr, fc=fc, ec=edge, lw=0.35))
    if kind == 'oocyte':
        cell(x, y, r * 0.8)
        c.add_patch(Circle((x + r * 0.28, y - r * 0.2), r * 0.24, fc='white', ec=edge, lw=0.35))
    elif kind == 'zygote':
        cell(x, y, r * 0.8)
        for dx in (-0.24, 0.24):
            c.add_patch(Circle((x + dx * r, y), r * 0.17, fc='white', ec=edge, lw=0.35))
    elif kind in ('e2c', 'l2c'):
        for dx in (-0.42, 0.42):
            cell(x + dx * r, y, r * 0.44)
            c.add_patch(Circle((x + dx * r, y), r * 0.14, fc=RED if kind == 'l2c' else 'white', ec=edge, lw=0.3))
    elif kind == '4cell':
        for dx in (-0.36, 0.36):
            for dy in (-0.36, 0.36):
                cell(x + dx * r, y + dy * r, r * 0.34)
    elif kind == '8cell':
        for a in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            cell(x + 0.52 * r * np.cos(a), y + 0.52 * r * np.sin(a), r * 0.27)
    elif kind == '16cell':
        for a in np.linspace(0, 2 * np.pi, 11, endpoint=False):
            cell(x + 0.62 * r * np.cos(a), y + 0.62 * r * np.sin(a), r * 0.18)
        for a in np.linspace(0, 2 * np.pi, 5, endpoint=False):
            cell(x + 0.25 * r * np.cos(a), y + 0.25 * r * np.sin(a), r * 0.17)
    elif kind == 'morula':
        cell(x, y, r * 0.82, fc='#e6dcc6')
        for a in np.linspace(0, 2 * np.pi, 9, endpoint=False):
            c.add_patch(Circle((x + 0.5 * r * np.cos(a), y + 0.5 * r * np.sin(a)), r * 0.07, fc=edge, ec='none'))
    elif kind == 'blastocyst':
        c.add_patch(Circle((x, y), r * 0.86, fc=cyto, ec=edge, lw=0.35))
        c.add_patch(Circle((x, y), r * 0.64, fc='white', ec=edge, lw=0.35))
        c.add_patch(Wedge((x, y), r * 0.64, 200, 340, width=r * 0.34, fc='#d8c9a6', ec=edge, lw=0.35))


def forest(ax, rows, xlim, xlabel, zero=0.0, fmt='{:+.3f}'):
    """rows: dict(label, est, lo, hi, color, hollow, marker, dashed, gap). Values printed right of the axis."""
    y, ys = 0.0, []
    for r in rows:
        y -= r.get('gap', 0)
        c = r.get('color', INK)
        ax.plot([r['lo'], r['hi']], [y, y], color=c, lw=0.9, ls=(0, (2, 1.4)) if r.get('dashed') else '-',
                solid_capstyle='butt', zorder=3)
        ax.scatter([r['est']], [y], s=r.get('s', 15), marker=r.get('marker', 'o'),
                   facecolor='white' if r.get('hollow') else c, edgecolor=c, lw=0.8, zorder=4)
        ax.text(1.02, y, fmt.format(r['est']), transform=ax.get_yaxis_transform(), ha='left', va='center',
                fontsize=5.4, color=INK2)
        ys.append(y)
        y -= 1
    ax.axvline(zero, color=MUTED, lw=0.5, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([r['label'] for r in rows])
    ax.tick_params(axis='y', length=0, pad=2)
    ax.set_xlim(*xlim)
    ax.set_ylim(min(ys) - 0.7, 0.7)
    ax.set_xlabel(xlabel)
    ax.spines['left'].set_visible(False)


# =====================================================================================================
def figure1():
    fig = newfig(176)
    # ---- (a) schematic ------------------------------------------------------------------------------
    c = canvas(fig, 0, 1, WMM, 62)
    letter(fig, 1, 1.5, 'a')
    kinds = ['oocyte', 'zygote', 'e2c', 'l2c', '4cell', '8cell', '16cell', 'morula', 'blastocyst']
    names = ['Oocyte', 'Zygote', 'Early\n2-cell', 'Late\n2-cell', '4-cell', '8-cell', '16-cell', 'Morula', 'Blasto-\ncyst']
    xs = np.linspace(47, 121, len(kinds)); step = xs[1] - xs[0]; tw, th = step - 0.7, 3.6
    c.add_patch(FancyBboxPatch((xs[2] - 4.4, 1.2), xs[3] - xs[2] + 8.8, 17.6, boxstyle='round,pad=0,rounding_size=1.2',
                               fc='#fbe9e6', ec='none'))
    c.text((xs[2] + xs[3]) / 2, 19.8, 'no cell division', ha='center', va='top', fontsize=4.8, color=RED, style='italic')
    for x, k, n in zip(xs, kinds, names):
        embryo_icon(c, x, 6.2, k, r=3.3)
        c.text(x, 10.6, n, ha='center', va='top', fontsize=5.0, linespacing=0.95)
    c.text(5, 6.2, 'Mouse preimplantation\ndevelopment', ha='left', va='center', fontsize=5.8, fontweight='bold', linespacing=1.05)
    yz, ym = 25.0, 29.6
    c.text(5, yz, 'Zygotic genome activation', ha='left', va='center', fontsize=5.4, color=INK2)
    c.text(5, ym, 'Maternal RNA remaining', ha='left', va='center', fontsize=5.4, color=INK2)
    for i, x in enumerate(xs):
        if i in (1, 2):
            tile(c, x, yz, tw, th, LRED, text='minor', fs_=4.6)
        elif i == 3:
            tile(c, x, yz, tw, th, RED, text='major', tc='white', fs_=4.6, bold=True)
        else:
            tile(c, x, yz, tw, th, '#f3f3f3')
        alpha = [1.0, 0.85, 0.6, 0.32, 0.14, 0.05, 0, 0, 0][i]
        tile(c, x, ym, tw, th, matplotlib.colors.to_hex(np.array(matplotlib.colors.to_rgb(BLUE)) * alpha + (1 - alpha) * np.array([0.955] * 3)))
    rows = [('Timing across species', 'Gate 1, 2a · Fig. 1b, d', range(0, 7), 'GSE225056\nmouse, pig, cow, rabbit', '#9a9a9a'),
            ('Replication', 'Gate 2b · Fig. 1c', range(1, 9), 'GSE45719 · single embryos', '#9a9a9a'),
            ('Replication', 'Gate 2b · Fig. S1', range(0, 6), 'GSE66582 · bulk', '#9a9a9a'),
            ('Perturbation', 'Gates 3, 4 · Figs 2–5', (2, 3), 'six GEO series: A485 ± DUX,\nmaternal KOs, α-amanitin, Obox3, SCNT', INK)]
    for r_, (role, gate, cover, name, col) in enumerate(rows):
        y = 37 + r_ * 6.4
        c.text(5, y - 0.9, role, ha='left', va='center', fontsize=5.6, fontweight='bold')
        c.text(5, y + 1.6, gate, ha='left', va='center', fontsize=4.8, color=MUTED)
        for i, x in enumerate(xs):
            tile(c, x, y, tw, th, col if i in cover else '#f3f3f3')
        c.text(127, y, name, ha='left', va='center', fontsize=5.1, color=INK2, linespacing=1.0)
    c.text(xs[0] - step / 2, 33.2, 'Datasets (shaded = stages sampled)', ha='left', va='center', fontsize=4.8, color=MUTED)
    c.add_patch(FancyBboxPatch((127, 3), 38, 22, boxstyle='round,pad=0,rounding_size=1.2', fc='#f5f5f5', ec=GRID, lw=0.5))
    c.text(129, 5.5, 'Readouts', ha='left', va='top', fontsize=5.6, fontweight='bold')
    c.text(129, 9.3, '• multispecies transcriptomic clock\n  (Tyshkovskiy et al., 2026)\n• exact per-gene decomposition\n• SUPPA2 splicing (PSI)\n• plan frozen before every score',
           ha='left', va='top', fontsize=4.9, color=INK2, linespacing=1.15)

    # ---- (b) four species ----------------------------------------------------------------------------
    E = rd('gate1_embryo_tage', index_col=0)
    I = rd('gate1_intervals'); I = I[I.analysis == 'primary']
    COL = 'EN_Chronoage_Multispecies_Multitissue_scaleddiff'
    PANELS = [('Mouse', ['Oocyte', 'Zygote', 'Early-2-cell', 'Late-2-cell', '4-cell', '8-cell', '16-cell'],
               ['Oo', 'Zy', 'E2C', 'L2C', '4C', '8C', '16C'], ('Early-2-cell', 'Late-2-cell')),
              ('Pig', ['Oocyte', 'Zygote', '2-cell', '4-cell', 'Day2', 'Day3', 'Morula'],
               ['Oo', 'Zy', '2C', '4C', 'D2', 'D3', 'Mor'], ('2-cell', '4-cell')),
              ('Cow', ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
               ['Oo', 'Zy', '2C', '4C', '8C', '16C', 'Mor'], ('8-cell', '16-cell')),
              ('Rabbit', ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', '16-cell', 'Morula'],
               ['Oo', 'Zy', '2C', '4C', '8C', '16C', 'Mor'], ('8-cell', '16-cell'))]
    rng = np.random.default_rng(20260914)
    letter(fig, 1, 66, 'b')
    for j, (sp, order, labels, zga) in enumerate(PANELS):
        ax = axmm(fig, 16 + j * 37.5, 70, 34, 38)
        d = E[E.species == sp]
        a0, a1 = order.index(zga[0]), order.index(zga[1])
        ax.axvspan(a0, a1, color=BAND, lw=0, zorder=0)
        means, los, his = [], [], []
        for k, st in enumerate(order):
            v = d.loc[d.stage == st, COL].values
            ax.scatter(k + (rng.random(len(v)) - 0.5) * 0.35, v, s=2.2, color=MUTED, edgecolor='none', alpha=0.7, zorder=2)
            bm = np.array([rng.choice(v, len(v)).mean() for _ in range(2000)])
            means.append(v.mean()); los.append(np.percentile(bm, 2.5)); his.append(np.percentile(bm, 97.5))
        ax.vlines(range(len(order)), los, his, color=INK, lw=0.9, zorder=3)
        ax.scatter(range(len(order)), means, s=9, marker='D', color=INK, zorder=4)
        ax.axhline(0, color=GRID, lw=0.6, zorder=1)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(labels, fontsize=5.1)
        ax.tick_params(axis='x', length=0)
        ax.set_xlim(-0.6, len(order) - 0.4); ax.set_ylim(-0.75, 0.42)
        if j:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('Transcriptomic age\nrelative to oocyte')
        g = I[I.species == sp]
        r = g.iloc[int(g.delta.values.argmin())]
        ax.text(0.04, 0.04, f'{sp}\nlargest drop: {r.interval.replace("->", " → ").replace("Early-2-cell", "E2C").replace("Late-2-cell", "L2C").replace("Day3", "D3")}',
                transform=ax.transAxes, fontsize=5, color=INK2, va='bottom')
        ax.text((a0 + a1) / 2, 1.01, 'major ZGA', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=4.8, color=INK2)
        if j == 0:
            k = order.index('Late-2-cell')
            ax.text(k - 0.5, 0.3, f'Δ {r.delta:+.2f}', ha='center', va='center', fontsize=5.4, color=RED, fontweight='bold')
        ygrid(ax)

    # ---- (c) GSE45719 ------------------------------------------------------------------------------------
    rng = np.random.default_rng(20260917)
    E2 = rd('gate2b_R1_embryo_tage', index_col=0)
    R = rd('gate2a_intervals')
    S = rd('gate2a_simulation').set_index('species')
    ORDER = ['Zygote', 'Early-2C', 'Mid-2C', 'Late-2C', '4C', '8C', '16C', 'EarlyBl', 'MidBl', 'LateBl']
    LAB = ['Zy', 'E2C', 'M2C', 'L2C', '4C', '8C', '16C', 'EBl', 'MBl', 'LBl']
    letter(fig, 1, 118, 'c')
    ax = axmm(fig, 16, 124, 80, 40)
    ax.axvspan(0.8, 3.2, color=BAND, lw=0, zorder=0)
    ax.text(2.0, 1.01, 'two-cell stage', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=4.8, color=INK2)
    for col, colr, dx, lab in [('tAge_primary', INK, -0.12, 'all genes (V0)'), ('tAge_V2', BLUE, 0.12, 'maternal and zygotic genes removed (V2)')]:
        means, lo, hi = [], [], []
        for k, st in enumerate(ORDER):
            v = E2.loc[E2.stage == st, col].values
            ax.scatter(k + dx + (rng.random(len(v)) - 0.5) * 0.12, v, s=2.6, color=colr, alpha=0.35, edgecolor='none', zorder=2)
            bm = np.array([rng.choice(v, len(v)).mean() for _ in range(2000)])
            means.append(v.mean()); lo.append(np.percentile(bm, 2.5)); hi.append(np.percentile(bm, 97.5))
        x = np.arange(len(ORDER)) + dx
        ax.vlines(x, lo, hi, color=colr, lw=0.8, zorder=3)
        ax.scatter(x, means, s=9, marker='D', color=colr, zorder=4, label=lab)
    ax.axhline(0, color=GRID, lw=0.6)
    n = E2.groupby('stage').size().reindex(ORDER)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([f'{l}\n{v}' for l, v in zip(LAB, n.values)], fontsize=5.1, linespacing=1.1)
    ax.tick_params(axis='x', length=0)
    ax.text(-0.02, -0.075, 'n', transform=ax.transAxes, ha='right', va='top', fontsize=5.1, color=MUTED)
    ax.set_ylabel('Transcriptomic age\nrelative to zygote')
    ax.text(0.99, 0.98, 'GSE45719 · one point per embryo', transform=ax.transAxes, ha='right', va='top', fontsize=5.2, color=INK2)
    ax.legend(loc='lower left', fontsize=5.1)
    ygrid(ax)

    # ---- (d) composition ---------------------------------------------------------------------------------
    letter(fig, 104, 118, 'd')
    ax = axmm(fig, 118, 124, 47, 40)
    labels = ['V0', 'V1', 'V2', 'sim.']
    for j, (sp, colr, ecol) in enumerate([('Mouse', INK, '#8c8c8c'), ('Cow', '#9a9a9a', INK)]):
        z = R[(R.species == sp) & R.data_zga].set_index('variant')
        vals = [z.loc[v, ['delta', 'ci_lo', 'ci_hi']].values for v in ['V0_original', 'V1_drop_nonclock_dynamic', 'V2_drop_all_dynamic']]
        vals.append(S.loc[sp, ['delta_sim', 'ci_lo', 'ci_hi']].values)
        x = np.arange(4) + (j - 0.5) * 0.34
        d = np.array([v[0] for v in vals]); lo_ = np.array([v[1] for v in vals]); hi_ = np.array([v[2] for v in vals])
        interval = R[(R.species == sp) & R.data_zga].interval.iloc[0].replace('->', ' → ').replace('Early-2-cell', 'E2C').replace('Late-2-cell', 'L2C')
        ax.bar(x, d, width=0.3, color=colr, zorder=2, label=f'{sp}, {interval}')
        ax.vlines(x, lo_, hi_, color=ecol, lw=0.7, zorder=3)
    ax.axhline(0, color=INK2, lw=0.5)
    ax.set_xticks(range(4)); ax.set_xticklabels(labels)
    ax.tick_params(axis='x', length=0)
    ax.set_ylim(-0.4, 0.0)
    ax.set_ylabel('Clock change in the ZGA interval')
    ax.legend(loc='lower right', fontsize=5.0, handlelength=0.9)
    ax.text(0.5, -0.2, 'V0 all genes · V1 non-clock dynamic genes removed\nV2 all dynamic genes removed · sim. maternal clearance only',
            transform=ax.transAxes, ha='center', va='top', fontsize=4.8, color=MUTED)
    ygrid(ax)
    fs.save(fig, 'Figure1', outdir=OUT)


# =====================================================================================================
def figure2():
    L = rd('fig2b_clock_per_library')
    A = rd('gate3_arm_drops')
    fig = newfig(124)
    # ---- (a) schematic --------------------------------------------------------------------------------
    c = canvas(fig, 0, 1, WMM, 56)
    letter(fig, 1, 1.5, 'a')
    cols = [('oocyte', 'Oocyte'), ('zygote', 'Zygote'), ('e2c', 'Early 2-cell'), ('l2c', 'Late 2-cell')]
    xs = [52, 66, 80, 94]; tw, th = 12.6, 3.3
    for x, (k, n) in zip(xs, cols):
        embryo_icon(c, x, 5.2, k, r=3.0)
        c.text(x, 9.4, n, ha='center', va='top', fontsize=5.0)
    blocks = [('P1 · GSE280522', 'drug added to zygotes', [('DMSO', 'vehicle', 1), ('A485', 'A485', 1), ('A485+Dux', 'A485 + DUX', 1)]),
              ('P2 · GSE221985', 'maternal knockout', [('control', 'Control', 0), ('Tardbp_matKO', 'Maternal Tardbp KO', 0)]),
              ('P3 · GSE300734', 'maternal knockout', [('control', 'WT', 0), ('Brg1_matKO', 'Maternal Brg1 KO', 0)])]
    y = 16.5
    for name, how, arms in blocks:
        yc = y + 1.95 * (len(arms) - 1)
        c.text(5, yc - 1.1, name, ha='left', va='center', fontsize=5.6, fontweight='bold')
        c.text(5, yc + 1.5, how, ha='left', va='center', fontsize=4.8, color=MUTED)
        for k, (arm, lab, start) in enumerate(arms):
            yy = y + k * 3.9
            col = '#d9d9d9' if arm in ('DMSO', 'control') else ROLE[arm]
            for i, x in enumerate(xs):
                tile(c, x, yy, tw, th, col if i >= start else '#f3f3f3')
                if i >= 2:
                    c.plot(x, yy, 'o', ms=2.6, color='white' if col not in ('#d9d9d9', '#f3f3f3') else INK, mec=INK, mew=0.4)
            c.text(102, yy, lab, ha='left', va='center', fontsize=5.1)
        y += 3.9 * len(arms) + 2.6
    c.plot(46.5, 51.5, 'o', ms=2.6, color=INK)
    c.text(48.3, 51.5, 'library collected', ha='left', va='center', fontsize=4.8, color=MUTED)
    tile(c, 72, 51.5, 5, 2.4, AMBER); c.text(75.3, 51.5, 'perturbation present', ha='left', va='center', fontsize=4.8, color=MUTED)
    tile(c, 100, 51.5, 5, 2.4, '#d9d9d9'); c.text(103.3, 51.5, 'no perturbation', ha='left', va='center', fontsize=4.8, color=MUTED)
    c.add_patch(FancyBboxPatch((135, 5), 30, 31, boxstyle='round,pad=0,rounding_size=1.2', fc='#f5f5f5', ec=GRID, lw=0.5))
    c.text(137, 7.5, 'Drop per arm', fontsize=5.5, fontweight='bold', va='top')
    c.text(137, 11.2, 'D = mean(L2C) − mean(E2C)', fontsize=5.0, va='top')
    c.text(137, 16.5, 'Interaction', fontsize=5.5, fontweight='bold', va='top')
    c.text(137, 20.2, 'I = D(perturbed) − D(control)', fontsize=5.0, va='top')
    c.text(137, 25.5, 'I > 0: smaller decrease', fontsize=5.0, va='top', color=INK2)
    c.text(137, 29.5, 'R = D(A485 + DUX) − D(A485)', fontsize=5.0, va='top', color=INK2)
    c.text(135, 41, 'Libraries per group: 2–4 after\nprespecified QC (3 of 76 excluded)', fontsize=4.8, va='center', color=MUTED)

    # ---- (b) library points with D and I --------------------------------------------------------------
    letter(fig, 1, 60, 'b')
    fig.legend(handles=stage_handles() + arm_handles([('control', 'control'), ('A485', 'ZGA-blocking arm'), ('A485+Dux', 'A485 + DUX')]),
               loc='center', bbox_to_anchor=(0.55, 1 - 62 / fig._hmm), ncol=5, fontsize=5.1, handlelength=1.0)
    specs = [('GSE280522', ['control', 'A485', 'A485+Dux'], 17, 44), ('GSE221985', ['control', 'Tardbp_matKO'], 77, 30),
             ('GSE300734', ['control', 'Brg1_matKO'], 124, 30)]
    IKEY = {'A485': 'INTERACTION A485 - control', 'A485+Dux': 'INTERACTION A485+Dux - control',
            'Tardbp_matKO': 'INTERACTION Tardbp_matKO - control', 'Brg1_matKO': 'INTERACTION Brg1_matKO - control'}
    for j, (gse, arms, x0, w) in enumerate(specs):
        ax = axmm(fig, x0, 70, w, 34)
        arm_strip(ax, L[L.gse == gse], 'tAge', arms, d_text=True)
        ax.axhline(0, color=GRID, lw=0.6, zorder=0)
        ax.set_title(DS[gse], fontsize=5.6, fontweight='bold', pad=9, loc='left')
        for i, arm in enumerate(arms):
            if arm == 'control':
                continue
            r = A[(A.gse == gse) & (A.variant == 'V0') & (A.arm == IKEY[arm])].iloc[0]
            ax.text(i, -0.2, f"I {r['drop']:+.2f}\n[{r.ci_lo:.2f}, {r.ci_hi:.2f}]", transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=4.8, color=ROLE[arm], linespacing=1.1)
        if j == 0:
            ax.set_ylabel('Transcriptomic age\n(relative to control E2C)')
        ygrid(ax)
    rr = A[(A.gse == 'GSE280522') & (A.variant == 'V0') & (A.arm == 'RESCUE A485+Dux - A485')].iloc[0]
    fig.text(17 / WMM, 1 - 116 / fig._hmm,
             f"D, drop per arm (late − early two-cell mean). I, interaction against the control arm of the same dataset, with 95% bootstrap "
             f"interval (2,000 replicates).\nR = D(A485 + DUX) − D(A485) = {rr['drop']:+.3f} [{rr.ci_lo:.3f}, {rr.ci_hi:.3f}]. "
             f"Values with dynamic genes removed (V2) are in Table 2.", fontsize=4.8, color=MUTED, va='top', linespacing=1.3)
    fs.save(fig, 'Figure2', outdir=OUT)


# =====================================================================================================
def figure3():
    """Figure 3: decomposition of the two-cell decrease (a-c) and the DUX rescue arm (d-f)."""
    T = rd('posthoc_gate3_contribution_genes', index_col=0)
    T3 = rd('posthoc_gate3_contribution_genes_P3', index_col=0)
    CC = rd('fig3a_4c_cumulative_contributions')
    H = rd('fig4a_heatmap_z')
    Mc = rd('fig4a_heatmap_columns')
    Zs = rd('fig4b_zygotic_score_per_library')
    RC = rd('posthoc_rescue_contributions', index_col=0)
    fig = newfig(224)
    # ---- (a) ranked contributions and their running sum ------------------------------------------------
    letter(fig, 1, 1, 'a')
    d1 = CC[CC.series == 'P1_control'].reset_index(drop=True)
    d3 = CC[CC.series == 'P3_control'].reset_index(drop=True)
    ax0 = axmm(fig, 17, 6, 63, 9)
    ax0.bar(d1['rank'], d1.contribution, width=1.0, color=[DOWNC if v < 0 else UPC for v in d1.contribution], lw=0)
    ax0.axhline(0, color=INK2, lw=0.4)
    ax0.set_xlim(0, len(d1) + 1); ax0.set_xticks([])
    ax0.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(2))
    ax0.set_ylabel('Contribution\nper gene', fontsize=5.0)
    ax0.spines['bottom'].set_visible(False)
    ax0.text(0.99, 0.92, 'P1 GSE280522, control arm · 1,839 clock genes', transform=ax0.transAxes, ha='right', va='top',
             fontsize=4.9, color=INK2)
    ax = axmm(fig, 17, 17, 63, 37)
    m1, e1 = d1.loc[d1.cumulative.idxmin()], d1.iloc[-1]
    m3, e3 = d3.loc[d3.cumulative.idxmin()], d3.iloc[-1]
    ax.plot(d1['rank'], d1.cumulative, color=INK, lw=1.3, zorder=3, label=f'P1 GSE280522 ({m1.cumulative:.2f} → {e1.cumulative:.2f})')
    ax.plot(d3['rank'], d3.cumulative, color=MUTED, lw=1.1, ls=(0, (3, 1.5)), zorder=2, label=f'P3 GSE300734 ({m3.cumulative:.2f} → {e3.cumulative:.2f})')
    ax.axhline(0, color=GRID, lw=0.6)
    ax.annotate(f'{m1.cumulative:.2f}\nsum of the downward\ncontributions', (m1['rank'], m1.cumulative), xytext=(30, 12),
                textcoords='offset points', ha='left', va='bottom', fontsize=4.8, color=INK,
                arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED))
    ax.annotate(f'{e1.cumulative:.2f}\nreported change', (e1['rank'], e1.cumulative), xytext=(-10, 9),
                textcoords='offset points', ha='right', va='bottom', fontsize=4.8, color=INK,
                arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED))
    ax.set_xlim(0, len(d1) + 1); ax.set_ylim(-1.75, 0.3)
    ax.set_xlabel('Clock genes, ranked from the most negative to the most positive contribution')
    ax.set_ylabel('Running sum of contributions\n(= clock change, late − early two-cell)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.56, 0.99), fontsize=4.8, handlelength=1.6)
    ygrid(ax)

    # ---- (b) largest contributors, control vs A485 -----------------------------------------------------
    letter(fig, 95, 1, 'b')
    ax = axmm(fig, 110, 6, 46, 90)
    top = T.nsmallest(20, 'c_control').iloc[::-1]
    yy = np.arange(len(top))
    ax.barh(yy + 0.19, top.c_control, height=0.36, color=INK, zorder=2, label='control')
    ax.barh(yy - 0.19, top.c_perturbed, height=0.36, color=AMBER, zorder=2, label='A485')
    ax.set_yticks(yy); ax.set_yticklabels(top.symbol, fontsize=5.2, style='italic'); ax.tick_params(axis='y', length=0)
    ax.axvline(0, color=INK2, lw=0.5)
    ax.set_xlim(-0.036, 0.006); ax.set_ylim(-0.7, len(top) - 0.3)
    ax.set_xlabel('Contribution to the clock change')
    ax.spines['left'].set_visible(False)
    ax.legend(loc='lower left', fontsize=5.0)
    ax.text(0.0, 1.008, 'P1 GSE280522 · 20 largest downward contributors', transform=ax.transAxes, fontsize=5, color=INK2, va='bottom')
    for y, v in zip(yy, top.log2FC_control_L2C_vs_E2C):
        ax.text(1.02, y, f'{v:+.1f}', transform=ax.get_yaxis_transform(), fontsize=4.7, color=INK2, va='center')
    ax.text(1.02, len(top) - 0.35, 'log2FC', transform=ax.get_yaxis_transform(), fontsize=4.6, color=MUTED, va='bottom')

    # ---- (c) P1 vs P3 ------------------------------------------------------------------------------------
    letter(fig, 1, 58, 'c')
    ax = axmm(fig, 17, 63, 40, 40)
    j = T[['c_control', 'symbol']].join(T3[['c_control']], rsuffix='_P3', how='inner')
    r_all = j[['c_control', 'c_control_P3']].corr().iloc[0, 1]
    j = j[(j.c_control != 0) | (j.c_control_P3 != 0)]
    ax.scatter(j.c_control, j.c_control_P3, s=2.4, color=MUTED, alpha=0.5, edgecolor='none', zorder=2)
    lab = j.nsmallest(4, 'c_control')
    ax.scatter(lab.c_control, lab.c_control_P3, s=8, color=INK, edgecolor='white', lw=0.3, zorder=3)
    offs = {'Klf9': (4, -5), 'Neto2': (-22, 2), 'Pi4k2a': (4, -8), 'Smyd2': (4, 3)}
    for _, r in lab.iterrows():
        ax.annotate(r.symbol, (r.c_control, r.c_control_P3), xytext=offs.get(r.symbol, (4, 2)), textcoords='offset points',
                    fontsize=5.0, style='italic', arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED))
    lim = 0.038
    ax.axhline(0, color=GRID, lw=0.5, zorder=1); ax.axvline(0, color=GRID, lw=0.5, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=0.5, ls='--', zorder=1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xticks([-0.03, 0, 0.03]); ax.set_yticks([-0.03, 0, 0.03])
    ax.text(0.04, 0.97, f'Pearson r = {r_all:.2f}\n1,839 clock genes\n28 of the top 50 shared', transform=ax.transAxes, va='top', fontsize=5.0, color=INK2)
    ax.set_xlabel('Contribution, P1 (A485 study)')
    ax.set_ylabel('Contribution, P3 (Brg1 study)')
    fig.text(63 / WMM, 1 - 66 / fig._hmm,
             'Contribution of a gene =\nclock coefficient × change in its\npreprocessed expression between\nlate and early two-cell libraries.\nThe contributions of all genes\nadd up exactly to the clock\nchange of the arm.',
             fontsize=4.8, color=MUTED, va='top', linespacing=1.25)

    # ---- (d) heat map of the zygotic set and 2C/DUX markers ------------------------------------------------
    letter(fig, 1, 108, 'd')
    order, groups = [], []
    for arm in ['control', 'A485', 'A485+Dux']:
        for st in ['E2C', 'L2C']:
            runs = list(Mc[(Mc.arm == arm) & (Mc.stage == st)].run)
            groups.append((arm, st, len(runs)))
            order += runs
    ax = axmm(fig, 50, 118, 100, 56)
    M = H[order].values
    im = ax.imshow(np.clip(M, -2.5, 2.5), aspect='auto', cmap=CMAP, vmin=-2.5, vmax=2.5, interpolation='nearest')
    ax.set_yticks(range(len(H))); ax.set_yticklabels([short(s) for s in H.symbol], fontsize=4.8, style='italic')
    ax.tick_params(axis='y', length=0, pad=1.5)
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    nz = int((H.group == 'zygotic set').sum())
    ax.axhline(nz - 0.5, color='white', lw=1.8)
    pos = 0
    for i, (arm, st, n) in enumerate(groups):
        ax.text(pos + n / 2 - 0.5, -0.8, st, ha='center', va='bottom', fontsize=5.2)
        if i < len(groups) - 1:
            ax.axvline(pos + n - 0.5, color='white', lw=2.4 if st == 'L2C' else 0.9)
        pos += n
    pos = 0
    for arm in ['control', 'A485', 'A485+Dux']:
        n = int((Mc.arm == arm).sum())
        ax.add_patch(Rectangle((pos - 0.45, -3.3), n - 0.1, 1.15, color=ROLE[arm], clip_on=False))
        ax.text(pos + n / 2 - 0.5, -2.7, LABEL[arm], ha='center', va='center', fontsize=5.4, color='white', fontweight='bold', clip_on=False)
        pos += n
    ax.set_xlim(-0.5, len(order) - 0.5); ax.set_ylim(len(H) - 0.5, -0.5)
    c = canvas(fig, 0, 118, 30, 56)
    c.set_xlim(0, 30)
    for y0, y1, lab in [(0, nz - 1, f'Zygotic set\n(n = {nz}, defined\nin the control arm)'),
                        (nz, len(H) - 1, f'Two-cell and\nDUX-target genes\n(n = {len(H) - nz})')]:
        ya, yb = 56 * (y0 / len(H)) + 0.3, 56 * ((y1 + 1) / len(H)) - 0.3
        c.plot([27.5, 27.5], [ya, yb], color=INK2, lw=0.8)
        c.text(26, (ya + yb) / 2, lab, ha='right', va='center', fontsize=5.0, color=INK2, linespacing=1.15)
    cax = axmm(fig, 153, 120, 2.2, 20)
    cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=5, length=1.5); cb.outline.set_visible(False)
    cb.set_label('z-score, log2(CPM + 1)', fontsize=5.0)
    fig.text(50 / WMM, 1 - 176 / fig._hmm, 'P1 GSE280522 · one column per library · each gene z-scored across the 23 libraries',
             fontsize=4.8, color=MUTED, va='top')

    # ---- (e) zygotic score ----------------------------------------------------------------------------
    letter(fig, 1, 180, 'e')
    ax = axmm(fig, 18, 184, 44, 31)
    arm_strip(ax, Zs, 'zygotic_score', ['control', 'A485', 'A485+Dux'], ms=8)
    ax.set_ylabel('Zygotic-set score,\nmean log2(CPM + 1)')
    ax.set_ylim(0, 4.3)
    ax.text(0.0, 1.03, 'A485 + DUX − A485 (L2C): +1.67 (1.14–2.04)', transform=ax.transAxes, fontsize=4.8, color=INK2, va='bottom')
    ax.legend(handles=stage_handles(), loc='upper right', fontsize=4.6, handlelength=0.8, borderaxespad=0.2)
    ygrid(ax)

    # ---- (f) running sums of the same genes in the three arms --------------------------------------------
    letter(fig, 80, 180, 'f')
    ax = axmm(fig, 94, 184, 66, 31)
    nneg = int((RC['control'] < 0).sum())
    offs = {'control': -0.09, 'A485': 0.13, 'A485+Dux': 0.0}
    for series, arm, lab in [('P1_control', 'control', 'Control'), ('P1_A485_ctrlorder', 'A485', 'A485'),
                             ('P1_A485+Dux_ctrlorder', 'A485+Dux', 'A485 + DUX')]:
        d = CC[CC.series == series]
        ax.plot(d['rank'], d.cumulative, color=ROLE[arm], lw=1.2, label=lab, zorder=3)
        v = float(d[d['rank'] == nneg].cumulative.iloc[0]); e = float(d.cumulative.iloc[-1])
        ax.text(nneg + 22, v + 0.05, f'{v:.2f}', ha='left', va='bottom', fontsize=4.7, color=ROLE[arm])
        ax.text(len(d) + 14, e + offs[arm], f'{e:.2f}', ha='left', va='center', fontsize=4.7, color=ROLE[arm])
    ax.axvline(nneg, color=GRID, lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.axhline(0, color=GRID, lw=0.6, zorder=1)
    ax.text(nneg + 18, 0.22, f'the {nneg} genes with a negative\ncontribution in the control arm', fontsize=4.6, color=MUTED, va='top')
    ax.set_xlim(0, 1839 + 150); ax.set_ylim(-1.6, 0.3)
    ax.set_xlabel('Clock genes in control order (most negative contribution first)')
    ax.set_ylabel('Running sum of\ncontributions')
    ax.legend(loc='lower right', fontsize=4.8, handlelength=1.6)
    ygrid(ax)
    fs.save(fig, 'Figure3', outdir=OUT)



# =====================================================================================================
def figure4():
    G = rd('gate4_progression')
    F = rd('posthoc_gate4_crossfit_folds')
    PS = rd('fig5a_psi_zsa_events_P1', index_col=0)
    Mc = rd('fig4a_heatmap_columns')
    DP = rd('fig5b_dpsi_per_arm_P1')
    specs = [('GSE280522', ['control', 'A485', 'A485+Dux']), ('GSE221985', ['control', 'Tardbp_matKO']),
             ('GSE300734', ['control', 'Brg1_matKO'])]
    fig = newfig(152)
    # ---- (a) PSI heat map of the control-defined events -----------------------------------------------
    letter(fig, 1, 2, 'a')
    order, groups = [], []
    for arm in ['control', 'A485', 'A485+Dux']:
        for st in ['E2C', 'L2C']:
            runs = list(Mc[(Mc.arm == arm) & (Mc.stage == st)].run)
            groups.append((arm, st, len(runs)))
            order += runs
    X = PS[order].values.astype(float)
    Zm = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, ddof=1, keepdims=True)
    ax = axmm(fig, 26, 14, 94, 56)
    im = ax.imshow(np.clip(Zm, -2.5, 2.5), aspect='auto', cmap=CMAP, vmin=-2.5, vmax=2.5, interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    npos = int((PS.dPSI_control > 0).sum())
    ax.axhline(npos - 0.5, color='white', lw=1.5)
    pos = 0
    for i, (arm, st, n) in enumerate(groups):
        ax.text(pos + n / 2 - 0.5, -0.5 - 0.004 * len(PS), st, ha='center', va='bottom', fontsize=5.3)
        if i < len(groups) - 1:
            ax.axvline(pos + n - 0.5, color='white', lw=2.4 if st == 'L2C' else 0.9)
        pos += n
    pos = 0
    hdr = 0.075 * len(PS)
    for arm in ['control', 'A485', 'A485+Dux']:
        n = int((Mc.arm == arm).sum())
        ax.add_patch(Rectangle((pos - 0.45, -0.5 - hdr), n - 0.1, 0.42 * hdr, color=ROLE[arm], clip_on=False))
        ax.text(pos + n / 2 - 0.5, -0.5 - 0.79 * hdr, LABEL[arm], ha='center', va='center', fontsize=5.5, color='white',
                fontweight='bold', clip_on=False)
        pos += n
    ax.set_xlim(-0.5, len(order) - 0.5); ax.set_ylim(len(PS) - 0.5, -0.5)
    c = canvas(fig, 0, 14, 25, 56)
    c.set_xlim(0, 25)
    for y0, y1, lab in [(0, npos - 1, f'PSI rises in\ncontrol\n(n = {npos})'), (npos, len(PS) - 1, f'PSI falls in\ncontrol\n(n = {len(PS) - npos})')]:
        ya, yb = 56 * (y0 / len(PS)) + 0.3, 56 * ((y1 + 1) / len(PS)) - 0.3
        c.plot([23, 23], [ya, yb], color=INK2, lw=0.8)
        c.text(21.5, (ya + yb) / 2, lab, ha='right', va='center', fontsize=5.2, color=INK2, linespacing=1.15)
    cax = axmm(fig, 122, 16, 2.2, 22)
    cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=5, length=1.5); cb.outline.set_visible(False)
    cb.set_label('PSI, z-score per event', fontsize=5.1)
    fig.text(26 / WMM, 1 - 72 / fig._hmm, f'P1 GSE280522 · {len(PS)} splicing-activation events defined on the control arm · one column per library',
             fontsize=4.9, color=MUTED, va='top')

    # ---- (b) distribution of the per-event change per arm ------------------------------------------------
    letter(fig, 132, 2, 'b')
    ax = axmm(fig, 141, 14, 23, 52)
    meds = []
    for arm, lab in [('control', 'Control'), ('A485', 'A485'), ('A485+Dux', 'A485 + DUX')]:
        v = np.sort(DP[DP.arm == arm].dPSI_oriented.values)
        ax.step(v, np.arange(1, len(v) + 1) / len(v), where='post', color=ROLE[arm], lw=1.1, label=lab, zorder=3)
        meds.append((lab, float(np.median(v)), ROLE[arm]))
    ax.axvline(0, color=MUTED, lw=0.6, zorder=1); ax.axhline(0.5, color=GRID, lw=0.5, zorder=1)
    ax.set_xlim(-0.6, 1.0); ax.set_ylim(0, 1.0)
    ax.set_xticks([-0.5, 0, 0.5, 1.0])
    ax.set_xlabel('ΔPSI, late − early two-cell\n(oriented by the control change)')
    ax.set_ylabel('Fraction of events')
    bb = dict(fc='white', ec='none', pad=0.6)
    ax.text(0.05, 0.985, 'median ΔPSI', transform=ax.transAxes, fontsize=4.6, color=MUTED, va='top', bbox=bb, zorder=5)
    for k, (lab, m, col) in enumerate(meds):
        ax.text(0.05, 0.925 - 0.068 * k, f'{m:.2f}  {lab}', transform=ax.transAxes, fontsize=4.7, color=col, va='top', bbox=bb, zorder=5)
    ygrid(ax)

    # ---- (c) per-library progression ------------------------------------------------------------------
    letter(fig, 1, 80, 'c')
    ax = axmm(fig, 17, 86, 86, 44)
    fig.legend(handles=stage_handles() + arm_handles([('control', 'control'), ('A485', 'ZGA-blocking arm'), ('A485+Dux', 'A485 + DUX')]),
               loc='center', bbox_to_anchor=(60 / WMM, 1 - 141 / fig._hmm), ncol=5, fontsize=5.0, handlelength=1.0)
    x, centers = 0.0, []
    for gse, arms in specs:
        d = G[G.gse == gse]
        start = x
        for arm in arms:
            for st, mk, dx in [('E2C', 'o', -0.2), ('L2C', 'o', 0.2)]:
                v = d[(d.arm == arm) & (d.stage == st)].progression.values
                xx = x + dx + np.linspace(-0.05, 0.05, len(v)) if len(v) > 1 else np.array([x + dx])
                ax.scatter(xx, v, s=8, marker=mk, facecolor=ROLE[arm] if st == 'L2C' else 'white', edgecolor=ROLE[arm], lw=0.7, zorder=3)
                ax.plot([x + dx - 0.12, x + dx + 0.12], [v.mean()] * 2, color=ROLE[arm], lw=1.1, zorder=4)
            x += 1.0
        centers.append((start + x - 1.0) / 2)
        if gse != specs[-1][0]:
            ax.axvline(x - 0.5 + 0.35, color=GRID, lw=0.6)
        x += 0.7
    ax.axhline(0, color=GRID, lw=0.5); ax.axhline(1, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    ax.set_xticks(centers); ax.set_xticklabels([DS[g] for g, _ in specs])
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.6, x - 0.8)
    ax.set_ylabel('Splicing progression\n(control E2C = 0, L2C = 1)')
    ygrid(ax)

    # ---- (d) cross-fitted folds -----------------------------------------------------------------------
    letter(fig, 108, 80, 'd')
    ax = axmm(fig, 128, 86, 36, 42)
    rng5 = np.random.default_rng(1)
    ylabs = []
    for yq, q, col, lab in [(1, 'I_A485', AMBER, 'I: A485 − held-out\ncontrol'), (0, 'R', TEAL, 'R: A485 + DUX\n− A485')]:
        v = F[q].values
        yy = yq + (rng5.random(len(v)) - 0.5) * 0.28
        ax.scatter(v, yy, s=9, facecolor=col, edgecolor=[RED if (q == 'I_A485' and x >= 0) else 'white' for x in v],
                   lw=[0.9 if (q == 'I_A485' and x >= 0) else 0.3 for x in v], zorder=3)
        ax.plot([v.mean()] * 2, [yq - 0.25, yq + 0.25], color=INK, lw=1.4, zorder=4)
        npos_, nneg_ = int((v > 0).sum()), int((v < 0).sum())
        ylabs.append(f'{lab}\nmean {v.mean():+.2f}; {nneg_} of 16 < 0'.replace('-', '−'))
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_yticks([1, 0]); ax.set_yticklabels(ylabs)
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.5, 1.5); ax.set_xlim(-0.62, 0.62)
    ax.set_xlabel('Difference in progression, per fold (P1)')
    ax.spines['left'].set_visible(False)
    ax.text(-0.55, -0.3, 'one point per cross-fitted fold (n = 16); bar = mean\nred outline: folds in which A485 was not below control',
            transform=ax.transAxes, fontsize=4.7, color=MUTED, va='top')
    ax.grid(axis='x', color=GRID, lw=0.4); ax.set_axisbelow(True)
    fs.save(fig, 'Figure4', outdir=OUT)


SUPP = f'{B}/figures/supp'


def figureS3():
    """Supplementary Figure S3 (was Figure 2d): secondary single-stage contrasts, one point per library."""
    L = rd('fig2b_clock_per_library')
    S2 = rd('gate3_secondary_single_stage')
    fig = newfig(74)
    panels = [('GSE162345', [('control_45hpi', 'control\n45 hpi'), ('alpha-amanitin_45hpi', 'α-amanitin\n45 hpi'),
                             ('control_54hpi', 'control\n54 hpi'), ('alpha-amanitin_54hpi', 'α-amanitin\n54 hpi')], 14, 50),
              ('GSE235547', [('control_ICSI', 'ICSI'), ('SCNT_EGFP', 'SCNT\nEGFP'), ('SCNT_Obox3', 'SCNT\nObox3'),
                             ('control_siCtrl', 'siCtrl'), ('siObox3', 'siObox3')], 74, 46),
              ('GSE248499', [('control', 'IVF'), ('SCNT_control', 'SCNT'), ('SCNT_Kdm3a', 'SCNT\nKdm3a'), ('SCNT_Kdm4d', 'SCNT\nKdm4d')], 128, 37)]
    for gse, arms, x0, w in panels:
        ax = axmm(fig, x0, 12, w, 40)
        d = L[L.gse == gse]
        for i, (arm, lab) in enumerate(arms):
            v = d[d.arm == arm].tAge.values
            ctl = arm.startswith('control')
            ax.scatter(i + np.linspace(-0.08, 0.08, len(v)), v, s=11, facecolor='white' if ctl else INK, edgecolor=INK, lw=0.8, zorder=3)
            ax.plot([i - 0.18, i + 0.18], [v.mean()] * 2, color=INK, lw=1.3, zorder=4)
            r = S2[(S2.gse == gse) & (S2.arm == arm)]
            if len(r):
                ax.text(i, 1.01, f"{r['diff'].iloc[0]:+.2f}", transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                        fontsize=4.9, color=INK2)
        ax.set_xticks(range(len(arms))); ax.set_xticklabels([l for _, l in arms], fontsize=4.9, linespacing=1.0)
        ax.tick_params(axis='x', length=0)
        ax.set_xlim(-0.6, len(arms) - 0.4)
        ax.axhline(0, color=GRID, lw=0.6, zorder=0)
        ax.set_title(gse, fontsize=5.6, fontweight='bold', pad=8, loc='left')
        if gse == 'GSE162345':
            ax.set_ylabel("Transcriptomic age\n(relative to the dataset's control libraries)")
        ygrid(ax)
    fig.text(14 / WMM, 1 - 64 / fig._hmm,
             'open circles, control libraries; filled circles, perturbed libraries; bars, means. The number above each perturbed arm is its '
             'difference from the matched control\n(α-amanitin against the control of the same collection window; SCNT-EGFP and SCNT-Obox3 '
             'against ICSI; siObox3 against siCtrl; SCNT arms against IVF); 95% bootstrap intervals are given in Section 2.3.',
             fontsize=4.7, color=MUTED, va='top')
    fs.save(fig, 'FigureS3', outdir=SUPP)


def figureS4():
    """Supplementary Figure S4 (was Figure 3d): expression of the largest contributors."""
    K = rd('fig3d_key_genes_per_library')
    fig = newfig(76)
    # ---- (d) key genes ------------------------------------------------------------------------------------
    for k, g in enumerate(['Klf9', 'Neto2', 'Pi4k2a', 'Gpatch4', 'Psmb5']):
        ax = axmm(fig, 22 + k * 29.5, 10, 20, 46)
        arm_strip(ax, K[K.symbol == g], 'log2cpm', ['control', 'A485', 'A485+Dux'], ms=5.5, labels=False)
        ax.set_title(g, style='italic', fontsize=5.8, pad=2)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        if k == 0:
            ax.set_ylabel('log2(CPM + 1)')
        ygrid(ax)
    fig.legend(handles=stage_handles() + arm_handles([('control', 'control'), ('A485', 'A485'), ('A485+Dux', 'A485 + DUX')]),
               loc='center', bbox_to_anchor=(85 / WMM, 1 - 70 / fig._hmm), ncol=5, fontsize=5.0, handlelength=1.0)
    fs.save(fig, 'FigureS4', outdir=SUPP)


def figureS5():
    """Supplementary Figure S5 (was Figure 4d): per-gene contributions under A485 and A485 + DUX against control."""
    RC = rd('posthoc_rescue_contributions', index_col=0)
    fig = newfig(62)
    # ---- (d) contribution profiles ----------------------------------------------------------------------
    for k, (arm, lab) in enumerate([('A485', 'A485'), ('A485+Dux', 'A485 + DUX')]):
        ax = axmm(fig, 44 + k * 48, 10, 36, 36)
        ax.scatter(RC['control'], RC[arm], s=1.8, color=ROLE[arm], alpha=0.55, edgecolor='none')
        lim = 0.04
        ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=0.5, ls='--')
        ax.axhline(0, color=GRID, lw=0.5); ax.axvline(0, color=GRID, lw=0.5)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_title(f'{lab}\nr = {RC["control"].corr(RC[arm]):.2f}', fontsize=5.2, pad=2, color=ROLE[arm] if arm != 'A485' else INK2)
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([-0.03, 0, 0.03]))
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([-0.03, 0, 0.03]))
        ax.tick_params(labelsize=4.8)
        if k:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('Perturbed arm')
    fig.text(86 / WMM, 1 - 55 / fig._hmm, 'Contribution, control arm', ha='center', va='top', fontsize=6)
    fs.save(fig, 'FigureS5', outdir=SUPP)


def figureS6():
    """Supplementary Figure S6 (was Figure 5a): PCA display of PSI."""
    PCA = rd('fig5a_psi_pca')
    specs = [('GSE280522', ['control', 'A485', 'A485+Dux']), ('GSE221985', ['control', 'Tardbp_matKO']),
             ('GSE300734', ['control', 'Brg1_matKO'])]
    fig = newfig(56)
    # ---- (a) PSI PCA -------------------------------------------------------------------------------------
    for j, (gse, arms) in enumerate(specs):
        ax = axmm(fig, 8 + j * 33, 12, 27, 30)
        d = PCA[PCA.gse == gse]
        for arm in arms:
            for st, mk in [('E2C', 'o'), ('L2C', 'o')]:
                q = d[(d.arm == arm) & (d.stage == st)]
                ax.scatter(q.PC1, q.PC2, s=12, marker=mk, facecolor=ROLE[arm] if st == 'L2C' else 'white',
                           edgecolor=ROLE[arm], lw=0.8, zorder=3)
        ax.set_xlabel(f'PC1 ({100 * d.var_PC1.iloc[0]:.0f}%)', labelpad=1)
        ax.set_ylabel(f'PC2 ({100 * d.var_PC2.iloc[0]:.0f}%)', labelpad=1)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'{DS[gse]}\n{int(d.n_events.iloc[0]):,} events', fontsize=5.2, pad=2, color=INK2)
    handles = [plt.Line2D([], [], marker='o', ls='', mfc='white', mec=INK, ms=3.6, label='early 2-cell'),
               plt.Line2D([], [], marker='o', ls='', mfc=INK, mec=INK, ms=3.6, label='late 2-cell'),
               Rectangle((0, 0), 1, 1, color=INK, label='control'),
               Rectangle((0, 0), 1, 1, color=AMBER, label='ZGA-blocking arm'),
               Rectangle((0, 0), 1, 1, color=TEAL, label='A485 + DUX')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(107 / WMM, 1 - 13 / fig._hmm), fontsize=5.2, labelspacing=0.5)
    fig.text(107 / WMM, 1 - 40 / fig._hmm, 'PCA of per-event PSI on the Gate 4\nfiltered events; display only',
             fontsize=4.8, color=MUTED, va='top')

    fs.save(fig, 'FigureS6', outdir=SUPP)


def figureS1():
    """Supplementary Figure S1: GSE66582 library clock values (Gate 2b R2; direction only)."""
    D = rd('gate2b_R2_library_tage')
    ORDER = ['MII', 'Zygote', 'Early-2C', '2C', '4C', '8C', 'ICM']
    LAB = ['MII', 'Zy', 'E2C', '2C', '4C', '8C', 'ICM']
    fig = newfig(62)
    ax = axmm(fig, 22, 8, 92, 40)
    ax.axvspan(1.6, 3.4, color=BAND, lw=0, zorder=0)
    ax.text(2.5, 1.01, 'two-cell stage', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=4.8, color=INK2)
    for col, colr, dx, lab in [('tAge_primary', INK, -0.15, 'all genes (V0)'), ('tAge_V2', BLUE, 0.15, 'maternal and zygotic genes removed (V2)')]:
        for k, st in enumerate(ORDER):
            v = D.loc[D.stage == st, col].values
            ax.scatter(k + dx + np.linspace(-0.04, 0.04, len(v)), v, s=9, color=colr, edgecolor='white', lw=0.3, zorder=3,
                       label=lab if k == 0 else None)
            ax.plot([k + dx - 0.1, k + dx + 0.1], [v.mean()] * 2, color=colr, lw=1.2, zorder=4)
    ax.axhline(0, color=GRID, lw=0.6)
    n = D.groupby('stage').size().reindex(ORDER)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([f'{l}\n{v}' for l, v in zip(LAB, n.values)], fontsize=5.1, linespacing=1.1)
    ax.tick_params(axis='x', length=0)
    ax.text(-0.02, -0.075, 'n', transform=ax.transAxes, ha='right', va='top', fontsize=5.1, color=MUTED)
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.set_ylabel('Transcriptomic age\nrelative to MII oocyte')
    ygrid(ax)
    fig.legend(handles=[plt.Line2D([], [], marker='o', ls='', mfc=c, mec=c, ms=3.4, label=l)
                        for c, l in [(INK, 'all genes (V0)'), (BLUE, 'maternal and zygotic\ngenes removed (V2)')]],
               loc='upper left', bbox_to_anchor=(120 / WMM, 1 - 10 / fig._hmm), fontsize=5.1, labelspacing=0.8)
    fig.text(120 / WMM, 1 - 30 / fig._hmm, 'GSE66582 · one point per library;\nbar = mean; E2C → 2C change\n'
             'V0 −0.26, V2 −0.18 (direction only)', fontsize=4.9, color=INK2, va='top')
    fs.save(fig, 'FigureS1', outdir=SUPP)


def figureS2():
    """Supplementary Figure S2: post-hoc zygotic timing (was FigGZ_2), redrawn without bars."""
    I = rd('gate1alt_intervals')
    SH = {'Oocyte': 'Oo', 'Zygote': 'Zy', 'Early-2-cell': 'E2C', 'Late-2-cell': 'L2C', '2-cell': '2C', '4-cell': '4C',
          '8-cell': '8C', '16-cell': '16C', 'Morula': 'Mor', 'Day2': 'D2', 'Day3': 'D3'}
    SPECIES = ['Mouse', 'Pig', 'Cow', 'Rabbit']
    ROWS = [('d_zygotic', 'Zygotic-gene share\n(change)'), ('d_maternal', 'Maternal-gene share\n(change)'),
            ('d_clock', 'Transcriptomic age\n(change)')]
    fig = newfig(118)
    for j, sp in enumerate(SPECIES):
        g = I[I.species == sp].reset_index(drop=True)
        labels = [f'{SH[a]}–{SH[b]}' for a, b in (s.split('->') for s in g.interval)]
        kz = int(g.d_zygotic.values.argmax())
        kl = int(np.where(g.literature_zga.values)[0][0])
        x = np.arange(len(g))
        letter(fig, 12 + j * 29 if j else 1, 2, 'abcd'[j])
        for i, (col, ylab) in enumerate(ROWS):
            ax = axmm(fig, 21 + j * 29, 9 + i * 30, 22, 24)
            ax.axvspan(kz - 0.45, kz + 0.45, color=BAND, lw=0, zorder=0)
            ax.axvspan(kl - 0.47, kl + 0.47, fill=False, ec=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)
            ax.vlines(x, g[col + '_lo'], g[col + '_hi'], color=INK, lw=0.8, zorder=3)
            ax.scatter(x, g[col], s=8, marker='D', color=INK, zorder=4)
            ax.axhline(0, color=MUTED, lw=0.5, zorder=1)
            ax.set_xticks(x)
            ax.set_xticklabels(labels if i == 2 else [], rotation=60, ha='right', rotation_mode='anchor', fontsize=4.8)
            ax.tick_params(axis='x', length=0)
            ax.tick_params(axis='y', labelsize=4.8)
            ax.set_xlim(-0.6, len(g) - 0.4)
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
            ygrid(ax)
            if j == 0:
                ax.set_ylabel(ylab, fontsize=5.4)
            if i == 0:
                ax.set_title(sp, fontsize=6, pad=3)
    # ---- (e) pooled ---------------------------------------------------------------------------------------
    letter(fig, 136, 2, 'e')
    ax = axmm(fig, 143, 9, 22, 54)
    for sp, mk in zip(SPECIES, ['o', 's', '^', 'D']):
        g = I[I.species == sp]
        ax.scatter(g.d_zygotic.rank() / len(g), g.d_clock, s=9, marker=mk, facecolor='white' if sp in ('Pig', 'Rabbit') else INK,
                   edgecolor=INK, lw=0.7, zorder=3, label=sp)
    ax.axhline(0, color=MUTED, lw=0.5)
    ax.set_xlabel('Zygotic-gene gain\n(rank within species)')
    ax.set_ylabel('Change in transcriptomic age')
    ax.set_xticks([0.25, 0.5, 0.75, 1.0]); ax.set_xticklabels(['.25', '.50', '.75', '1'])
    ax.tick_params(labelsize=4.8)
    ygrid(ax)
    ax.legend(loc='upper left', bbox_to_anchor=(-0.55, -0.28), fontsize=5.0, ncol=2, handletextpad=0.2, columnspacing=0.8)
    ax.text(-0.55, -0.62, 'pooled Spearman ρ = 0.52\nwithin-species permutation\np = 0.0034 (post hoc)', transform=ax.transAxes,
            fontsize=4.9, color=INK2, va='top')
    fig.text(21 / WMM, 1 - 112 / fig._hmm, 'shaded: largest measured zygotic gain in these embryos · dashed outline: ZGA interval '
             'from the literature · diamonds and bars: change and 95% bootstrap interval', fontsize=4.8, color=MUTED, va='center')
    fs.save(fig, 'FigureS2', outdir=SUPP)


def figureS7():
    """Supplementary Figure S7: library quality control of the 76 Gate 3 libraries (prespecified rule)."""
    Q = rd('figS7_library_qc')
    ORDER = [('GSE280522', 'P1'), ('GSE221985', 'P2'), ('GSE300734', 'P3'),
             ('GSE162345', 'secondary'), ('GSE248499', 'secondary'), ('GSE235547', 'secondary')]
    fig = newfig(96)
    for k, (gse, role) in enumerate(ORDER):
        d = Q[Q.gse == gse]
        ax = axmm(fig, 18 + (k % 3) * 52, 8 + (k // 3) * 42, 40, 28)
        ok, bad = d[d.qc_pass], d[~d.qc_pass]
        ax.scatter(ok.p_pseudoaligned, ok.detected_genes / 1000, s=8, color=INK, edgecolor='white', lw=0.3, zorder=3)
        ax.scatter(bad.p_pseudoaligned, bad.detected_genes / 1000, s=16, marker='x', color=INK, lw=0.9, zorder=4)
        ax.axvline(30, color=MUTED, lw=0.6, ls=(0, (2, 2)))
        ax.axhline(10 ** d.detected_cut_log10.iloc[0] / 1000, color=MUTED, lw=0.6, ls=(0, (2, 2)))
        ax.set_xlim(0, 100)
        lo, hi = d.detected_genes.min() / 1000, d.detected_genes.max() / 1000
        cut = 10 ** d.detected_cut_log10.iloc[0] / 1000
        ax.set_ylim(min(lo, cut) - 0.8, hi + 0.8)
        ax.set_title(f'{role} · {gse}\n{int(d.qc_pass.sum())} of {len(d)} passed', fontsize=5.2, pad=2, color=INK2)
        ax.tick_params(labelsize=4.8)
        ygrid(ax)
        if k % 3 == 0:
            ax.set_ylabel('Detected genes (×1,000)')
        if k // 3 == 1:
            ax.set_xlabel('Reads pseudoaligned (%)')
    fig.legend(handles=[plt.Line2D([], [], marker='o', ls='', mfc=INK, mec=INK, ms=3.2, label='passed'),
                        plt.Line2D([], [], marker='x', ls='', mec=INK, ms=3.6, mew=0.9, label='excluded'),
                        plt.Line2D([], [], ls=(0, (2, 2)), color=MUTED, lw=0.6, label='cut-offs: 30% pseudoaligned; '
                                   'detected genes ≥ dataset median − 3 MAD (log10)')],
               loc='center', bbox_to_anchor=(0.5, 1 - 91 / fig._hmm), ncol=3, fontsize=5.0, handlelength=1.6)
    fs.save(fig, 'FigureS7', outdir=SUPP)


if __name__ == '__main__':
    for w in (sys.argv[1:] or ['1', '2', '3', '4', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7']):
        {'1': figure1, '2': figure2, '3': figure3, '4': figure4, 'S1': figureS1, 'S2': figureS2,
         'S3': figureS3, 'S4': figureS4, 'S5': figureS5, 'S6': figureS6, 'S7': figureS7}[w]()
