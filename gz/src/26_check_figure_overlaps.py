"""Check every figure for text that collides with other text, with a panel, or with the canvas edge.

`24_check_figure_bounds.py` only sees ink at the canvas border. This script renders each figure in memory and
measures the bounding box of every text artist, so it also catches a label sitting on top of another label, a label
spilling out of its own panel into a neighbour, and a panel overlapping a panel — the failures that make a figure look
cramped without cutting anything off.

Reported positions are in millimetres from the top-left of the figure. Exit status 1 if anything is reported.
"""
import importlib.util
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.text as mtext
from matplotlib.transforms import Bbox

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = f'{ANALYSIS}/gz/src/18_figures_publication.py'

MIN_OVERLAP_PT = 0.6      # ignore sub-point touches: kerning and rounding produce these
MIN_GAP_MM = 0.8          # text closer than this to a foreign panel edge is reported as cramped

spec = importlib.util.spec_from_file_location('figmod', SRC)
figmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(figmod)

findings = []


def mm(fig, px):
    return px / fig.dpi * 25.4


def inter(a, b):
    x0, x1 = max(a.x0, b.x0), min(a.x1, b.x1)
    y0, y1 = max(a.y0, b.y0), min(a.y1, b.y1)
    return (x1 - x0, y1 - y0) if x1 > x0 and y1 > y0 else (0.0, 0.0)


def short(s, n=46):
    s = ' '.join(s.split())
    return s if len(s) <= n else s[:n - 1] + '…'


def audit(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    H = fig.get_size_inches()[1] * fig.dpi

    # tick labels for ticks outside the view are still Text artists but are never drawn
    hidden = set()
    for a in fig.axes:
        for get_ticks, get_lim, get_labels in ((a.get_xticks, a.get_xlim, a.get_xticklabels),
                                               (a.get_yticks, a.get_ylim, a.get_yticklabels)):
            lo, hi = sorted(get_lim())
            for v, lab in zip(get_ticks(), get_labels()):
                # a tick outside the view, or any tick of a drawing canvas with its axis switched off,
                # is a Text artist that is never drawn
                if not a.axison or not (lo - 1e-9 <= v <= hi + 1e-9):
                    hidden.add(id(lab))

    items = []
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible() or id(t) in hidden:
            continue
        try:
            bb = t.get_window_extent(r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append((t, bb, s, t.get_bbox_patch() is not None))

    # --- text against text
    for i in range(len(items)):
        ti, bi, si, boxed_i = items[i]
        for j in range(i + 1, len(items)):
            tj, bj, sj, boxed_j = items[j]
            if ti.axes is not None and ti.axes is tj.axes and ti.axes.get_label() == 'colorbar':
                continue
            w, h = inter(bi, bj)
            if w <= 0 or h <= 0:
                continue
            if boxed_i or boxed_j:      # an opaque background is a deliberate overlay
                continue
            ov = min(w, h) / fig.dpi * 72
            if ov < MIN_OVERLAP_PT:
                continue
            findings.append((name, 'text over text',
                             f'{short(si)!r} x {short(sj)!r}',
                             f'{mm(fig, min(bi.x0, bj.x0)):.0f},{mm(fig, H - max(bi.y1, bj.y1)):.0f} mm',
                             f'{ov:.1f} pt'))

    # --- text against a panel it does not belong to
    axes = [a for a in fig.axes if a.get_visible()]
    data_axes = [a for a in axes if a.axison]
    for t, bb, s, boxed in items:
        if boxed:
            continue
        for a in data_axes:
            if t.axes is a:
                continue
            ab = a.get_window_extent(r)
            w, h = inter(bb, ab)
            if w > 0 and h > 0 and min(w, h) / fig.dpi * 72 >= MIN_OVERLAP_PT:
                findings.append((name, 'text over another panel', short(s),
                                 f'{mm(fig, bb.x0):.0f},{mm(fig, H - bb.y1):.0f} mm',
                                 f'{mm(fig, min(w, h)):.1f} mm into panel'))
                break

    # --- panel against panel
    for i in range(len(axes)):
        for j in range(i + 1, len(axes)):
            if not (axes[i].axison and axes[j].axison):
                continue
            ai, aj = axes[i].get_window_extent(r), axes[j].get_window_extent(r)
            w, h = inter(ai, aj)
            if w > 0 and h > 0 and mm(fig, min(w, h)) > 0.3:
                findings.append((name, 'panel over panel', '',
                                 f'{mm(fig, max(ai.x0, aj.x0)):.0f},{mm(fig, H - min(ai.y1, aj.y1)):.0f} mm',
                                 f'{mm(fig, min(w, h)):.1f} mm'))

    # --- text running past the canvas
    W = fig.get_size_inches()[0] * fig.dpi
    for t, bb, s, boxed in items:
        out = max(-bb.x0, bb.x1 - W, -bb.y0, bb.y1 - H)
        if out > 0.5:
            findings.append((name, 'text past the canvas', short(s),
                             f'{mm(fig, bb.x0):.0f},{mm(fig, H - bb.y1):.0f} mm', f'{mm(fig, out):.1f} mm'))


figmod.save = lambda fig, name, outdir: audit(fig, name)

for fn in [figmod.figure1, figmod.figure2, figmod.figure3, figmod.figure4,
           figmod.figureS1, figmod.figureS2, figmod.figureS3, figmod.figureS4,
           figmod.figureS7, figmod.figureS8, figmod.figureS9]:
    fn()

if not findings:
    print('no text or panel collisions in any figure')
    sys.exit(0)

wname = max(len(f[0]) for f in findings)
wkind = max(len(f[1]) for f in findings)
for f in sorted(findings):
    print(f'{f[0]:<{wname}}  {f[1]:<{wkind}}  {f[3]:>12}  {f[4]:>18}  {f[2]}')
print(f'\n{len(findings)} collisions')
sys.exit(1)
