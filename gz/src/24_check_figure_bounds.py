"""Check every rendered figure for content cut off at the canvas border.

Matplotlib silently clips text that extends past the figure; the clipped words are simply absent from the exported
file, so the defect is invisible unless the border is inspected. Every panel in these figures has a margin, so any
non-white pixel within a few pixels of the edge means something was cut. Run after src/18_figures_publication.py.

Exit status 1 if any figure fails, so it can gate a re-render.
"""
import glob
import os
import sys

import numpy as np
from PIL import Image

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
B = f'{ANALYSIS}/gz'
BAND = 3      # pixels from the edge to inspect
INK = 245     # 8-bit grey below this counts as ink (the figures are drawn on white)


def border_ink(path):
    a = np.array(Image.open(path).convert('L'))
    sides = {'top': a[:BAND, :], 'bottom': a[-BAND:, :], 'left': a[:, :BAND], 'right': a[:, -BAND:]}
    return {k: int((v < INK).sum()) for k, v in sides.items() if (v < INK).any()}


def main():
    files = sorted(glob.glob(f'{B}/figures/pub/Figure*.png')) + sorted(glob.glob(f'{B}/figures/supp/FigureS*.png'))
    if not files:
        sys.exit('no rendered figures found; run src/18_figures_publication.py first')
    failed = []
    for p in files:
        bad = border_ink(p)
        name = os.path.basename(p)
        if bad:
            failed.append(name)
            print(f'{name:16s} CLIPPED  {bad}')
        else:
            print(f'{name:16s} clean')
    if failed:
        print(f'\n{len(failed)} of {len(files)} figures have content cut off at the canvas border: {", ".join(failed)}')
        sys.exit(1)
    print(f'\nall {len(files)} figures clear of the canvas border')


if __name__ == '__main__':
    main()
