#!/usr/bin/env python3
"""Cut the published figures into single panels for the general-reader page.

A whole six-panel figure is unreadable for someone meeting the subject for the first time; each panel is
cropped by its millimetre coordinates in the figure script so the explanation can sit next to the one panel
it is about.
"""
import os

from PIL import Image

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
ASSETS = f'{B}/figures/review_assets'
OUT = f'{ASSETS}/pubcrop'
os.makedirs(OUT, exist_ok=True)
WMM = 167.0

SRC = {'F1': (f'{B}/figures/pub/Figure1.png', 254),
       'F2': (f'{B}/figures/pub/Figure2.png', 132),
       'F3': (f'{B}/figures/pub/Figure3.png', 228),
       'F4': (f'{B}/figures/pub/Figure4.png', 186),
       'F5': (f'{B}/figures/pub/Figure5.png', 72),
       'S7': (f'{B}/figures/supp/FigureS7.png', 156)}

# name: (figure, x0, y0, x1, y1 in mm, output width in px)
CROPS = {
    'f1b_mouse': ('F1', 0, 63, 57, 117, 900),
    'f1b_all':   ('F1', 0, 63, 167, 117, 1500),
    'f1c':       ('F1', 0, 116, 101, 181, 1200),
    'f1d':       ('F1', 99, 116, 167, 181, 900),
    'f1e':       ('F1', 0, 180, 167, 246, 1500),   # 180, not 174: panel d's footnote runs to y 178
    'f2a':       ('F2', 0, 0, 167, 58, 1500),
    'f2b':       ('F2', 0, 58, 167, 132, 1500),
    'f3a':       ('F3', 6, 0, 90, 53, 1100),       # x from 6: panel c's letter sits at x 1, y 54
    'f3b':       ('F3', 90, 0, 167, 102, 1000),
    'f3c':       ('F3', 0, 52, 90, 104, 1000),
    'f3d':       ('F3', 0, 100, 167, 182, 1500),
    'f3e':       ('F3', 0, 182, 78, 228, 950),
    'f3f':       ('F3', 78, 184, 167, 228, 1100),
    'f4a':       ('F4', 0, 0, 132, 80, 1300),
    'f4b':       ('F4', 128, 0, 167, 80, 600),
    'f4c':       ('F4', 0, 78, 104, 140, 1250),    # 104: panel d's row labels reach left to about 103
    'f4d':       ('F4', 105, 78, 167, 140, 800),
    'f4e':       ('F4', 0, 140, 167, 186, 1500),
    'f5a':       ('F5', 0, 0, 60, 61, 900),   # 61: the n row sits below the tick labels
    'f5b':       ('F5', 57, 0, 118, 56, 950),
    'f5c':       ('F5', 118, 4, 167, 52, 800),
    's7a':       ('S7', 0, 0, 90, 70, 1000),
    's7b':       ('S7', 74, 0, 167, 146, 1150),
}

for name, (fig, x0, y0, x1, y1, w) in CROPS.items():
    path, hmm = SRC[fig]
    im = Image.open(path).convert('RGB')
    sx, sy = im.width / WMM, im.height / hmm
    box = (int(x0 * sx), int(y0 * sy), int(min(x1, WMM) * sx), int(min(y1, hmm) * sy))
    c = im.crop(box)
    c = c.resize((w, max(1, round(c.height * w / c.width))), Image.LANCZOS)
    c.save(f'{OUT}/{name}.jpg', 'JPEG', quality=84, optimize=True)
    print(f'{name:12s} {x1 - x0:5.0f} x {y1 - y0:4.0f} mm -> {c.width}x{c.height}px'
          f'  {os.path.getsize(f"{OUT}/{name}.jpg") // 1024:4d} KB')

# whole figures, small, for the reference section at the end
for tag, (path, hmm) in SRC.items():
    im = Image.open(path).convert('RGB')
    w = 1100
    im.resize((w, round(im.height * w / im.width)), Image.LANCZOS).save(
        f'{OUT}/whole_{tag}.jpg', 'JPEG', quality=80, optimize=True)
print('\nwhole figures written:', ', '.join(f'whole_{t}' for t in SRC))
print('total', sum(os.path.getsize(f'{OUT}/{f}') for f in os.listdir(OUT)) // 1024, 'KB')
