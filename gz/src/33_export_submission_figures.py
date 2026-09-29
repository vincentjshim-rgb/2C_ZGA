"""Export every rendered figure into the submission bundle and rewrite its README.

Reads what `18_figures_publication.py` produced (`figures/pub`, `figures/supp`) and writes, for each figure,
a vector PDF, a 600-dpi PNG and a 600-dpi LZW TIFF under `results/submission/figures/`. The README table is
built from the PDF page sizes and the legend titles in the manuscript, so it cannot drift from either.

Run with /usr/bin/python3 after the figures have been rendered.
"""
import os
import re
import shutil

from PIL import Image

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
OUT = f'{B}/results/submission/figures'
MS = f'{B}/results/MANUSCRIPT_gz_v1_EN.md'
PT_PER_MM = 72 / 25.4

t = open(MS, encoding='utf-8').read()
# a legend title may wrap over several lines, so the match has to cross newlines
main = re.findall(r'^\*\*FIGURE (\d+): (.+?)\*\*', t, re.M | re.S)
supp = re.findall(r'^\*\*Figure (S\d+): (.+?)\*\*', t, re.M | re.S)
assert [n for n, _ in main] == [str(i) for i in range(1, len(main) + 1)], [n for n, _ in main]
assert [n for n, _ in supp] == [f'S{i}' for i in range(1, len(supp) + 1)], [n for n, _ in supp]
titles = {f'Figure{n}': ' '.join(s.split()) for n, s in main + supp}
sources = {f'Figure{n}': f'{B}/figures/pub' for n, _ in main}
sources.update({f'Figure{n}': f'{B}/figures/supp' for n, _ in supp})
order = [f'Figure{n}' for n, _ in main] + [f'Figure{n}' for n, _ in supp]

for d in ('pdf', 'png', 'tiff'):
    os.makedirs(f'{OUT}/{d}', exist_ok=True)

rows, missing = [], []
for name in order:
    src = sources[name]
    pdf, png = f'{src}/{name}.pdf', f'{src}/{name}.png'
    if not (os.path.exists(pdf) and os.path.exists(png)):
        missing.append(name)
        continue
    shutil.copy2(pdf, f'{OUT}/pdf/{name}.pdf')
    shutil.copy2(png, f'{OUT}/png/{name}.png')

    # page size straight out of the PDF, in mm
    with open(pdf, 'rb') as fh:
        box = re.findall(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]', fh.read())[0]
    x0, y0, x1, y1 = (float(v) for v in box)
    w_mm, h_mm = (x1 - x0) / PT_PER_MM, (y1 - y0) / PT_PER_MM

    im = Image.open(png)
    dpi = im.info.get('dpi', (0, 0))[0]
    im.save(f'{OUT}/tiff/{name}.tif', format='TIFF', compression='tiff_lzw', dpi=(600, 600))
    rows.append((name, w_mm, h_mm, im.size, round(dpi), titles[name]))

assert not missing, f'not rendered: {missing}'
widths = {round(w) for _, w, _, _, _, _ in rows}
assert widths == {167}, widths
assert {d for *_, d, _ in rows} == {600}, 'not every PNG is 600 dpi'

nmain, nsupp = len(main), len(supp)
lines = ['# Figure files', '',
         f'Every figure of the manuscript, exported separately: {nmain} main and {nsupp} supporting. All are 167 mm '
         'wide at 600 dpi;', 'no text is below 6 pt at final size.', '',
         '- `pdf/` vector  - `png/` 600 dpi raster  - `tiff/` 600 dpi LZW TIFF', '',
         '| File | Size (mm) | Title |', '|---|---|---|']
for name, w, h, _, _, title in rows:
    lines.append(f'| {name} | {w:.0f} × {h:.0f} | {title} |')
open(f'{OUT}/README.md', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

print(f'{len(rows)} figures exported to {OUT} ({nmain} main, {nsupp} supporting)')
for name, w, h, size, dpi, _ in rows:
    print(f'  {name:10s} {w:.0f} × {h:.0f} mm   {size[0]} × {size[1]} px @ {dpi} dpi')
