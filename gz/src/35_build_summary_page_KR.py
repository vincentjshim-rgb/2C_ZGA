"""Render the Korean manuscript summary as a reading page, with the five main figures in place.

Input:  results/MANUSCRIPT_SUMMARY_KR.md  (the summary; this script never edits it)
Output: results/review/manuscript_review_KR.html

The page is self-contained: every figure is embedded, re-encoded from `figures/pub` on each build, so it cannot
show a figure that has since been redrawn. Figures are placed after the section that first discusses them.

Run with /usr/bin/python3 (needs markdown and PIL).
"""
import base64
import io
import os
import re

import markdown
from PIL import Image

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
SRC = f'{B}/results/MANUSCRIPT_SUMMARY_KR.md'
OUT = f'{B}/results/review/manuscript_review_KR.html'

# section heading that each main figure follows, and its one-line caption
PLACE = {
    '### 2.1': ('Figure1', '그림 1 — 하락은 있는가, 언제인가, 청소 효과인가'),
    '### 2.3': ('Figure2', '그림 2 — ZGA를 막으면 각 arm의 하락이 어떻게 되는가'),
    '### 2.5': ('Figure3', '그림 3 — 시계 값의 유전자별 분해와 DUX 구제'),
    '### 2.6': ('Figure4', '그림 4 — 시계 가중치를 쓰지 않는 splicing 진행 지표'),
    '### 2.7': ('Figure5', '그림 5 — 사람 배아에서의 같은 분해'),
}


def embed(name, width=1500):
    im = Image.open(f'{B}/figures/pub/{name}.png').convert('RGB')
    im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, 'JPEG', quality=82, optimize=True)
    return 'data:image/jpeg;base64,' + base64.b64encode(buf.getvalue()).decode()


md = open(SRC, encoding='utf-8').read()
title = md.splitlines()[0].lstrip('# ').strip()

# split into blocks at the headings that carry a figure, so the figure lands after that whole section
parts, order = [], list(PLACE)
chunks = re.split(r'\n(?=### 2\.\d)', md)
html_parts = []
for ch in chunks:
    key = next((k for k in order if ch.startswith(k)), None)
    html_parts.append(markdown.markdown(ch, extensions=['tables']))
    if key:
        name, cap = PLACE[key]
        html_parts.append(f'<figure class="fig"><img src="{embed(name)}" alt="{cap}" loading="lazy">'
                          f'<figcaption>{cap}</figcaption></figure>')
body = '\n'.join(html_parts)

HTML = f'''<title>{title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Gowun+Batang:wght@400;700&family=Noto+Sans+KR:wght@300;400;500;700&family=IBM+Plex+Mono:wght@400&display=swap">
<style>
:root {{
  --paper:#FAF8F4; --card:#FFFDF9; --ink:#211E1A; --ink2:#413C34; --muted:#6F6A61;
  --rule:#E4DFD5; --rule2:#D6D0C3; --accent:#7A4A1E; --band:#F3ECDD;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper:#16150F; --card:#1E1C16; --ink:#EDE8DE; --ink2:#CFC8BA; --muted:#9A9387;
    --rule:#332F26; --rule2:#433E33; --accent:#D9A05B; --band:#2A2519;
  }}
}}
:root[data-theme="dark"] {{
  --paper:#16150F; --card:#1E1C16; --ink:#EDE8DE; --ink2:#CFC8BA; --muted:#9A9387;
  --rule:#332F26; --rule2:#433E33; --accent:#D9A05B; --band:#2A2519;
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:'Noto Sans KR',system-ui,sans-serif; font-weight:300; line-height:1.78;
  font-size:15.5px; letter-spacing:-0.002em;
}}
.wrap {{ max-width:860px; margin:0 auto; padding:0 20px 90px; }}
header {{ border-bottom:1px solid var(--rule); background:var(--card); }}
header .inner {{ max-width:860px; margin:0 auto; padding:14px 20px; font-size:12.5px; color:var(--muted); }}
h1 {{
  font-family:'Gowun Batang',serif; font-weight:700; font-size:30px; line-height:1.4;
  margin:40px 0 6px; text-wrap:balance; letter-spacing:-0.01em;
}}
h2 {{
  font-family:'Gowun Batang',serif; font-weight:700; font-size:21px; margin:52px 0 4px;
  padding-top:22px; border-top:1px solid var(--rule2);
}}
h3 {{ font-weight:700; font-size:16.5px; margin:34px 0 2px; color:var(--accent); }}
p {{ margin:12px 0; }}
hr {{ border:0; border-top:1px solid var(--rule); margin:40px 0; }}
strong {{ font-weight:700; }}
blockquote {{
  margin:22px 0; padding:14px 18px; background:var(--band); border-left:3px solid var(--accent);
  border-radius:0 6px 6px 0; font-size:14.5px;
}}
blockquote p {{ margin:4px 0; }}
ol, ul {{ padding-left:22px; }}
li {{ margin:7px 0; }}
.tablewrap, table {{ width:100%; }}
table {{
  border-collapse:collapse; margin:20px 0; font-size:13.5px; display:block; overflow-x:auto;
  white-space:normal;
}}
th, td {{ border-bottom:1px solid var(--rule); padding:8px 11px; text-align:left; vertical-align:top; }}
th {{ font-weight:500; color:var(--muted); font-size:12.5px; border-bottom:1px solid var(--rule2); }}
td:first-child {{ font-weight:500; }}
figure.fig {{
  margin:30px 0; padding:16px; background:var(--card); border:1px solid var(--rule); border-radius:8px;
}}
figure.fig img {{ width:100%; display:block; border-radius:3px; }}
figcaption {{ margin-top:10px; font-size:12.5px; color:var(--muted); }}
a {{ color:var(--accent); text-decoration-thickness:1px; text-underline-offset:2px; }}
code {{ font-family:'IBM Plex Mono',monospace; font-size:0.88em; background:var(--band); padding:1px 5px; border-radius:3px; }}
@media (max-width:620px) {{
  h1 {{ font-size:24px; }} h2 {{ font-size:19px; }} body {{ font-size:15px; }}
}}
</style>
<header><div class="inner">원고 한국어 요약 · 투고본의 모든 수치를 그대로 옮긴 것 · 차이가 있으면 영문본이 맞다</div></header>
<div class="wrap">
{body}
</div>
'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, 'w', encoding='utf-8').write(HTML)
print(f'wrote {OUT} ({os.path.getsize(OUT) // 1024} KB, {len(PLACE)} figures embedded)')
