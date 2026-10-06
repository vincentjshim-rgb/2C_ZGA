"""Build the Aging Cell submission draft from MANUSCRIPT_gz_v1_EN.md.
Changes made here (nothing else is edited):
  - title page with the journal's required elements (placeholders where the author must decide),
  - references rewritten in APA 7th with full author lists and volume(issue) from Crossref,
  - the data/code availability text moved out of Methods into a statement section,
  - the statement sections Aging Cell prints after the Methods added as templates.
Output: results/MANUSCRIPT_gz_v1_EN_submission.md. Run with /usr/bin/python3; then 40, then 34."""
import json, re

import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
R = f'{B}/results'
t = open(f'{R}/MANUSCRIPT_gz_v1_EN.md', encoding='utf-8').read()
cr = json.load(open(f'{B}/metadata/refs_crossref.json'))

# ---- references: APA 7th, full author lists, volume(issue) -------------------------------------------
body, refs_block = t.split('## References', 1)
refs_block, tail_note = refs_block.split('\n---', 1) if '\n---' in refs_block else (refs_block, '')
entries = [e.strip() for e in refs_block.split('\n\n') if e.strip() and re.match(r'^[A-ZÁ-Ý]', e.strip())]
note = [e.strip() for e in refs_block.split('\n\n') if e.strip().startswith('*Reference check')]


def initials(given):
    return ' '.join(p[0].upper() + '.' for p in re.split(r'[\s.\-]+', given.strip()) if p)


def apa_authors(auth):
    names = [f'{f}, {initials(g)}' if g else f for f, g in auth if f]
    if not names:
        return None
    if len(names) == 1:
        return names[0]
    if len(names) <= 20:
        return ', '.join(names[:-1]) + ', & ' + names[-1]
    return ', '.join(names[:19]) + ', … ' + names[-1]


new_refs = []
for e in entries:
    doi = re.search(r'https://doi\.org/(\S+)', e).group(1)
    ym = re.search(r'\((\d{4})\)\.', e)
    meta = cr.get(doi, {})
    rest = e[ym.end():].strip()
    # Aging Cell follows APA 6th and states: "issue numbers should not be included unless each
    # issue in the volume begins with page one" — so no issue numbers are inserted.
    a = apa_authors(meta.get('authors', []))
    new_refs.append(f'{a} ({ym.group(1)}). {rest}' if a else e)

# ---- move the data/code availability text out of the Methods ----------------------------------------
m = re.search(r'### 4\.1[0-9] Data and code availability\n\n(.*?)\n\n---', body, re.S)
avail = m.group(1).strip()
body = body[:m.start()] + '---' + body[m.end():]

TITLE = re.search(r'^# (.+)$', t, re.M).group(1).strip()
TITLE_PAGE = """**Title:** __TITLE__

**Running title:** __RUNNING__

**Authors:** [TO BE COMPLETED — full names in submission order, with superscript affiliation numbers and ORCID iDs]

**Affiliations:** [TO BE COMPLETED — numbered affiliations]

**Correspondence:** [TO BE COMPLETED — name, department, institution, postal address, e-mail]

**Article type:** Original Research Article

**Counts:** abstract {ABS} words; main text {MAIN} words (Introduction, Results, Discussion and Methods), with a
further {LEG} words of figure legends and {TAB} words of tables; 5 figures; 2 tables; 7 supporting figures;
5 supporting tables; {NREF} references.

---

"""

STATEMENTS = """
## Author contributions

[TO BE COMPLETED — CRediT roles per author. Draft for a single-author submission: conceptualisation, methodology,
software, formal analysis, data curation, writing — original draft, writing — review and editing, visualisation.]

## Funding information

[TO BE COMPLETED — grant numbers and funders, or the statement that the work received no specific funding.]

## Conflict of interest statement

[TO BE COMPLETED — draft: The author declares no conflict of interest.]

## Acknowledgements

[TO BE COMPLETED — acknowledge the depositors of the public datasets and any computational support. The analysis used
only public data; no new experiments were performed.]

## Ethics statement

This study analysed publicly deposited sequencing data only. No new animal or human experiments were performed and no
ethical approval was required.

## Data availability statement

{AVAIL}

## References

"""


def wordcount(text, a, b):
    seg = text.split(a, 1)[1].split(b, 1)[0]
    return len(re.findall(r"[A-Za-z0-9α-ω][\w'’.\-/]*", seg))


head, rest = body.split('## Abstract', 1)
draft_note = re.search(r'\*Draft v1.*?\*\n', head, re.S)
out = TITLE_PAGE + '## Abstract' + rest
out += STATEMENTS.replace('{AVAIL}', avail) + '\n\n'.join(new_refs) + '\n\n'
out += ('\n\n'.join(note) + '\n' if note else '')
abs_wc = wordcount(out, '## Abstract', '**Keywords')
main = sum(wordcount(out, a, b) for a, b in [('## 1. Introduction', '## 2. Results'), ('## 2. Results', '## 3. Discussion'),
                                             ('## 3. Discussion', '## 4. Methods'), ('## 4. Methods', '## Tables')])
leg = wordcount(out, '## Figure legends', '## Supporting Information figure legends')
tab = wordcount(out, '## Tables', '## Figure legends')
running = re.search(r'\*\*Running title:\*\* (.+)', t).group(1).strip()
for k, v in [('{ABS}', str(abs_wc)), ('__TITLE__', TITLE), ('__RUNNING__', running), ('{MAIN}', f'{main:,}'),
             ('{LEG}', f'{leg:,}'), ('{TAB}', f'{tab:,}'), ('{NREF}', str(len(new_refs)))]:
    out = out.replace(k, v)
open(f'{R}/MANUSCRIPT_gz_v1_EN_submission.md', 'w', encoding='utf-8').write(out)
print(f'written; main text {main:,} + legends {leg:,} + tables {tab:,}; abstract {abs_wc}; {len(new_refs)} references')
print('sections:', [l for l in out.split('\n') if l.startswith('## ')])
