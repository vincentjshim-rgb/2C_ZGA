# Human decomposition plan — addendum 1 (2026-09-27, before any human clock score)

**Status.** Written before `src/31_human_decomposition.py` produced any clock value: the first run failed inside the
preprocessing with an empty matrix and no score was computed or seen. The questions, the prespecified interval and
all three reading rules of the frozen plan are unchanged.

**Subject.** The frozen plan states: "Human is one of the clock's four training species, so no ortholog mapping is
needed; gene symbols are mapped to Entrez identifiers." That is wrong about the data format.

## What is actually the case

The clock's 10,487 input features are **mouse** Entrez identifiers. Checked directly against the tables the package
ships: 10,487 of 10,487 features are found in `Gene_table_mouse.csv`, and 0 of 10,487 in the human, rat or monkey
tables. Human expression must therefore be carried onto mouse identifiers before it can be scored, exactly as the
model's own training pipeline did.

## Correction

Human gene symbols are mapped in three declared steps, all with tables shipped inside the clock package:

1. symbol → human Entrez (`Gene_table_human.csv`),
2. human Entrez → mouse Entrez (`Table_of_orthologs.csv`, columns `Entrez.Human` and `Entrez.Mouse`),
3. mouse Entrez → mouse Ensembl (`Gene_table_mouse.csv`), so the matrix enters `tage_py.preprocess` in the same form
   as every other dataset in this study and the validated preprocessing path runs unchanged.

Where several human genes map to one mouse gene their counts are summed, as for the other cross-species datasets.

## What this changes about the reading

Nothing in the rules, but one sentence of the plan's rationale no longer holds and is withdrawn: the human analysis
is **not** free of ortholog mapping. It is, however, better placed than the cow, pig and rabbit analyses of Gate 1,
which used one-to-one Ensembl orthologs assembled here; the human mapping uses the correspondence table distributed
with the clock itself, which is the mapping the model was built with. The manuscript must say this rather than repeat
the claim that human needs no mapping.

The number of weighted clock genes detected, which reading 3 of the frozen plan already requires, is now also the
number surviving this mapping, and is reported as such.
