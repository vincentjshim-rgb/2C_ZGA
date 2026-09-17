# Gate 1 plan — addendum 1 (2026-09-14, written before any embryo clock score was produced)

Loading bugs found on first run (the run stopped before prediction for cow and rhesus; no tAge value was printed
or saved for any species). Fixes change identifier mapping only; stages, intervals, clock, statistics and decision
rules of `GATE1_PLAN_frozen.md` are unchanged.

1. Cow matrix uses RefSeq transcript accessions, not Ensembl gene IDs (28,906 NM_/XM_ rows; 37 ENS rows).
   Mapping: NM_ (version stripped) → cow Ensembl gene via Ensembl BioMart `refseq_mrna` xref (97.3% of NM_ rows);
   XM_ accessions are all "suppressed" at NCBI and absent from Ensembl xrefs, so XM_ → gene symbol parsed from the
   NCBI nuccore title (last parenthesised token; LOC* symbols dropped) → cow Ensembl gene by exact `external_gene_name`
   match (unique names only). Counts of transcripts mapping to the same cow gene are summed; cow gene → mouse via
   Ensembl one-to-one orthologs as for pig and rabbit. Cow coverage is reported and flagged as a species-specific
   limitation; a cow-excluded decision is reported as a sensitivity row.
2. Rhesus: the tAge package ortholog table is Macaca fascicularis (ENSMFAG), whereas the matrix is M. mulatta
   (ENSMMUG). Rhesus → mouse uses Ensembl BioMart `mmulatta_gene_ensembl` one-to-one orthologs. Rhesus stays descriptive.
3. A species is scored only if it has >= 1 library in every ordered stage after QC.
