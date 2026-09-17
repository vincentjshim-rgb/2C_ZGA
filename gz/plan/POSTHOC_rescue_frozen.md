# POST HOC — why the A485+Dux arm does not restore the clock decrease
(frozen 2026-09-17, before running `src/14_posthoc_rescue.py`; written after Gate 3 and after the contribution
decomposition were seen)

**Status.** POST HOC and descriptive, like `POSTHOC_gate3_contribution_frozen.md`. It does not change Gate 3
(MIXED) and does not re-open the rescue rule. It asks one question with two possible answers that are distinguishable
in the data.

## Question

In P1 GSE280522 the rescue term was R = D(A485+Dux) − D(A485) = −0.0071 [−0.0629, 0.0572]: adding Dux to A485 left the
drop where A485 put it. Two explanations are separable:

- **(a) The rescue did not happen biologically.** Dux co-expression did not restore zygotic transcription in this arm,
  so the clock had nothing to respond to. The Gate 3 rescue term is then uninformative about the clock.
- **(b) The rescue happened but not on the genes the clock reads.** Zygotic transcription was restored while the
  contributions identified in the decomposition stayed absent.

## Prespecified measurements (GSE280522 only, QC-passing libraries, same preprocessing as Gate 3)

1. **Zygotic-set expression, clock-independent.** Zygotic genes defined by the Gate 3 V2 rule in the control arm
   (CPM < 1 in ≥ 80% of control E2C and CPM ≥ 10 in ≥ 50% of control L2C), Ensembl space. Per library, the score is the
   mean log2(CPM+1) over that set. Compared at L2C between control, A485 and A485+Dux; maternal-set score reported
   beside it as a contrast. Bootstrap of libraries, 2000 replicates, seed 20260918.
2. **Known 2-cell / Dux-target markers**, reported per arm at L2C as log2(CPM+1): all genes whose symbol matches
   `Zscan4*`, `Tcstv*`, `Tdpoz*`, `Zfp352`, and `Dux*` itself. Fixed list, written before looking.
3. **Clock contributions of the same genes as before.** The decomposition of `13_posthoc_gate3_contribution.py`
   extended to the A485+Dux arm: for the 20 genes with the largest control down-contribution, `c` in control, A485 and
   A485+Dux, and their sums.

## Prespecified reading

- Zygotic-set score in A485+Dux near control and clearly above A485 → restoration happened → explanation (b), and the
  gap is then between what Dux restores and what the clock weights.
- Zygotic-set score in A485+Dux near A485 and clearly below control → restoration did not happen in this arm →
  explanation (a); the Gate 3 rescue term says nothing about the clock, and that is what will be reported.
- Anything between the two → reported as partial, with both numbers, and neither explanation is asserted.

## Limits

Descriptive; bootstrap intervals are reported for the arm differences and nothing else. A transgene delivered as an
expression construct need not be visible in a kallisto quantification against the Ensembl transcriptome, so a low
endogenous `Dux` value is not evidence that the construct was absent; it is reported, not interpreted. No causal verb.
No gene is promoted to a claim. The result does not change the Gate 3 verdict either way.
