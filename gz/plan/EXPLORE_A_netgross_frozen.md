# Exploration A — reading the clock twice: net and gross (frozen 2026-10-12, before any gross value was computed)

**Status.** First of four exploratory analyses (A–D) run to choose the paper's core logic. Exploratory means: the
outcome decides direction, not a gate verdict; anything from here that enters the paper is labelled as planned on this
date, after all earlier gate verdicts were known. No new data; every dataset below has already been scored.

## Question

A linear clock difference is D = Σᵢ cᵢ Δᵢ. The same terms give a second number, G = Σᵢ |cᵢ Δᵢ|, which cannot cancel.
We have shown D to be a small residual of opposing sums in the two-cell window. Does G — the clock's own
coefficient-weighted measure of how much its genes were remodelled — locate zygotic genome activation where D does
not, and does it register the perturbation and the rescue that D conceals?

## Definitions (fixed)

For a comparison of two groups of libraries (late vs early stage; or within an arm), with the 1,839 non-zero
coefficients cᵢ of the multi-species chronological clock and group means x̄ᵢ of the preprocessed feature
(`scaled_diff`, computed once per dataset exactly as in the earlier gates, reference as before):

- termᵢ = cᵢ (x̄ᵢ,late − x̄ᵢ,early); genes absent from the dataset contribute exactly 0
- **net** D = Σ termᵢ (the clock difference already reported)
- **gross** G = Σ |termᵢ|
- down = Σ_{termᵢ<0} termᵢ, up = Σ_{termᵢ>0} termᵢ (so D = down + up, G = up − down)
- **cancellation** C = |D| / G (1 = no cancellation, 0 = complete)

G scales with the number of clock genes a dataset detects, so G is compared only within a dataset, never across.

Uncertainty: 2,000 bootstrap resamples of libraries within stage (or arm × stage), seed 20261012, percentile 95%
intervals; for "which interval is maximal", the bootstrap probability P_max that the named interval has the largest G.

## Datasets and intervals (as in the frozen gates; nothing added)

- GSE225056, four species, the 24 consecutive-stage intervals of Gate 1 (6 per species).
- GSE45719 (mouse, embryo pseudobulks) and GSE66582 (mouse bulk, n = 2 per stage): consecutive-stage intervals.
- GSE280522 (control, A485, A485+DUX), GSE221985 (control, maternal Tardbp KO), GSE300734 (control, maternal Brg1
  KO): early → late two-cell within each arm, all genes (V0).
- GSE36552 (human, 20 embryo pseudobulks): the six consecutive-stage intervals from oocyte to late blastocyst.

## Predictions, written before computing

**P1 — species.** In each species the interval of maximal G is the literature major-ZGA interval fixed in the Gate 1
plan (mouse early → late two-cell; pig 2-cell → 4-cell; cow and rabbit 8-cell → 16-cell), with P_max ≥ 0.6.
*Supported* if ≥ 3 of 4 species; *unsupported* otherwise. Secondary, reported only: coincidence with the interval of
maximal measured zygotic gain (Gate 1-alt). If G is maximal at fertilisation (oocyte → zygote) or at a division
interval, that is the result and is reported as such.

**P2 — perturbation and rescue.** In each primary dataset, G is lower in the perturbed arm than in the control arm
(ΔG = G(perturbed) − G(control) < 0). *Supported* if the sign holds in all three; the bootstrap interval of ΔG is
reported for each. Rescue: G(A485+DUX) − G(A485) > 0 with a bootstrap interval excluding zero. The cancellation C is
reported per arm; no rule is attached to it.

**P3 — human.** In GSE36552 the interval of maximal G is 4-cell → 8-cell, the literature interval of human major
embryonic genome activation (Braude et al., 1988), whereas the net decrease was at 8-cell → morula. P_max ≥ 0.6.

**P4 — mouse replication series.** In GSE45719 the interval of maximal G lies within the two-cell stage; in GSE66582
(direction only, n = 2) it is early two-cell → two-cell.

## How the outcome is read for the direction decision

- P1 supported and P2 supported → A is a candidate core logic ("the scalar misses the switch; the components find
  it"), to be weighed against B–D.
- One of the two supported → A is supporting material, not a core logic.
- Neither → the net/gross reading is dropped; the cancellation remains what it is now, a property of D.

## Outputs

`results/explore_A_netgross_intervals.tsv` (dataset, species, interval, n_early, n_late, D, G, down, up, C, CIs,
P_max), `results/explore_A_netgross_arms.tsv` (dataset, arm, D, G, down, up, C, ΔG vs control with CI),
`results/explore_A_netgross.md` (readings P1–P4 against the rules above), figure data for a candidate Figure 6.
Script `src/42_explore_A_netgross.py`.
