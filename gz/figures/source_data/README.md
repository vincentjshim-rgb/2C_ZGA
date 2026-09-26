# Figure source data

File names carry the panel numbering of the original figure plan (five main figures, supplementary figures in plan
order). The manuscript merged the plan's Figures 3 and 4 into one figure and renumbered the supplementary figures by
first citation, so the final panels are:

| File | Final panel |
|---|---|
| `gate1_embryo_tage.tsv`, `gate1_intervals.tsv` | Figure 1b |
| `gate2b_R1_embryo_tage.tsv` | Figure 1c |
| `gate2a_intervals.tsv`, `gate2a_simulation.tsv` | Figure 1d |
| `fig2b_clock_per_library.tsv`, `gate3_arm_drops.tsv` | Figure 2b |
| `gate3_secondary_single_stage.tsv` | Figure S4 (secondary arms) |
| `fig3a_4c_cumulative_contributions.tsv` | Figure 3a (P1/P3 control) and Figure 3f (three arms, control order) |
| `posthoc_gate3_contribution_genes.tsv`, `posthoc_gate3_contribution_genes_P3.tsv` | Figure 3b, 3c |
| `fig4a_heatmap_z.tsv`, `fig4a_heatmap_columns.tsv` | Figure 3d |
| `fig4b_zygotic_score_per_library.tsv` | Figure 3e |
| `posthoc_rescue_contributions.tsv` | Figure 3f; Figure S6 |
| `fig5a_psi_zsa_events_P1.tsv` | Figure 4a |
| `fig5b_dpsi_per_arm_P1.tsv` | Figure 4b |
| `gate4_progression.tsv`, `gate4_arm_progression.tsv` | Figure 4c |
| `posthoc_gate4_crossfit_folds.tsv` | Figure 4d |
| `gate1alt_intervals.tsv` | Figure S1 |
| `gate2b_R2_library_tage.tsv`, `gate2b_R2_intervals.tsv` | Figure S2 |
| `figS7_library_qc.tsv` | Figure S3 |
| `fig3d_key_genes_per_library.tsv` | Figure S5 |
| `figS8_volcano_events.tsv` | Figure S7 (splicing-event selection) |
| `fig5a_psi_pca.tsv` | Figure S8 |

Written by `src/17_figdata.py`, `src/19_figdata_supp.py`, `src/22_figdata_v2.py` and `src/25_figdata_volcano.py`;
every re-extracted value is checked against the stored result before it is written. Figures are drawn by
`src/18_figures_publication.py` from these files only, and two checks run after every render:
`src/24_check_figure_bounds.py` (no ink at the canvas border) and `src/26_check_figure_overlaps.py` (no text or
panel collisions).
