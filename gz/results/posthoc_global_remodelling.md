# Post hoc: global transcriptome remodelling vs gene-specific change (plan: plan/POSTHOC_global_remodelling_frozen.md)

VERDICT (control arms, prespecified rule): **EXPLAINED** by global remodelling

P1 interaction (A485 - control): **NOT SUPPORTED** by global remodelling

D = G + S; G = (sum of coefficients) x (mean feature change) is the drop expected if every clock gene moved by the average amount; S is the gene-specific remainder. The null re-assigns the coefficient vector to genes 2,000 times (seed 20260920), holding the observed feature changes fixed.

quantity                                           D_observed  G_global_component  S_specific_component  abs_S_over_abs_G  null_mean  null_p2.5  null_p97.5  fraction_of_null_below_observed  D_features_down_only  D_features_up_only  D_expression_rises  D_expression_falls
dataset      arm                                                                                                                                                                                                                                                              
P1_GSE280522 A485                                     -0.0566             -0.0094               -0.0471            4.9956    -0.0054    -0.2289      0.2193                           0.3230               -0.0092             -0.0474              0.0107             -0.0673
             A485+Dux                                 -0.0637              0.0047               -0.0683           14.6748     0.0053    -0.2467      0.2653                           0.3005                0.0181             -0.0818             -0.0508             -0.0128
             control                                  -0.2375              0.0072               -0.2447           34.0666     0.0042    -0.2399      0.2505                           0.0260               -0.0208             -0.2167             -0.2149             -0.0225
             interaction_A485_minus_control            0.1809             -0.0166                0.1975           11.8903    -0.0126    -0.2141      0.1889                           0.9675                0.1677              0.0132              0.2257             -0.0447
P3_GSE300734 Brg1_matKO                               -0.0623              0.0088               -0.0711            8.0657     0.0096    -0.2075      0.2284                           0.2635                0.0385             -0.1008             -0.0914              0.0291
             control                                  -0.1714              0.0122               -0.1836           15.0729     0.0138    -0.2558      0.2645                           0.0805               -0.0219             -0.1496             -0.1483             -0.0231
             interaction_Brg1_matKO_minus_control      0.1091             -0.0034                0.1125           33.4620    -0.0071    -0.1848      0.1546                           0.9150                0.0663              0.0428              0.0570              0.0522
