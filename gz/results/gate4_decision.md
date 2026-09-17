# Gate 4 decision

VERDICT: **ZSA ZGA-DEPENDENT**

ZSA events per dataset: {'GSE280522': 511, 'GSE221985': 1376, 'GSE300734': 77}; ZSA PRESENT: True

Interactions (progression, blocking arm): {'GSE280522': -0.4431907824989859, 'GSE221985': -0.42143038632274177, 'GSE300734': -0.45115706849284676}; P1 rescue: 0.23013972148736372

## ZSA summary

      gse event_type  n_filtered  n_zsa  frac_zsa  n_events_total  n_libraries
GSE280522         A3        1868     58    0.0310             NaN          NaN
GSE280522         A5        1741     63    0.0362             NaN          NaN
GSE280522         AF        4162    166    0.0399             NaN          NaN
GSE280522         AL         795     37    0.0465             NaN          NaN
GSE280522         MX         241      8    0.0332             NaN          NaN
GSE280522         RI         832     24    0.0288             NaN          NaN
GSE280522         SE        4076    155    0.0380             NaN          NaN
GSE280522        ALL       13715    511    0.0373         88526.0         23.0
GSE221985         A3        1726    157    0.0910             NaN          NaN
GSE221985         A5        1687    194    0.1150             NaN          NaN
GSE221985         AF        4014    527    0.1313             NaN          NaN
GSE221985         AL         769     75    0.0975             NaN          NaN
GSE221985         MX         212     21    0.0991             NaN          NaN
GSE221985         RI         757     41    0.0542             NaN          NaN
GSE221985         SE        3888    361    0.0928             NaN          NaN
GSE221985        ALL       13053   1376    0.1054         88526.0          8.0
GSE300734         A3         635      4    0.0063             NaN          NaN
GSE300734         A5         590      6    0.0102             NaN          NaN
GSE300734         AF        2394     37    0.0155             NaN          NaN
GSE300734         AL         311      6    0.0193             NaN          NaN
GSE300734         MX         103      0    0.0000             NaN          NaN
GSE300734         RI         241      0    0.0000             NaN          NaN
GSE300734         SE        1458     24    0.0165             NaN          NaN
GSE300734        ALL        5732     77    0.0134         88526.0         11.0

## Arm progression

      gse                           quantity   value   ci_lo   ci_hi  n_E2C  n_L2C
GSE280522                          P control  1.0000  0.8138  1.1828    4.0    4.0
GSE280522                             P A485  0.5568  0.4245  0.6752    4.0    3.0
GSE280522                         P A485+Dux  0.7869  0.5674  0.9986    4.0    4.0
GSE280522         INTERACTION A485 - control -0.4432 -0.6694 -0.2196    NaN    NaN
GSE280522     INTERACTION A485+Dux - control -0.2131 -0.5042  0.0717    NaN    NaN
GSE280522             RESCUE A485+Dux - A485  0.2301 -0.0161  0.4734    NaN    NaN
GSE221985                          P control  1.0000  0.9762  1.0238    2.0    2.0
GSE221985                     P Tardbp_matKO  0.5786  0.5190  0.6381    2.0    2.0
GSE221985 INTERACTION Tardbp_matKO - control -0.4214 -0.4880 -0.3548    NaN    NaN
GSE300734                          P control  1.0000  0.8176  1.1824    2.0    3.0
GSE300734                       P Brg1_matKO  0.5488  0.4574  0.6403    3.0    3.0
GSE300734   INTERACTION Brg1_matKO - control -0.4512 -0.6604 -0.2320    NaN    NaN

## Clock coupling (descriptive, no inference)

Spearman, arm x stage group means, progression vs Gate 3 V0 tAge: {'GSE280522': np.float64(-0.942857142857143), 'GSE221985': np.float64(-0.39999999999999997), 'GSE300734': np.float64(-1.0), 'pooled': np.float64(-0.6028606508323641)}

      gse                        interaction  progression_I  clock_I_V0
GSE280522         INTERACTION A485 - control        -0.4432      0.1809
GSE280522     INTERACTION A485+Dux - control        -0.2131      0.1738
GSE221985 INTERACTION Tardbp_matKO - control        -0.4214      0.0415
GSE300734   INTERACTION Brg1_matKO - control        -0.4512      0.1091

Spearman of interactions: 0.0
