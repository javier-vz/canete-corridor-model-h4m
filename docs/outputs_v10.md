# Tables and figures of the simplified manuscript (v10 to v12 candidate) and the files behind them

The simplified version of the article renumbers its tables and figures and moves several of them to the Supplementary
Material (S1–S6). This page gives the file behind each of them. [outputs.md](outputs.md) keeps the numbering of the
previous version. Paths of `eval/`, `ppm/`, `fig/` and `fig_etapa3/` are relative to `out/final/`; paths of
`revision/out/` are relative to the repository root. Files marked *restricted* need the restricted data
([README](../README.md#restricted-data)).

## Main text

| Element | Previous number | File | Written by |
|---|---|---|---|
| Table 1 (data sources), Table 2 (nominal weights) | Tables 1–2 | no file: constants of `code/final_h4m_state.py` | — |
| Table 3, the five fields at two scales | Table 10 (reformulated) | `revision/out/table3/TABLA10_REFORMULADA.csv`; inputs `eval/tabla10_predictores.csv`, `eval/captura_celda.csv`, `revision/out/S3_ablation/P1A_ABLACION.csv`, `revision/out/S5_counts_densities/RHO_CONTEOS_DENSIDADES.csv` | `revision/code/tabla10_reformulada.py` |
| Table 4, hydrographic treatment | Table 11 | `eval/tabla11_hidrologia.csv` | `code/final_evaluate.py` |
| Figure 7, Ward-4 zones, mean *F* and site-labelled cells | new | `revision/out/figure7/fig7_ward4_F_sitios.pdf` (data: `fig7_datos.csv`) | `revision/code/figura7_ward4.py` |
| Figures 1–6 | Figures 1–6 | as in [outputs.md](outputs.md) | — |
| FETE values quoted in Section 4.2 | Table 8 | `eval/tabla8_fete.csv` | `code/final_evaluate.py` |
| Ablation values quoted in Section 4.4 | new | `revision/out/S3_ablation/P1A_ABLACION_RESUMEN.csv` | `revision/code/s1_ablacion_campos.py` |

## Supplementary Material

| Element | Previous number | File | Written by |
|---|---|---|---|
| Table S1, weight sensitivity | Table 3 | `eval/tabla3_pesos.csv` | `code/final_evaluate.py` |
| Table S2, sensitivity to *T* | Table 4 | `eval/tabla4_T.csv` | `code/final_evaluate.py` |
| Table S3, attractor quantile | Table 5 | `eval/tabla5_cuantil.csv` | `code/final_evaluate.py` |
| Table S4, ablation of *L*; Figure S1 | new | `revision/out/S3_ablation/S1_ABLACION_TABLA.csv`, `S1_ABLACION_FIGURA.pdf` | `revision/code/s1_ablacion_campos.py`, `s1_tabla_figura.py` |
| Table S5, the eight unverified reaches; Figure S2 | new | `revision/out/S4_reaches/S2_OCHO_TRAMOS.csv`, `S2_MAPA.pdf`, `S2_ocho_tramos.gpkg` | `revision/code/s2_ocho_tramos.py` (needs the HydroRIVERS shapefile of South America and the review form, not distributed) |
| Table S6, partitions | Table 12 | `eval/tabla12_particiones.csv`, `eval/kmeans100_resumen.csv` | `code/final_evaluate.py` |
| Table S7, geographic criteria of the partitions | new | `revision/out/S5_ward4/S3_COMPARACION_PARTICIONES.csv`, `revision/out/S5_counts_densities/ZONAS_AREA_SITIOS.csv` | `revision/code/s3_ward4.py`, `zonal_conteos_densidades.py` |
| Table S8, the 20 Ward-4 zones | new | `revision/out/figure7/fig7_datos.csv`, `revision/out/S5_ward4/S3_WARD4_ZONAS.csv` | `revision/code/figura7_ward4.py`, `s3_ward4.py` |
| Table S9, counts, densities and the check without the zone with most site-labelled cells | new | `revision/out/S5_counts_densities/RHO_RESUMEN.csv`, `RHO_SIN_ZONA_DOMINANTE.csv` | `revision/code/zonal_conteos_densidades.py` |
| Table S10, Kolmogorov–Smirnov | Table 6 | `eval/tabla6_ks.csv` | `code/final_evaluate.py` |
| Table S11, point-process models | Table 7 | `ppm/ppm_H4m_multi_scale_aic.csv` | `code/final_ppm.R` (input *restricted*) |
| Tables S12–S13, FETE | Tables 8–9 | `eval/tabla8_fete.csv`, `eval/tabla9_fete_presupuestos.csv` | `code/final_evaluate.py` |
| Table S14, routing budget and run time | new | `revision/out/S6_fete_budget/S4_FETE_PRESUPUESTO.csv` | `revision/code/s4_fete_s6_mincul.py` |
| Figure S3, K-means zones | Figure 7 | layers in `qgis/` (drawn in QGIS) | `code/final_exports.py` |
| Figures S4–S5 | Figures 8–9 | `fig_etapa3/fig8_captura_zonas_ref_uniforme.pdf`, `fig/fig9_correlaciones.pdf` | `code/b9_etapa3_fig8.py`, `code/final_figures.py` |

The scripts of `revision/code/` read the stored runs and evaluation tables of `out/` and the site-labelled cells of the
analysis grid, so they need the restricted data; their outputs contain only aggregated values (by zone, field or realisation).
