# Outputs and where each result of the article comes from

All paths are relative to `out/final/`. Every CSV in `eval/` except the three routing files is written by
[`code/final_evaluate.py`](../code/final_evaluate.py) from the stored runs; numbers are rounded only when they are
written into the article. Spanish words in names are translated in [glossary.md](glossary.md).

## Article → file

| Article | File | Written by |
|---|---|---|
| Table 1 (data sources) and Table 2 (nominal weights) | no file: Table 2 repeats the constants of `code/final_h4m_state.py` | — |
| Table 3, weight sensitivity | `eval/tabla3_pesos.csv` (per seed: `eval/sensibilidades_por_semilla.csv`) | `final_evaluate.py` |
| Table 4, sensitivity to *T* | `eval/tabla4_T.csv` | `final_evaluate.py` |
| Table 5, attractor quantile | `eval/tabla5_cuantil.csv` | `final_evaluate.py` |
| Table 6, Kolmogorov–Smirnov | `eval/tabla6_ks.csv` | `final_evaluate.py` |
| Table 7, point-process models | `ppm/ppm_H4m_multi_scale_aic.csv` | `final_ppm.R` |
| Table 8, FETE baseline | `eval/tabla8_fete.csv` | `final_evaluate.py` |
| Table 9, FETE budgets | `eval/tabla9_fete_presupuestos.csv` | `final_evaluate.py` |
| Table 10, simple predictors | `eval/tabla10_predictores.csv` | `final_evaluate.py` |
| Table 11, hydrographic treatment | `eval/tabla11_hidrologia.csv` (per seed: `eval/hidrologia_por_semilla.csv`) | `final_evaluate.py` |
| Table 12, partitions | `eval/tabla12_particiones.csv`, `eval/kmeans100_resumen.csv` | `final_evaluate.py` |
| Figure 8, capture by zones | `fig_etapa3/fig8_captura_zonas_ref_uniforme.pdf` (data: `eval/zonas_prefijo.csv`, `eval/zonal_P2.csv`) | `b9_etapa3_fig8.py` |
| Figure 9, zonal correlations | `fig/fig9_correlaciones.pdf` (data: `eval/zonal_P2.csv`) | `final_figures.py` |
| Figures 4–7, maps | layers in `qgis/`; the maps were drawn in QGIS | `final_exports.py` |
| Figures 1–3 | study-area and elevation maps drawn in QGIS from `data/` and other sources, and a diagram; not produced by code | — |
| Numbers quoted in the text | `eval/cifras_texto.csv` and the other files of `eval/`, `ppm/` and `informe/` | see below |
| Run time of one realisation | `runtime_s` in `runs/S2/*.json` (median of the 30 canonical runs) | `final_runs.py` |

Figure 8 differs from `fig/fig8_captura_zonas.pdf` only in the legend label of the diagonal ("Uniform-density
reference" instead of "Random expectation"); `fig_etapa3/fig8_comprobacion.json` records that the rest of the image is
pixel-identical.

## `eval/` — evaluation tables

| File | Content | Main columns |
|---|---|---|
| `tabla3_pesos.csv` | 17 weight configurations (nominal and ±20 % of each weight; e.g. `a_w_-20`) | `config`, `attractors`, `pairs`, `n_seeds`, `auc_mean`, `auc_sd`, `rho20_mean`, `rho20_sd` |
| `tabla4_T.csv` | 7 configurations of *T* (`r`, `s`, `w`) | as above |
| `tabla5_cuantil.csv` | 4 attractor quantiles | as above |
| `sensibilidades_por_semilla.csv` | every seed of S3–S5 | `config`, `seed`, `attractors`, `pairs`, `auc`, `rho20`, `bloque` (block) |
| `tabla6_ks.csv` | KS tests of *E′* and *L*, all 84 site-labelled cells and the 20 outside the dominant zone | `site_set`, `field`, `D`, `p`, `method_used`, `n_sites`, `n_background` |
| `tabla8_fete.csv` | corridor model against FETE on *c′* and on 1 − *E′* (90 pairs, 30 seeds) | `method`, `n_pairs`, `n_seeds`, `auc_mean`, `auc_sd`, `lcc_mean`, `lcc_sd`, `m2_mean`, `m2_sd`, `all_paths_found` |
| `tabla9_fete_presupuestos.csv` | FETE on *c′* with 50–1,000 pairs | as above |
| `tabla10_predictores.csv` | AUC and ρ₂₀ of 1 − *h*, *T*, *E*, *E′*, *F* and *L* | `predictor`, `auc`, `rho20`, `p20` |
| `tabla11_hidrologia.csv` | four hydrographic treatments, mean of seeds 42–44 | `scenario`, `auc_mean`, `rho20_mean`, `capture10_mean`, …; columns marked `(B5)` repeat the decision metrics of Phase B5 |
| `hidrologia_por_semilla.csv` | the same per seed | `scenario`, `seed`, `auc`, `rho20`, `capture10_pct` |
| `tabla12_particiones.csv` | K-means, Ward (4 and 8 neighbours) and SKATER with 20, 30 and 40 zones | `algorithm`, `n`, `rho`, `p`, `rho_loo`, `p_loo`, zone-shape and size columns |
| `kmeans100_resumen.csv` | ρ over 100 K-means initialisations | `n`, `n_inits`, `rho_median`, `rho_p05`, `rho_p95`, `rho_loo_median`, `rho_loo_p05` |
| `zonal_P2.csv` | zonal correlations of the reference partition, n = 20, 30, 40 | `rho`, `p`, `rho_max`, `p_max`, `rho_loo`, `p_loo`, `dominant_zone`, `dominant_sites` |
| `zonas_prefijo.csv` | zones ranked by mean *L* with cumulative area and sites | `zone`, `L_mean`, `n_sites`, `n_cells`, `rank`, `cum_area_pct`, `cum_sites_pct`, `lift`, `n` |
| `captura_zonas.csv` | sites captured by complete zones within area budgets | `partition`, `ranking_field`, `budget_pct`, `zones_used`, `area_consumed_pct`, `sites_pct` |
| `captura_celda.csv` | cell capture at budgets of 5–30 % | `predictor`, `field`, `budget_pct`, `n_cells`, `sites`, `sites_pct`, `tie_min_pct`, `tie_max_pct` |
| `captura_celda_resumen.csv` | cell capture of *L* over the 30 seeds of S2 | `budget_pct`, `sites_pct_mean`, `sites_pct_sd`, `sites_pct_min`, `sites_pct_max` |
| `solapamiento.csv` | buffered overlap of the top-*K* cells with the site-labelled cells | `b` (buffer, cells), `K`, `K_mult`, `overlap` |
| `descomposicion_L.csv` | contribution of each term of *L* under three scalings | `scaling`, `term`, `mean`, `max`, `var`, `share_of_mean_sum` |
| `normalizacion_E5.csv` | *L* under three normalisations of *F* | `surface`, `AUC`, `rho20_P2`, `spearman_vs_N1`, … |
| `pendiente_E4_H4m.csv` | slope transformations (block S9) against the final configuration | `variant`, `attractors`, `attractors_shared_with_S1`, `auc`, `rho20_P2`, Spearman correlations |
| `barrido_S6.csv` | ε–κ sweep: mean and SD over 10 seeds for each of the 121 points | `eps`, `kappa`, `n`, `auc_media`, `auc_de`, `rho20_media`, `rho20_de`, `lcc_media`, `m2_media` |
| `cifras_texto.csv` | numbers of the text that are not in a table, with an identifier and an English description | `id`, `descripcion`, `valor` |
| `pares_rutas.csv` | the 90 pairs of the canonical run: length, sinuosity, cost, share of traffic | written by `final_diagnostics.py` |
| `resolucion_E2A_pares.csv`, `resolucion_E2A_resumen.csv` | the 90 pairs routed on *c′* at 250 m and 50 m with 8 and 16 directions (A1 = 250 m/8, A2 = 50 m/8, A3 = 50 m/16, A4 = 250 m/16) | written by `final_diagnostics.py` |
| `indice.json` | sha256 of the procedures, of `final_evaluate.py` and of every file it wrote, and the signature of code and data | written by `final_evaluate.py` |

## `informe/` — diagnostics behind numbers of the text

| File | Content | Written by |
|---|---|---|
| `diagnostico_E_atractores_sitios.csv` | affinity at attractors and site-labelled cells and their distance to the nearest reach, previous version (v2: *E*, 125 m buffer) and final (*E′*) | `b9_diag_E.py` |
| `rangos_campos_figura4.csv` | minimum and maximum of the fields mapped in Figure 4 | `b9_diag_E.py` |
| `atractores_sensibilidades.csv` | attractor set of every configuration of Tables 3–5 compared with the final set; the rebuilt state matches the `state_hash` of its run | `b9_diag_atractores.py` |
| `conteo_atractores_8_lineas.csv` | distance of each attractor to the eight unverified reaches, and the attractors added when they are removed | `b9_diag_atractores.py` |
| `diagnostico_L_vs_T.csv` | Spearman correlation between the zone means of *L* and *T* (and *F*) for 20, 30 and 40 zones, and between *L* and *T* over cells | `diag_L_vs_T.py` |

## `ppm/` — point-process models

`fields_H4m_for_ppm.csv` (restricted) is the input (cell centres `X`, `Y`; `F`, `L` and `E` = *E′* of seed 42; `is_site`), written by
`final_export_ppm.py`. `ppm_H4m_multi_scale_aic.csv` (AIC, deviance and slope by covariate and radius) and
`ppm_H4m_summary.txt` are written by `final_ppm.R`, copied here with its log `ppm_run.log`. The first lines of
`ppm_H4m_summary.txt` still say "MODEL V2" because they were copied from the script of the previous version; the numbers
are those of the final configuration.

## `runs/` — stored runs (blocks S1–S5, S7–S9, CTRL)

One JSON per run with the job specification (`spec`), the signature of code and data (`firma`), the hash of the model
state (`state_hash`), the numbers of attractors and pairs, the run metrics, the run time and the sha256 of the field
file. Runs with fields also have a `.npz` with `F` and `L` on the 93 × 144 lattice; the two runs of S1 also have
`_paths.pkl` with every routed path. `runs_S1-S9.log` is the log of the final computation.

## `S6/` — ε–κ sweep

The 1,210 runs were computed on a second machine with
[`out/b6/ryzen/barrido_H4m/final_sweep_eps_kappa.py`](../out/b6/ryzen/barrido_H4m/final_sweep_eps_kappa.py) and the same
code and data. `resultados_barrido_H4m/` has one JSON per run (`e<ε>_k<κ>_s<seed>.json`); `resumen_barrido_H4m.csv` is
the summary computed on that machine; `prueba_OK.json` records its environment test; `MANIFEST.sha256` and
`para_enviar_barrido_H4m.zip` are the manifest and the archive received from it; `verificacion_cruzada_2_trabajos.csv`
compares two runs repeated on the main machine (largest difference 4.3 × 10⁻¹⁸).

## `qgis/`, `fig/`, `previa/` and records

- `qgis/`: GeoTIFF and GeoPackage layers (EPSG:32718) of the input fields, *E′*, *c′*, *F* and *L* of seed 42, the 25
  attractors, the 84 site-labelled cells (cell centres; `sitios_84.gpkg`, restricted), the routes of S1 and the zones of P2.
- `fig/`: Figures 8 and 9 as first drawn from the evaluation CSVs; `fig_etapa3/`: Figure 8 as published.
- `previa/`: tests run before the final computation (reproduction of the previous version, KS provenance, environment
  test of the sweep machine).
- `PROCEDIMIENTOS.md`: evaluation procedures written before the final values (in Spanish; summary in
  [model.md](model.md#evaluation)). `MANIFEST_congelado.sha256`: sha256 of the code, data and procedures frozen before
  the final runs. `REGISTRO_CAMBIOS_POST_CONGELACION.md`: every change made after the freeze (in Spanish; summary in
  [provenance.md](provenance.md)). `diagnostics.log`: log of `final_diagnostics.py`.
