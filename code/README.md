# Scripts

Run every script from the root of the repository (`python code/<script>.py`). Paths are relative to the root; set the
environment variable `JCAA_ROOT` to use another root. Names are those of the frozen computation
([docs/glossary.md](../docs/glossary.md)).

## Model and final configuration

| Script | Role |
|---|---|
| `ca_kernel_v2.py` | kernel: fields, attractors, pairs, routing graph, one realisation, run metrics |
| `final_h4m_state.py` | final state: the hydrological rules on top of the kernel, and the designs of the sensitivity configurations |
| `final_common.py` | shared helpers of the final computation: state cache, partitions, zonal statistics, signature of code and data |
| `canonical_partition.py` | reference K-means partition (P2) |
| `routing.py` | least-cost routing on a regular lattice with 8 or 16 directions (routing comparison) |
| `river_geometry.py` | HydroRIVERS lines of the model window and river-cell masks |
| `h_scenarios.py`, `scenario_eval.py` | states of the hydrographic scenarios of the revision; `test_final_h4m_state.py` compares the final state with them |
| `paths.py` | root of the repository |

## Final computation

| Script | Role | Output |
|---|---|---|
| `final_runs.py` | runs blocks S1–S5, S7–S9 and CTRL; `--verificar` validates the stored runs | `out/final/runs/` |
| `../out/b6/ryzen/barrido_H4m/final_sweep_eps_kappa.py` | ε–κ sweep (block S6), run on a second machine | `out/final/S6/` |
| `final_evaluate.py` | every evaluation table, from the stored runs | `out/final/eval/` |
| `final_diagnostics.py` | routing comparison at 250 m and 50 m and route statistics | `out/final/eval/pares_rutas.csv`, `resolucion_E2A_*.csv` |
| `final_figures.py` | Figures 8 and 9 | `out/final/fig/` |
| `final_exports.py` | map layers of Figures 4–7 | `out/final/qgis/` |
| `final_export_ppm.py`, `final_ppm.R` | input and fit of the point-process models | `out/final/ppm/` |
| `final_ks_provenance.py` | reproduces the historical KS implementation, only to document where earlier values came from | `out/final/previa/ks_procedencia.csv` |
| `final_freeze.py` | writes the freeze manifest; `--verificar` checks it | `out/final/MANIFEST_congelado.sha256` |

## Diagnostics and figures added after the freeze

| Script | Output |
|---|---|
| `b6_design_counts.py` | `out/b6/design_counts_H4m.csv` |
| `b9_diag_E.py` | `out/final/informe/diagnostico_E_atractores_sitios.csv`, `rangos_campos_figura4.csv` |
| `b9_diag_atractores.py` | `out/final/informe/atractores_sensibilidades.csv`, `conteo_atractores_8_lineas.csv` |
| `diag_L_vs_T.py` (added in v2.0.0) | `out/final/informe/diagnostico_L_vs_T.csv` |
| `b9_etapa3_fig8.py` | `out/final/fig_etapa3/` (Figure 8 as published) |

## Tests and checks

| Script | Role |
|---|---|
| `test_final_h4m_state.py` | rebuilds the final state and compares it with the stored states (no model runs) |
| `test_final_pipeline.py` | pre-tests run before the final computation |
| `check_grid_columns.py` (added in v2.0.0) | checks the columns of the analysis grid against the DEM and the HydroRIVERS lines |

`../verify_repository.py` runs the checks in order ([docs/reproduce.md](../docs/reproduce.md)).
