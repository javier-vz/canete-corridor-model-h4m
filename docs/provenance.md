# Provenance of the final results

This page summarises in English the records kept in `out/final/` (in Spanish). Dates are 2026, times are Lima time
(UTC−5).

## Order of the final computation

1. **Procedures fixed first.** `out/final/PROCEDIMIENTOS.md` (5 November) defines every evaluation procedure before any
   final value existed. Its sha256 is recorded in `out/final/eval/indice.json`.
2. **Freeze.** `out/final/MANIFEST_congelado.sha256` (6 November, 09:41) records the sha256 of the 54 files of code,
   data and procedures used by the final computation, written by `code/final_freeze.py` after the pre-tests passed.
3. **Signature in every run.** Each run JSON records the signature of the code and data that produce runs
   (`final_common.py`, `final_runs.py`, `final_h4m_state.py`, `ca_kernel_v2.py`, `canonical_partition.py`, the grid, the
   river cells and the slope variants); `final_runs.py --verificar` rejects a run whose signature differs.
4. **Runs.** Blocks S1–S5, S7–S9 and CTRL were run on the main machine (309 of 309 valid on 6 November, 10:57), and block
   S6 on a second machine.
5. **Evaluation.** `code/final_evaluate.py` computed every table from the stored runs.

## Changes after the freeze

`out/final/REGISTRO_CAMBIOS_POST_CONGELACION.md` lists every change made after the freeze. The ones that concern files of
this repository are:

- **Change 1 (6 November, 11:10), `code/test_final_pipeline.py` only.** Its temporary folders became portable and it
  reads the published input of the previous point-process analysis from `data/ref/fields_canonical_v2_for_ppm.csv` (restricted), an
  unchanged copy of the original. No file of the run signature changed; the pre-tests and the controls were repeated
  with the same results (`out/final/previa/pruebas_previas_post_cambio_2026-11-06.csv`). This is why
  `final_freeze.py --verificar` reports exactly this one file as changed.
- **Integration of S6 (6 November, 17:55).** The archive received from the second machine
  (`out/final/S6/para_enviar_barrido_H4m.zip`) contains the same manifest as the package that was sent and its 1,210
  runs are valid. Two runs chosen at random were repeated on the main machine; the largest difference was 4.3 × 10⁻¹⁸
  (`out/final/S6/verificacion_cruzada_2_trabajos.csv`). `final_evaluate.py` was run again without code changes: it added
  `barrido_S6.csv` and the X7 rows of `cifras_texto.csv`, and the 24 existing tables were reproduced byte for byte.
  Figures 8 and 9 were redrawn from the same CSVs and are pixel-identical.
- **Diagnostics for the revised text (6 November, 18:50).** New scripts that only read frozen states and outputs:
  `code/b9_diag_E.py` and `code/b9_diag_atractores.py`, which write four files of `out/final/informe/`.
- **Figure 8 as published (6 November, 23:21).** `code/b9_etapa3_fig8.py` redraws Figure 8 with the frozen function
  `final_figures.fig8` and the same CSVs, changing only the legend label of the diagonal (`out/final/fig_etapa3/`).
- **Release v2.0.0.** Two scripts were added for the repository: `code/diag_L_vs_T.py`, which recomputes
  `out/final/informe/diagnostico_L_vs_T.csv` with the computation of the internal report script that first wrote it, and
  `code/check_grid_columns.py`, a read-only check of the grid. No stored file was modified.

The kernel, the final state, the runner, the evaluation, the figures of `out/final/fig/`, the map layers and the
point-process analysis were not changed after the freeze.

## Inputs prepared in earlier phases

The final computation reads the following files as frozen inputs (their sha256 are in the freeze manifest or in
`MANIFEST.sha256`). They were produced in earlier phases of the revision by scripts that are not part of this release.

| File | Content | Produced by |
|---|---|---|
| `out/state_v2_repro.pkl` | model state of the previous version (v2), used as a reference and for the lattice geometry; restricted | `repro_v2.py` |
| `out/river_cell_mask.npy`, `out/dist_river_grid.npy` | river cells and distance to the nearest reach on the lattice | `make_river_mask.py` |
| `out/hydro/hydrorivers_model_lines_classified.gpkg`, `out/hydro/river_cells_main.npy` | the 87 HydroRIVERS reaches that intersect the grid extent plus 5 km (the 81 lines of the clip and 6 more at least 4 km from every cell), with Strahler order, upstream length and a main-channel flag derived from the topology of the network, and the cells of the main channel; used by the main-channel and eight-reach scenarios of Table 11 | `hydro_classify.py` |
| `out/b4/t1/hydrorivers_official_match.csv` | line-by-line match of four line files with the records of the HydroRIVERS shapefile; only the rows of dataset C87 (the 87 features above) are used, to identify the eight unverified reaches | `hydro_trace.py` |
| `out/b4/runs/`, `out/b5/runs/` | runs of the hydrographic sensitivity (Table 11), seeds 42–44, and the state of the final configuration (`def_H4m.pkl`) | `b4_runs_with_paths.py`, `b5_exclusion_runs.py` |
| `out/b5/t5/` | decision metrics and attractor changes of the hydrographic scenarios | `b5_decision_metrics.py` |
| `out/e3/partitions.npz`, `out/e3/regionalization_comparison.csv` | Ward and SKATER partitions with 20, 30 and 40 zones | `e3_regionalization.py` |
| `out/kmeans_labels_100.npy`, `out/rowmajor_*.npy` | 100 K-means partitions (random states 0–99) and the row-major cell order | `kmeans_ensemble.py` |
| `out/e4/scenario_fields.pkl` | slope fields of the four slope transformations | `e4_slope.py` |
| `out/b6/` | sweep package sent to the second machine (`ryzen/`), expected attractor and pair counts, fingerprint of the final state | `b6_ryzen_package.py`, `code/b6_design_counts.py`, `code/test_final_h4m_state.py` |

`code/test_final_h4m_state.py` checks that the river cells, the distance grid and the final state rebuilt from the
grid agree with these files (verification step 4).
