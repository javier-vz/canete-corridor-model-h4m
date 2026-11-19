# Reproduce and verify

## Environment

- Python 3.13 with the versions pinned in [`requirements.txt`](../requirements.txt). The final computation, the
  evaluation and the verification of this release were run on Linux x86-64.
- The ε–κ sweep (block S6) was run on a Windows 11 machine with Python 3.12.13 and the same code and data
  (`out/final/S6/prueba_OK.json`).
- R 4.3.3 with spatstat 3.0.7, only for the point-process models (`code/final_ppm.R`).

Run every command from the root of the repository.

## Verification

```bash
python verify_repository.py            # about 3 minutes
python verify_repository.py --keep     # keeps the temporary copy and prints its path
```

The script copies `code/`, `data/` and `out/` to a temporary folder, runs the checks there and never writes to the
repository. It does not run the model.

| Step | Check |
|---|---|
| 1 | sha256 of every file listed in `MANIFEST.sha256` (all files of `code/`, `data/` and `out/`) |
| 2 | `code/final_freeze.py --verificar`: the code, data and procedures match the freeze manifest written before the final runs, except `code/test_final_pipeline.py`, whose post-freeze change is documented (Change 1 in [provenance.md](provenance.md)) |
| 3 | `code/check_grid_columns.py`: the columns of the analysis grid against the DEM and the HydroRIVERS lines ([data/README.md](../data/README.md)) |
| 4 | `code/test_final_h4m_state.py`: the final state (attractors, pairs, *E′*, *c′*, *T* and every edge weight) is rebuilt and compared with the stored states |
| 5 | `code/final_runs.py --verificar`: each of the 309 stored runs has the expected specification, code-and-data signature, metrics and field files |
| 6 | `code/final_evaluate.py` recomputes the 23 evaluation CSVs and `indice.json` from the stored runs; each must equal the stored file |
| 7 | `code/b9_diag_E.py`, `code/b9_diag_atractores.py`, `code/diag_L_vs_T.py` and `code/b6_design_counts.py` recompute the diagnostics and the expected attractor and pair counts; each output must equal the stored file |
| 8 | no other file of the temporary copy was created or changed |

**Equality.** On the platform used for this release every recomputed file is byte-identical. On another processor or
operating system floating-point results can differ in the last digits: on the sweep machine *E′* differed by at most
1.1 × 10⁻¹⁶ and *L* by 2.2 × 10⁻¹⁶, with identical AUC and ρ₂₀. The script therefore accepts a recomputed table when its
text cells are identical and its numeric cells differ by at most 10⁻⁹ in relative terms, and it reports which files
were byte-identical and the largest difference of the others. In step 4 a failure of an exact comparison of
floating-point fields is followed by the same comparison with a tolerance of 10⁻¹²; attractors and pairs must always be
identical.

**Not covered.** The point-process models (R), the figures and the map layers are checked by sha256 only; Figure 8 can
be redrawn with `code/b9_etapa3_fig8.py`, which compares it pixel by pixel with the stored figure.

## Re-run the computation

The stored runs are what the article reports. To recompute them, work in a separate copy of the repository so that the
stored files stay as they are.

```bash
# 309 runs of S1-S5, S7-S9 and CTRL (about 3 CPU-hours; one canonical run takes about 1 minute)
python code/final_runs.py --bloques S1 S2 S3 S4 S5 S7 S8 S9 CTRL --workers 4 --salida out/rerun/runs

# epsilon-kappa sweep, block S6 (1,210 runs, about 117 CPU-hours), from its own folder
cd out/b6/ryzen/barrido_H4m
python final_sweep_eps_kappa.py --prueba     # environment test: reproduces the reference state and run
python final_sweep_eps_kappa.py --lanzar     # runs missing or invalid jobs; --workers N to limit processes
python final_sweep_eps_kappa.py --resumen    # summary of the sweep
cd ../../../..
```

`final_runs.py` only computes runs that are missing or invalid in the output folder, so it can be stopped and
restarted. To recompute the tables from new runs, replace `out/final/runs/` in the copy and run the scripts below.

```bash
python code/final_diagnostics.py      # routing at 250 m and 50 m and route statistics (eval/pares_rutas.csv, eval/resolucion_E2A_*.csv)
python code/final_evaluate.py         # evaluation tables (out/final/eval/)
python code/final_figures.py          # Figures 8 and 9 (out/final/fig/)
python code/b9_etapa3_fig8.py         # Figure 8 as published (out/final/fig_etapa3/)
python code/final_exports.py          # map layers of Figures 4-7 (out/final/qgis/)
python code/final_export_ppm.py       # input of the point-process models
Rscript code/final_ppm.R              # point-process models (out/final/ppm/)
python code/b9_diag_E.py              # diagnostics of out/final/informe/
python code/b9_diag_atractores.py
python code/diag_L_vs_T.py
```

The hydrographic sensitivity runs of Table 11 (`out/b4/runs/`, `out/b5/runs/`) and the stored partitions are inputs
prepared in earlier phases ([provenance.md](provenance.md)); they are not recomputed by these commands.

## Pre-tests

`code/test_final_pipeline.py` holds the pre-tests run before the final computation (stops at the first failure);
`--ctrl` compares the CTRL runs with the published values of the previous version. Their results are in
`out/final/previa/`.
