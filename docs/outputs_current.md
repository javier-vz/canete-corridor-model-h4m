# Outputs corresponding to the revised JCAA-80 manuscript

**Source of truth:** the final manuscript supplied to this staging review:
`main_v12_2_para_revision_Dina(6).tex` and its 26-page PDF.
The data and scripts below are frozen outputs of H4m, not newly computed results.

| Manuscript item | Repository path / explanation |
|---|---|
| Table 1 (data) | `data/README.md` and `code/final_h4m_state.py` |
| Table 2 (weights) | `code/final_h4m_state.py` and `docs/model.md` |
| Table 3 (AUC, capture, F/L/T/E′, Ward-4) | `revision/out/table3/TABLA10_REFORMULADA.csv` and `revision/code/tabla10_reformulada.py` |
| Table 4 (hydrological sensitivity) | `out/final/eval/tabla11_hidrologia.csv` |
| Figure 7 (Ward-4 and simulated traffic) | `revision/out/figure7/fig7_ward4_F_sitios.pdf` and `fig7_datos.csv` (aggregate zone counts) |
| Figure 5 / Figure 6 (84 marked cells) | **NOT redistributed here:** point-labelled figures require archaeological data release approval |
| Figure 1 / Figure 2 | maps in the manuscript; publication of external imagery subject to source terms |
| Figure 3 / Figure 4 | manuscript workflow and derived environmental fields; source code in `code/` |
| Traffic F (seed 42) | `out/final/qgis/F_canon_s42.tif` and frozen model runs under `out/final/runs/` |
| Simpler comparator fields | `out/final/qgis/T.tif`, `E_prime.tif`, `h_norm.tif`, etc. |
| FETE diagnostics | `out/final/eval/tabla8_fete.csv`, Supplementary S6 and `revision/out/S6_fete_budget/` |
| Model ablations | `revision/out/S3_ablation/` and Supplementary S3 |
| Regionalisation | `out/e3/partitions.npz`, `revision/out/S5_counts_densities/` and Supplementary S5 |

`docs/outputs_v10.md` and `docs/outputs.md` document previous editorial numbering; use this
file for the revised manuscript. **The full archaeological reference set is withheld** and its
evaluation cannot be completely reproduced from public data alone.
