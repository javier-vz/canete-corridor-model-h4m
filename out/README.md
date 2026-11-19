# Stored outputs

| Folder or file | Content |
|---|---|
| `final/` | everything the article reports: stored runs, evaluation tables, diagnostics, point-process results, figures, map layers and the records of the final computation ([docs/outputs.md](../docs/outputs.md)) |
| `b4/runs/`, `b5/runs/` | runs of the hydrographic sensitivity (Table 11), seeds 42–44, and the stored final state `b4/runs/def_H4m.pkl` |
| `b4/t1/` | match between HydroRIVERS lines and other line files, used to identify the eight unverified reaches |
| `b5/t5/` | decision metrics and attractor changes of the hydrographic scenarios |
| `b6/ryzen/barrido_H4m/` | package sent to the second machine for the ε–κ sweep: its script, copies of the kernel and of the state construction, data, reference run and manifest; `b6/ryzen/JCAA80_barrido_H4m_Ryzen.zip` is the archive that was sent (restricted, like the copy of the grid in its `data/` folder) |
| `b6/design_counts_H4m.csv` | expected attractor and pair counts of every sensitivity configuration (`code/b6_design_counts.py`) |
| `b6/fingerprint_H4m.json` | fingerprint of the final state (`code/test_final_h4m_state.py`) |
| `e3/` | Ward and SKATER partitions with 20, 30 and 40 zones |
| `e4/scenario_fields.pkl` | slope fields of the slope transformations (block S9) |
| `hydro/` | classified model lines and main-channel cells |
| `kmeans_labels_100.npy`, `rowmajor_iy_ix.npy`, `rowmajor_sites.npy` | 100 K-means partitions, the row-major cell order and the site-labelled cells in that order (`rowmajor_sites.npy`, restricted) |
| `river_cell_mask.npy`, `dist_river_grid.npy` | river cells and distance to the nearest reach on the 93 × 144 lattice |
| `state_v2_repro.pkl` | model state of the previous version (v2); restricted |

The files outside `final/` were prepared in earlier phases of the revision and are read as frozen inputs
([docs/provenance.md](../docs/provenance.md)). The `.pkl` files are Python pickles written by the authors; load them only
from this repository.
