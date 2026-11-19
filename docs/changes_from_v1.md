# Changes from v1.0.x

Versions v1.0.0 and v1.0.1 correspond to the earlier implementation used for the originally submitted article. Their historical code is not included in this repository.

## Model

- **Hydrological rules.** v1 had none. The final configuration excludes the 1,020 river cells and cells within 125 m of
  a reach from attractor selection, reduces the affinity next to reaches ($E'$) and adds a cost on river cells ($c'$)
  ([model.md](model.md)). An intermediate version (v2) excluded only cells within 125 m of a reach; it is used here as a
  reference and for the CTRL runs.
- **Kernel.** `code/ca_kernel_v2.py` is the kernel of v1.0.1 (`src/ca_kernel.py`) with an optional attractor-exclusion
  mask. The hydrological rules are applied on top of it by `code/final_h4m_state.py`.
- **Evaluation.** The procedures were fixed in writing before the final computation; the evaluation adds cell capture at
  fixed budgets, exact Kolmogorov–Smirnov tests, contiguity-constrained partitions, a hydrographic sensitivity, slope
  transformations, a routing comparison at 50 m and point-process models.

## Data

- `data/grilla_con_W_hidrica.csv` and `data/DEM_vallebajo.tif` are the files of v1.0.0. The grid is now stored with its
  original Windows line endings (CRLF), which the checksums of the final computation depend on; v1.0.x stored it with
  LF line endings, so Git shows the whole file as changed although no value changed.
- The data documentation of v1.0.x said that each cell takes the mean elevation of its pixels and a slope from a 3 × 3
  derivative of the aggregated grid. The grid actually holds the elevation of the single DEM pixel at the cell centre
  and the central-difference slope of that pixel; [data/README.md](../data/README.md) describes and checks this. The
  notebook `02_data_preparation.ipynb` of v1.0.x describes the same averaging and is not included.
- New: the list of river cells, the HydroRIVERS clip used for them and the frozen inputs of `out/`.

## Repository

- The v1 layout (`src/`, `scripts/`, `notebooks/`, `original_scripts/`, `docs/ARCHITECTURE.md`, `docs/FAQ.md`,
  `docs/REPRODUCE.md`) is not in v2.0.0; it is kept in the tags above.
- New: `code/` with the scripts of the final computation, `out/` with every stored output, `verify_repository.py`,
  `MANIFEST.sha256`, a GitHub Actions workflow that runs the verification, and these documents.
- `CITATION.cff`: version 2.0.0, the repository URL and the author order of the article in the preferred citation.
