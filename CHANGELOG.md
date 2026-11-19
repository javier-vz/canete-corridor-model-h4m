## Release candidate v2.2 (not published; 19 November 2026)

- Update public-facing documentation to match the revised paper: F is the primary simulated
  traffic output, L a secondary composite comparator, and evaluations use 84 *labelled cells*.
- Document exclusion of archaeological coordinates, site-distance columns and derived layers.
- Flag locations still displayed in the manuscript Figures 5/6 and Supplementary Figure S3.
- Link a page-by-page audit of the existing S1–S6 supplementary package.
- Preserve frozen H4m code, inputs and outputs byte-for-byte (no new simulations).
- Document why a clean public candidate does not by itself purge historical GitHub releases.

---

# Changelog

## v2.1.0 (candidate, not published) — simplified version of the revised article

- Restricted data: the three site columns of the analysis grid and six files with the site-labelled cells are left out of
  the public repository (`RESTRICTED_MANIFEST.sha256`); the grid is distributed without those columns
  (`data/grilla_con_W_hidrica_publica.csv`). `tools/restore_restricted.py` rebuilds the original files byte for byte from
  an authorised bundle, and `verify_repository.py` then runs unchanged.
- `verify_public.py` checks the public files; the GitHub Actions workflow runs it.
- `revision/`: scripts and outputs added for the simplified version (Table 3, Figure 7 and Supplementary Tables S4–S9 and
  S14, Figures S1–S2), with paths relative to the repository; `docs/outputs_v10.md` maps the new numbering.
- No frozen file of `code/`, `data/` or `out/` was modified, and no run was repeated.
- Data documentation: the description of the archaeological reference set was corrected (to be completed by the authors).
  The evaluation set is described as 84 site-labelled cells (cells whose centre lies within about 150 m of a survey
  record; at least 57 distinct records), not as 84 sites, and the restricted files are marked where the documentation
  mentions them.

## v2.0.0 — final configuration of the revised article

- The repository now holds the final configuration of the revised article: `code/` (scripts of the final computation),
  `data/` (inputs), `out/` (stored runs, evaluation tables, diagnostics, figures, map layers and records) and `docs/`.
- `verify_repository.py` checks every file and recomputes the evaluation tables and diagnostics from the stored runs,
  without new simulations; a GitHub Actions workflow runs it.
- New scripts `code/check_grid_columns.py` (checks of the grid columns) and `code/diag_L_vs_T.py` (recomputes a stored
  diagnostic). No stored result was modified.
- `data/grilla_con_W_hidrica.csv`: same values as in v1.0.0, now stored with its original CRLF line endings, on which
  the checksums depend. `.gitattributes` disables end-of-line conversion for every file.
- Data documentation corrected: elevation and slope are those of the DEM pixel at the cell centre, not cell averages.
- The v1 layout (`src/`, `scripts/`, `notebooks/`, `original_scripts/`, `docs/`) was removed from the main branch; it is
  kept in the tags v1.0.0 and v1.0.1 ([docs/changes_from_v1.md](docs/changes_from_v1.md)).
- `CITATION.cff`: version, repository URL and the author order of the article in the preferred citation.

## v1.0.1

Fixes to the kernel interface (`return_fields`), the command-line scripts, the background of the leave-dominant-out KS
test and the seeds of the FETE (1 − E) runs; ablation scripts rerun with the K-means seed of the paper
(`random_state=42`).

## v1.0.0

Code, grid and DEM of the originally submitted version.
