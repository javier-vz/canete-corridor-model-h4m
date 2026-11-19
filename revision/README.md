# Revision of the simplified version (v10)

Scripts and outputs added for the simplified version of the article and its Supplementary Material. They do not run
the model and do not change any frozen file of `code/`, `data/` or `out/`: they read the stored runs and evaluation
tables and compute aggregated values. [docs/outputs_v10.md](../docs/outputs_v10.md) maps each table and figure of the
article to these files.

| Script (`code/`) | Output (`out/`) | Content |
|---|---|---|
| `s1_ablacion_campos.py` | `S3_ablation/P1A_*.csv` | exploratory ablation of *L*, recomputed from the stored fields of the 31 canonical realisations (seed 42 and seeds 2000–2029); controls against the frozen tables |
| `s1_tabla_figura.py` | `S3_ablation/S1_ABLACION_TABLA.csv`, `S1_tabla_filas.tex`, `S1_ABLACION_FIGURA.pdf` | Table S4 and Figure S1 |
| `zonal_conteos_densidades.py` | `S5_counts_densities/*.csv` | zonal correlations with the counts and with the densities of site-labelled cells, zone areas, and the correlation without the zone with most site-labelled cells (Table S9) |
| `tabla10_reformulada.py` | `table3/*.csv`, `tabla10_filas.tex` | Table 3 of the main text |
| `figura7_ward4.py` | `figure7/` | Figure 7 of the main text and its data |
| `s3_ward4.py` | `S5_ward4/*.csv` | geographic criteria of the 20-zone partitions and the Ward-4 zones (Tables S7–S8) |
| `s4_fete_presupuesto.py` | `S6_fete_budget/S4_FETE_PRESUPUESTO.csv` | routing budget and run time of the corridor model and the FETE baselines (Table S14) |
| `s2_ocho_tramos.py` | `S4_reaches/` | the eight unverified HydroRIVERS reaches (Table S5, Figure S2); run with the environment variables `HYDRORIVERS_SA_SHP` and `REACH_REVIEW_FORM` |

## Running them

All scripts except `s2_ocho_tramos.py` run from the repository root with the restricted data restored
(`tools/restore_restricted.py`), in this order:

```bash
python revision/code/s1_ablacion_campos.py
python revision/code/s1_tabla_figura.py
python revision/code/zonal_conteos_densidades.py
python revision/code/tabla10_reformulada.py
python revision/code/figura7_ward4.py
python revision/code/s3_ward4.py
python revision/code/s4_fete_presupuesto.py
```

They were checked by running them in a restored copy of this repository and comparing their CSV outputs with the files
used for the manuscript (`REVISION_CHECK.csv`). PDF and PNG files can differ in their metadata. `s2_ocho_tramos.py` also needs the
HydroRIVERS v1.0 shapefile of South America (HydroSHEDS, not redistributed) and the imagery review form of the eight
reaches (not distributed; its content is summarised in `S4_reaches/S2_OCHO_TRAMOS.csv`).

The files of `out/` contain values aggregated by zone, field, reach or realisation; none gives the location of a site.
The reach layer `S4_reaches/S2_ocho_tramos.gpkg` derives from HydroRIVERS and keeps the HydroSHEDS licence
([THIRD_PARTY.md](../THIRD_PARTY.md)).
