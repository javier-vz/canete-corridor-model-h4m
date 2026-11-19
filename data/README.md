# Data

| File | Content | Terms |
|---|---|---|
| `grilla_con_W_hidrica_publica.csv` | analysis grid without its three site columns: one row per valid 250 m cell (12,826 rows, 10 columns) | CC BY 4.0; HydroSHEDS licence for the river columns |
| `grilla_con_W_hidrica.csv` | complete analysis grid (13 columns) used by the code; **restricted**, rebuilt by `tools/restore_restricted.py` | not distributed |
| `DEM_vallebajo.tif` | SRTM elevation model of the study area | public domain |
| `celdas_fluviales_1020.csv` | the 1,020 river cells *R* (`X`, `Y` of the cell centre; `main` = 1 for the 218 cells of the main channel) | HydroSHEDS licence |
| `HydroRIVERS_canete_UTM18S.gpkg` | 121 HydroRIVERS v1.0 lines around the study area (geometry only, EPSG:32718) | HydroSHEDS licence |
| `ref/fields_canonical_v2_for_ppm.csv` | input of the point-process models of the previous version (v2), used by the pre-tests to check the R environment; **restricted** (contains the site-labelled cells) | not distributed |

All coordinates are in WGS 84 / UTM zone 18S (EPSG:32718), in metres. The terms are explained in
[THIRD_PARTY.md](../THIRD_PARTY.md).

## Analysis grid (`grilla_con_W_hidrica.csv`)

| Column | Unit | Content | Read by the model |
|---|---|---|---|
| `X`, `Y` | m | centre of the cell; the centres form a regular 250 m lattice | yes |
| `distancia_sitios_nuevos` | m | **restricted**: distance from the cell centre to the nearest recorded archaeological location | no |
| `presencia_sitio_nuevo` | 0/1 | **restricted**: 1 for the 84 site-labelled cells (grid centres within approximately 150 m of recorded locations) | only for the evaluation |
| `distancia_mincul` | m | **restricted**: distance from the cell centre to the nearest site of a layer of the national heritage registry of the Ministerio de Cultura del Perú (the layer is not included) | no |
| `altitud` | m | elevation of the DEM pixel that contains the cell centre | yes, as `h_norm` |
| `pendiente` | degrees | slope of that pixel, computed on the DEM by central differences | yes, as `p_norm` |
| `dist_river_m` | m | distance from the cell centre to the nearest of the 81 HydroRIVERS lines that intersect the grid extent plus 5 km | yes, in *E′* and in the attractor exclusion |
| `dist_river_norm` | — | min–max normalisation of `dist_river_m` | no |
| `p_norm`, `h_norm` | — | min–max normalisations of `pendiente` and `altitud` | yes, as *p* and *h* |
| `W_river` | — | 1 − `dist_river_norm` | no |
| `W_corridor` | — | `W_river` × (1 − `p_norm`); the hydrological corridor index *W* of the article | yes, as *W* |

The min–max normalisation is (*x* − min *x*) / (max *x* − min *x* + 10⁻¹²) over the 12,826 cells. Each cell takes the
elevation and slope of the single DEM pixel at its centre; the other pixels of the cell are not used.

The archaeological reference dataset is described in the revised manuscript as a compilation of georeferenced
archaeological records from multiple sources: published archaeological investigations and their publicly available
spatial records (Williams and Merino, 2006; Serrudo, 2013; Salcedo Camacho, 2018; Fernandini Parodi, 2023), as well as
records from the Sistema de Información Geográfica de Arqueología (SIGDA), Ministerio de Cultura del Perú. The
records represent heterogeneous archaeological evidence, including sites with architectural remains, surface pottery
concentrations or scatters, and other archaeological materials. This provenance description follows Section 2 and
Table 1 of the revised article; it does not assign all records to a single survey or archaeological project.

The site-labelled cells are grid cells whose centres lie within approximately 150 m of a recorded location (rule
reconstructed from stored distances). A record may label multiple cells, and multiple records may label one cell;
therefore, 84 site-labelled cells represent at least 57 distinct records, not 84 independent archaeological sites.
Unlabelled cells form the background, not confirmed absences. Consistent chronology and documented survey coverage
or intensity are not available for the complete analytical dataset.

The restricted column `distancia_mincul` measures distances to a layer of the Ministerio de Cultura registry; its
presence in the computational grid does not justify attributing every evaluation label to that layer. The original
coordinates, positive-cell labels, and other information enabling reconstruction of archaeological locations are not
included in the publicly distributed datasets.

The site layers themselves (coordinates and attributes), the three site columns and the files derived from them are
not part of the public repository ([README](../README.md#public-restricted-and-third-party-materials)). Information on the location of
archaeological sites may be subject to Peruvian cultural-heritage law; use for purposes other than research may require
authorisation from the Ministerio de Cultura.

### How the grid was built

1. A table with `X`, `Y`, the three site columns, `altitud` and `pendiente` was prepared beforehand from the DEM and the
   site layers (`grilla_final_modelo_logistico.csv`, not included).
2. The script `original_scripts/CA1.py` of release
   [v1.0.0](https://github.com/javier-vz/canete-corridor-model/blob/v1.0.0/original_scripts/CA1.py) added the other six
   columns: it clipped the HydroRIVERS v1.0 lines of South America to the grid extent plus 5 km, computed the distance
   from each cell centre to their union and applied the normalisations above.

[`code/check_grid_columns.py`](../code/check_grid_columns.py) (verification step 3) confirms these rules on the files of
this folder: `altitud` equals the DEM pixel at every cell centre; `pendiente` equals the central-difference slope of that
pixel (largest difference 7.7 × 10⁻⁶ degrees); `dist_river_m` equals the distance to the 81 lines (largest difference
1.8 × 10⁻⁹ m); the normalised and *W* columns follow the formulas above; and in the 84 site-labelled cells the distance to the
nearest survey record is at most 145.8 m (less than half the cell diagonal), while in every other cell it is at least
150.1 m (more than half the cell width).

## DEM (`DEM_vallebajo.tif`)

SRTM 1 arc-second global elevation (NASA Shuttle Radar Topography Mission, distributed by the USGS EROS Center,
doi:[10.5066/F7PR7TFT](https://doi.org/10.5066/F7PR7TFT)), clipped to the study area and reprojected to EPSG:32718 with
30.96 m pixels: 1,944 × 1,579 pixels, 16-bit integers, no-data value −32768 outside the clipped area; extent
76.518–75.961° W and 13.172–12.727° S; elevations from −9 to 3,771 m. The release v1.0.0 documentation states that the
tile was downloaded from USGS EarthExplorer (SRTM version 3).

## River network

`HydroRIVERS_canete_UTM18S.gpkg` contains the geometry of 121 lines of HydroRIVERS v1.0 (Lehner and Grill 2013) around
the study area, without their attribute table. The 81 lines that intersect the grid extent plus 5 km
(`code/river_geometry.load_rivers`) give `dist_river_m` and the river cells: a cell belongs to *R* when its 250 m
footprint intersects one of them. `celdas_fluviales_1020.csv` lists those cells; `code/test_final_h4m_state.py` checks
that it reproduces the stored river-cell mask.

The official HydroRIVERS v1.0 shapefile of South America has 87 reaches that intersect the grid extent plus 5 km: the
81 of this clip and 6 more that lie at least 4 km from every grid cell. Computed from the 87 official reaches, the
distance to the nearest reach and the 1,020 river cells are the same (largest difference in `dist_river_m`
1.8 × 10⁻⁹ m); this was checked against the official file, which is not redistributed here.
