# Cañete corridor model — H4m reproducibility repository

Supporting code, environmental inputs, and saved computational outputs for
*An Unsupervised Corridor-Based Model for Archaeological Prospection in the Lower Cañete Valley, Peru*
(D. S. Cornejo Meza and J. Vera Zúñiga; submitted to the *Journal of Computer Applications in Archaeology*).

## Scientific result and scope

The **primary model output** is the simulated traffic field **F**, a reproducible representation of potential connectivity between **25 environmentally selected attractors**. The final H4m configuration uses **90 attractor pairs**, **four routes per pair**, an eight-neighbour routing grid, and the stochastic parameters documented in the manuscript.

The **composite field L** combines traffic, topographic suitability, and two smoothed traffic components. It is retained as a **secondary comparative field**, not as a calibrated probability of archaeological presence or evidence of predictive superiority.

The analysis domain comprises **12,826 cells of 250 m** (801.625 km²). The evaluation uses **84 site-labelled grid cells**, not 84 independent archaeological sites. The archaeological reference dataset is described in the revised manuscript as a compilation of georeferenced archaeological evidence from **multiple sources**, including published archaeological investigations and available spatial records (Williams and Merino, 2006; Serrudo, 2013; Salcedo, 2018; Fernandini, 2023), and SIGDA records of the Ministerio de Cultura del Perú. The evidence is heterogeneous, encompassing architectural remains, surface pottery concentrations or scatters, and other archaeological materials.

An evaluation cell was labelled when its centre lay within approximately 150 m of a recorded archaeological location. The stored distance information supports a lower bound of **57 distinct records**. The remaining **12,742 cells** constitute background, not confirmed archaeological absences. Archaeological reference records were used **only for evaluation**, never to construct the unsupervised model. The analytical dataset does not provide consistent chronology or documented survey coverage and intensity.

The simulated corridors express **potential connectivity** under the model assumptions; they are not identified historical routes.

## Public, restricted, and third-party materials

- **Public repository content:** H4m model code; non-sensitive environmental inputs; stored simulated traffic and environmental fields; aggregate evaluation tables; documentation; and integrity checks.
- **Not distributed in the current public data files:** exact archaeological point coordinates, the list of the 84 site-labelled cell centres, archaeological presence/absence labels, distances to archaeological records, and files enabling reconstruction of those locations. The environmental grid is published without its three restricted archaeological columns.
- **Third-party environmental data:** SRTM, HydroRIVERS, and HydroSHEDS-derived materials retain their original attribution and redistribution terms. The MIT licence for the authors' code does **not** override third-party terms. See [`THIRD_PARTY.md`](THIRD_PARTY.md).
- **Article and Supplementary Material:** these are separate scholarly outputs. The public repository does not provide the underlying archaeological coordinate dataset. The article's figures and supplementary maps are governed by the manuscript's editorial and cultural-heritage disclosure decisions; their appearance does not imply that underlying georeferenced records are available for download.

**Historical versions need separate attention.** The current public `main` branch contains the revised H4m materials, but the former `v1.0.0` and `v1.0.1` tags and earlier commits may still expose obsolete descriptions or location-derived data. Updating `main` does not erase those historical objects. The repository history has **not** been certified as sanitized. See [`docs/GITHUB_HISTORY_REMEDIATION.md`](docs/GITHUB_HISTORY_REMEDIATION.md).

## Repository structure

```text
code/        Frozen H4m computational code, evaluation scripts, and tests
data/        Public environmental grid, DEM, and river layers (original terms apply)
out/         Stored numerical runs, simulated fields, and aggregated results
revision/    Post-hoc evaluation, regionalisation, and figure outputs
docs/        Procedures, provenance, correspondence, and release audit
tools/       Tools for authorized local restoration of restricted inputs
paper/       Manuscript/supplement status and supporting documentation
```

The historical exploratory archive `pre_paper_v3` is **not** a source of public archaeological records and is not part of this release. Older versions with **24 attractors** and references to **84 sites** do not describe the revised H4m model.

## Quick verification

For checks that do **not** require the archaeological reference data:

```bash
python verify_public.py
```

This verifies the checksums of the public files and checks that the public environmental grid contains **12,826 rows and ten environmental columns**. Passing this check does **not** establish that analyses involving the restricted archaeological labels can be fully recomputed from the public download.

Re-evaluation against the archaeological records requires authorized access to those records. If an authorized restricted bundle is available to a collaborator, use the local restoration workflow documented in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) and [`tools/`](tools/). **Never commit restored archaeological inputs or outputs that reveal locations.** No new simulations are necessary to inspect the stored H4m results.

## Reproducibility and limitations

The repository contains saved H4m computational results and code associated with the revised manuscript. The checks and output correspondence are described in [`docs/REPRODUCE.md`](docs/REPRODUCE.md), [`docs/provenance.md`](docs/provenance.md), and [`docs/outputs_current.md`](docs/outputs_current.md).

The publicly available files support examination of the environmental model and its saved aggregate outputs. **Complete independent re-evaluation of statistics involving archaeological labels is not possible using the public files alone.** This restriction protects archaeological location information and should not be interpreted as the absence of computational documentation.

## Citation and contact

Please cite the accompanying article and software metadata in [`CITATION.cff`](CITATION.cff). No Zenodo DOI should be assumed until one has been assigned. For questions about restricted archaeological records, contact the corresponding author; any permitted access or reuse remains subject to the relevant data custodians and source conditions.

The authors retain responsibility for the interpretation of the model, its archaeological context, and the accuracy of reported results.
