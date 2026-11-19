# Public release gate — archaeology and licensing

State: **OPEN**. Not approved for publishing.

- [x] Public environmental grid has 12,826 cells and ten columns; it does not contain the three
  restricted archaeological columns.
- [x] The existing candidate excludes eight known restricted files, including the 84-cell GIS layer.
- [x] `python verify_public.py` passes before any additional publication.
- [ ] Project director / authorized data custodian confirms allowable public representation of
  archaeological locations, including derived raster cells and manuscript/supplement figures.
- [ ] Main PDF Figures 5/6 and Supplement Figure S3 pass the agreed disclosure policy.
- [ ] Historical `main`, tags `v1.0.0` and `v1.0.1`, releases and commit tree are remediated and
  independent access to previously exposed records has been assessed.
- [ ] Third-party HydroRIVERS/HydroSHEDS redistribution and required attribution are approved
  and documented. Author MIT licence does not cover them.
- [ ] Author confirms attribution and data access policy, and JCAA editorial requirements.
- [ ] A fresh clean-clone check and a re-run of `verify_public.py` pass.

Avoid sharing working archives with raw site distances, `.RData`, georeferenced positive-cell
layers or unreviewed maps. Do not upload the private Supplementary Material until Figures S3
and all author notes are resolved.
