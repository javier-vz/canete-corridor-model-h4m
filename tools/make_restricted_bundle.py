"""Build the bundle of restricted data from a complete copy of the repository (release v2.0.0 or a restored copy).

    python tools/make_restricted_bundle.py --full PATH_TO_COMPLETE_REPOSITORY --out PATH_TO_BUNDLE

The bundle holds grid_site_columns.csv (X, Y and the three site columns of the analysis grid) and the six other
restricted files at their paths in the repository. It must not be published: it gives the location of archaeological
sites. Every file is checked against RESTRICTED_MANIFEST.sha256 of this repository."""
import argparse, shutil, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import restricted as R

ap = argparse.ArgumentParser()
ap.add_argument("--full", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path)
a = ap.parse_args()
man = R.restricted_manifest()
grid = (a.full / R.GRID).read_bytes()
assert R.sha256_bytes(grid) == man[R.GRID], "the analysis grid of the complete copy is not the expected one"
pub, site = R.split_grid(grid)
if (R.ROOT / R.GRID_PUBLIC).exists():
    assert pub == (R.ROOT / R.GRID_PUBLIC).read_bytes(), "the public grid of this repository does not match the complete copy"
a.out.mkdir(parents=True, exist_ok=True)
(a.out / R.SITE_FILE).write_bytes(site)
for rel in R.WHOLE_FILES:
    src = a.full / rel
    assert R.sha256(src) == man[rel], f"{rel}: checksum differs"
    (a.out / rel).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, a.out / rel)
print(f"BUNDLE WRITTEN: {a.out} ({1 + len(R.WHOLE_FILES)} files). Do not publish it.")
