"""Restore the restricted files into this repository from an authorised bundle (see tools/make_restricted_bundle.py).

    python tools/restore_restricted.py --bundle PATH_TO_BUNDLE

Rebuilds data/grilla_con_W_hidrica.csv (and its identical copy in out/b6/) byte for byte from the public grid and the
site columns, copies the six other files, and checks all eight against RESTRICTED_MANIFEST.sha256. Afterwards,
python verify_repository.py runs the complete verification. Restored files must not be committed to a public repository
(they are listed in .gitignore)."""
import argparse, shutil, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import restricted as R

ap = argparse.ArgumentParser()
ap.add_argument("--bundle", required=True, type=Path)
a = ap.parse_args()
man = R.restricted_manifest()
full = R.merge_grid((R.ROOT / R.GRID_PUBLIC).read_bytes(), (a.bundle / R.SITE_FILE).read_bytes())
assert R.sha256_bytes(full) == man[R.GRID], "the rebuilt grid does not match RESTRICTED_MANIFEST.sha256"
for rel in R.GRID_COPIES:
    (R.ROOT / rel).parent.mkdir(parents=True, exist_ok=True)
    (R.ROOT / rel).write_bytes(full)
for rel in R.WHOLE_FILES:
    (R.ROOT / rel).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(a.bundle / rel, R.ROOT / rel)
bad = [rel for rel, h in man.items() if R.sha256(R.ROOT / rel) != h]
assert not bad, f"checksums differ: {bad}"
print(f"RESTRICTED FILES RESTORED: {len(man)} files match RESTRICTED_MANIFEST.sha256. Run: python verify_repository.py")
