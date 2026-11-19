"""Check the public files of the repository (no restricted data needed, standard library only, a few seconds).

    python verify_public.py

1. sha256 of every file listed in MANIFEST_PUBLIC.sha256 (all files of the public repository);
2. the public analysis grid has the 12,826 cells and the ten public columns, and no site column;
3. reports which restricted files (RESTRICTED_MANIFEST.sha256) are present. When all of them are present and valid,
   python verify_repository.py runs the complete verification (rebuilds the final model state and recomputes every
   evaluation table from the stored runs, without new simulations)."""
import hashlib, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tools"))
import restricted as R


def manifest(name):
    rows = {}
    for line in (ROOT / name).read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            h, p = line.split(maxsplit=1)
            rows[p.strip()] = h
    return rows


ok = True
pub = manifest("MANIFEST_PUBLIC.sha256")
bad = [p for p, h in pub.items() if not (ROOT / p).is_file() or R.sha256(ROOT / p) != h]
print(f"[{'PASS' if not bad else 'FAIL'}] 1. checksums of the {len(pub)} public files" + (f": differing or missing {bad[:5]}" if bad else ""))
ok &= not bad
lines = (ROOT / R.GRID_PUBLIC).read_bytes().splitlines()
head = lines[0].split(b",")
g_ok = len(lines) == 12827 and not set(head) & set(R.SITE_COLUMNS) and head == [c for c in R.ORIGINAL_HEADER if c not in R.SITE_COLUMNS]
print(f"[{'PASS' if g_ok else 'FAIL'}] 2. public grid: {len(lines) - 1} cells, {len(head)} columns, no site column")
ok &= g_ok
res = R.restricted_manifest(ROOT)
present = {p: (ROOT / p).is_file() and R.sha256(ROOT / p) == h for p, h in res.items() if (ROOT / p).exists()}
if not present:
    print(f"[INFO] 3. restricted files: none present ({len(res)} files, see README). The complete verification needs them.")
elif len(present) == len(res) and all(present.values()):
    print(f"[INFO] 3. restricted files: all {len(res)} present and valid. Run: python verify_repository.py")
else:
    print(f"[FAIL] 3. restricted files: {sum(present.values())} of {len(res)} valid; restore them with tools/restore_restricted.py")
    ok = False
print("PUBLIC CHECKS PASSED" if ok else "PUBLIC CHECKS FAILED")
sys.exit(0 if ok else 1)
