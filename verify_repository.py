"""Verify the stored results of the article without running new simulations.

Usage, from the root of the repository:

    python verify_repository.py          # about 3 minutes
    python verify_repository.py --keep   # keep the temporary copy

The checks run on a temporary copy of code/, data/ and out/, so the repository itself is never written.

 1. sha256 of every file listed in MANIFEST.sha256.
 2. Freeze manifest (code/final_freeze.py --verificar): only code/test_final_pipeline.py may differ (documented change).
 3. Columns of the analysis grid against the DEM and the HydroRIVERS lines (code/check_grid_columns.py).
 4. Final state rebuilt and compared with the stored states (code/test_final_h4m_state.py).
 5. Stored runs valid by content (code/final_runs.py --verificar).
 6. Evaluation tables recomputed from the stored runs (code/final_evaluate.py).
 7. Diagnostics recomputed (code/b9_diag_E.py, code/b9_diag_atractores.py, code/diag_L_vs_T.py,
    code/b6_design_counts.py).
 8. No other file of the copy was created or changed.

A recomputed file passes if it is byte-identical to the stored one or, as can happen on another processor or operating
system, if its text cells are identical and its numeric cells differ by at most 1e-9 in relative terms. Step 4 repeats a
failed exact comparison of floating-point fields with a tolerance of 1e-12. See docs/reproduce.md.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCOPE = ("code", "data", "out")
RTOL, ATOL, STATE_TOL = 1e-9, 1e-15, 1e-12
DOCUMENTED_FREEZE_CHANGES = ["code/test_final_pipeline.py"]
FLOAT_STATE_CHECKS = ("E' identical", "c' identical", "T identical", "graph edge weights identical", "same as def_H4m.pkl")
DIAGNOSTICS = {
    "code/b9_diag_E.py": ["out/final/informe/diagnostico_E_atractores_sitios.csv", "out/final/informe/rangos_campos_figura4.csv"],
    "code/b9_diag_atractores.py": ["out/final/informe/atractores_sensibilidades.csv", "out/final/informe/conteo_atractores_8_lineas.csv"],
    "code/diag_L_vs_T.py": ["out/final/informe/diagnostico_L_vs_T.csv"],
    "code/b6_design_counts.py": ["out/b6/design_counts_H4m.csv"],
}
FINGERPRINT = "out/b6/fingerprint_H4m.json"

STATE_TOLERANCE_CHECK = r"""
import pickle, sys
import numpy as np
sys.path.insert(0, "code")
import final_h4m_state as fs
from h_scenarios import make_state
CSV, RIV = "data/grilla_con_W_hidrica.csv", "data/celdas_fluviales_1020.csv"
def edge_diff(G1, G2):
    assert G1.number_of_edges() == G2.number_of_edges()
    return max(max(abs(G1.edges[e]["weight_base"] - G2.edges[e]["weight_base"]),
                   abs(G1.edges[e]["weight"] - G2.edges[e]["weight"])) for e in G1.edges)
same, worst = True, 0.0
for scen, kw in [("H4m", {}), ("H4m-main", {"cost_mask": "main"})]:
    st = fs.build_state_h4m(CSV, RIV, **kw); ref = make_state(scen); V = st.valid_mask
    same &= st.attractors == ref.attractors and st.pair_candidates == ref.pair_candidates
    for a, b in ((st.E_grid, ref.E_grid), (st.c_grid, ref.c_grid), (st.T_grid, ref.T_grid)):
        worst = max(worst, float(np.abs(a - b)[V].max()))
    worst = max(worst, edge_diff(st.G_template, ref.G_template))
    if scen == "H4m":
        b4 = pickle.load(open("out/b4/runs/def_H4m.pkl", "rb"))
        same &= st.attractors == b4["attractors"] and st.pair_candidates == b4["pairs"]
        worst = max(worst, float(np.abs(st.E_grid - b4["E"])[V].max()), float(np.abs(st.c_grid - b4["c"])[V].max()))
print("RESULT", same, repr(worst))
"""

results = []


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def report(name, ok, detail=""):
    results.append(bool(ok))
    print(f"[{'OK  ' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}", flush=True)


def note(text):
    print(f"       {text}", flush=True)


def read_manifest():
    entries = {}
    for line in (ROOT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            digest, rel = line.split("  ", 1)
            entries[rel.strip()] = digest
    return entries


def run(script_args, cwd):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", JCAA_ROOT=str(cwd), PYTHONIOENCODING="utf-8")
    proc = subprocess.run([sys.executable] + script_args, cwd=cwd, env=env, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")
    return proc.returncode, (proc.stdout + proc.stderr).strip().splitlines()


def numbers_close(a, b):
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False, 0.0
    if math.isnan(x) and math.isnan(y):
        return True, 0.0
    if math.isinf(x) or math.isinf(y):
        return x == y, 0.0
    diff = abs(x - y)
    scale = max(abs(x), abs(y))
    rel = diff / scale if scale > 0 else 0.0
    return diff <= ATOL or diff <= RTOL * scale, rel


def compare_csv(new, ref):
    """(equal within rounding, largest relative difference, reason)."""
    with open(new, newline="", encoding="utf-8") as f1, open(ref, newline="", encoding="utf-8") as f2:
        rows_new, rows_ref = list(csv.reader(f1)), list(csv.reader(f2))
    if len(rows_new) != len(rows_ref):
        return False, 0.0, f"{len(rows_new)} rows instead of {len(rows_ref)}"
    worst = 0.0
    for i, (r1, r2) in enumerate(zip(rows_new, rows_ref)):
        if len(r1) != len(r2):
            return False, worst, f"row {i}: {len(r1)} cells instead of {len(r2)}"
        for c1, c2 in zip(r1, r2):
            if c1 == c2:
                continue
            ok, rel = numbers_close(c1, c2)
            if not ok:
                return False, worst, f"row {i}: {c1!r} instead of {c2!r}"
            worst = max(worst, rel)
    return True, worst, ""


def compare_recomputed(copy, rel, digest):
    """Compare a recomputed file of the copy with the stored file; returns (status, detail)."""
    new, ref = copy / rel, ROOT / rel
    if not new.is_file():
        return "missing", "not written"
    if sha256(new) == digest:
        return "identical", ""
    if new.suffix == ".csv":
        ok, worst, why = compare_csv(new, ref)
        return ("rounding", f"largest relative difference {worst:.1e}") if ok else ("different", why)
    return "different", "bytes differ"


def compare_index(copy, statuses):
    """indice.json lists the sha256 of every table it wrote: identical tables must keep their hash, and the hash of a
    table that differs only by rounding must be that of the recomputed file."""
    rel = "out/final/eval/indice.json"
    new, old = json.loads((copy / rel).read_text(encoding="utf-8")), json.loads((ROOT / rel).read_text(encoding="utf-8"))
    if set(new) != set(old) or any(new[k] != old[k] for k in new if k != "archivos") or set(new["archivos"]) != set(old["archivos"]):
        return "different", "signature, procedure hash or file list differ"
    for name, digest in new["archivos"].items():
        status = statuses.get(f"out/final/eval/{name}", "identical")
        expected = old["archivos"][name] if status == "identical" else sha256(copy / "out/final/eval" / name)
        if digest != expected:
            return "different", f"hash of {name}"
    return ("identical", "") if sha256(copy / rel) == sha256(ROOT / rel) else ("rounding", "file hashes of rounded tables")


def environment():
    print(f"Python {platform.python_version()} on {platform.system()} {platform.machine()}")
    pins = {}
    for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines():
        if "==" in line and not line.lstrip().startswith("#"):
            name, version = line.split("#")[0].strip().split("==")
            pins[name.strip()] = version.strip()
    off = []
    for name, version in pins.items():
        try:
            have = metadata.version(name)
        except metadata.PackageNotFoundError:
            have = "missing"
        if have != version:
            off.append(f"{name} {have} (pinned {version})")
    if off:
        print("Note: versions differ from requirements.txt: " + "; ".join(off))
    print(flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--keep", action="store_true", help="keep the temporary copy")
    args = ap.parse_args()
    t0 = time.time()
    environment()
    manifest = read_manifest()
    in_scope = {p.relative_to(ROOT).as_posix() for d in SCOPE for p in (ROOT / d).rglob("*")
                if p.is_file() and p.name != "README.md" and "__pycache__" not in p.parts}

    # 1. checksums
    bad = [p for p, h in manifest.items() if not (ROOT / p).is_file() or sha256(ROOT / p) != h]
    report(f"1. checksums of the {len(manifest)} files of MANIFEST.sha256", not bad, f"differing or missing: {bad[:5]}" if bad else "")
    extra = sorted(in_scope - set(manifest))
    if extra:
        note(f"{len(extra)} files not in the manifest are ignored, e.g. {extra[:3]}")

    tmp = Path(tempfile.mkdtemp(prefix="canete_verify_"))
    copy = tmp / "repo"
    for rel in manifest:
        (copy / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, copy / rel)
    try:
        # 2. freeze manifest
        rc, out = run(["code/final_freeze.py", "--verificar"], copy)
        changed = None
        for line in out:
            if line.startswith("CHANGED:"):
                changed = sorted(json.loads(line.split(":", 1)[1].strip().replace("'", '"')))
        ok = changed == DOCUMENTED_FREEZE_CHANGES or (rc == 0 and changed is None)
        report("2. code, data and procedures match the freeze manifest", ok,
               "only the documented change to code/test_final_pipeline.py" if ok else f"({out[-1] if out else 'no output'})")

        # 3. grid columns
        rc, out = run(["code/check_grid_columns.py"], copy)
        report("3. grid columns match the DEM and the HydroRIVERS lines", rc == 0 and out and out[-1] == "GRID CHECKS PASSED",
               "" if rc == 0 else "; ".join(l for l in out if l.startswith("[FAIL]")))

        # 4. final state
        rc, out = run(["code/test_final_h4m_state.py"], copy)
        fails = [l for l in out if l.startswith("[FAIL]")]
        if rc == 0 and out and out[-1].strip() == "ALL OK":
            fp_same = sha256(copy / FINGERPRINT) == manifest[FINGERPRINT]
            report("4. final state rebuilt identical to the stored states", fp_same, "" if fp_same else "fingerprint differs")
        elif fails and all(any(k in l for k in FLOAT_STATE_CHECKS) for l in fails):
            rc2, out2 = run(["-c", STATE_TOLERANCE_CHECK], copy)
            res = [l for l in out2 if l.startswith("RESULT")]
            same, worst = (res[0].split()[1] == "True", float(res[0].split()[2])) if res else (False, float("nan"))
            report("4. final state rebuilt: attractors and pairs identical, fields equal within rounding",
                   rc2 == 0 and same and worst <= STATE_TOL, f"largest difference {worst:.1e}")
        else:
            report("4. final state rebuilt identical to the stored states", False, "; ".join(fails[:3]) or (out[-1] if out else ""))

        # 5. stored runs
        rc, out = run(["code/final_runs.py", "--verificar", "--bloques", "S1", "S2", "S3", "S4", "S5", "S7", "S8", "S9", "CTRL"], copy)
        summary = next((l for l in out if l.startswith("firma")), out[-1] if out else "")
        m = re.search(r"trabajos (\d+): válidos (\d+), inválidos (\d+), faltan (\d+)", summary)
        report("5. stored runs valid by content", rc == 0,
               f"{m.group(2)} of {m.group(1)} valid, {m.group(3)} invalid, {m.group(4)} missing" if m else summary)

        # 6. evaluation tables
        rc, out = run(["code/final_evaluate.py"], copy)
        index = json.loads((ROOT / "out/final/eval/indice.json").read_text(encoding="utf-8"))
        written = [f"out/final/eval/{name}" for name in index["archivos"]]
        statuses, details = {}, {}
        for rel in written:
            statuses[rel], details[rel] = compare_recomputed(copy, rel, manifest[rel])
        statuses["out/final/eval/indice.json"], details["out/final/eval/indice.json"] = compare_index(copy, statuses)
        n_id = sum(s == "identical" for s in statuses.values())
        n_round = [r for r, s in statuses.items() if s == "rounding"]
        bad = [f"{r} ({details[r]})" for r, s in statuses.items() if s not in ("identical", "rounding")]
        report(f"6. evaluation recomputed from the stored runs: {len(statuses)} files", rc == 0 and not bad,
               f"{n_id} byte-identical" + (f", {len(n_round)} equal within rounding" if n_round else "") + (f"; different: {bad[:3]}" if bad else ""))
        for r in n_round:
            note(f"{r}: {details[r]}")
        if rc != 0:
            note(out[-1] if out else "final_evaluate.py failed")

        # 7. diagnostics
        diag_status, bad = {}, []
        for script, outputs in DIAGNOSTICS.items():
            rc, out = run([script], copy)
            if rc != 0:
                bad.append(f"{script} failed ({out[-1] if out else ''})")
            for rel in outputs:
                diag_status[rel], d = compare_recomputed(copy, rel, manifest[rel])
                if diag_status[rel] not in ("identical", "rounding"):
                    bad.append(f"{rel} ({d})")
                elif diag_status[rel] == "rounding":
                    note(f"{rel}: {d}")
        n_id = sum(s == "identical" for s in diag_status.values())
        report(f"7. diagnostics recomputed: {len(diag_status)} files", not bad,
               f"{n_id} byte-identical" + (f", {len(diag_status) - n_id} equal within rounding" if n_id < len(diag_status) else "")
               + (f"; problems: {bad[:3]}" if bad else ""))

        # 8. nothing else changed
        rewritten = set(written) | {"out/final/eval/indice.json", FINGERPRINT} | set(diag_status)
        now = {p.relative_to(copy).as_posix() for p in copy.rglob("*") if p.is_file()}
        new_files = sorted(now - set(manifest))
        changed = sorted(r for r in manifest if r not in rewritten and (not (copy / r).is_file() or sha256(copy / r) != manifest[r]))
        report("8. no other file created or changed", not new_files and not changed,
               (f"new: {new_files[:5]} " if new_files else "") + (f"changed: {changed[:5]}" if changed else ""))
    finally:
        if args.keep:
            print(f"\nTemporary copy kept at {copy}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    print(f"\n{'ALL CHECKS PASSED' if all(results) else 'SOME CHECKS FAILED'} ({(time.time() - t0) / 60:.1f} min)")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
