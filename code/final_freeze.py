"""JCAA-80 final recompute: freeze manifest (validation point C). Written after the pre-tests pass and
before S1-S9; `--verificar` re-checks every hash (run before and after evaluation).
Output: out/final/MANIFEST_congelado.sha256"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT)
import argparse, hashlib, time
from pathlib import Path

FILES = [
    # code of the final recompute
    "code/final_common.py", "code/final_runs.py", "code/final_evaluate.py", "code/final_diagnostics.py", "code/final_figures.py",
    "code/final_exports.py", "code/final_export_ppm.py", "code/final_ppm.R", "code/final_h4m_state.py", "code/final_ks_provenance.py",
    "code/final_freeze.py", "code/test_final_pipeline.py", "code/test_final_h4m_state.py", "code/b6_design_counts.py",
    # frozen dependencies
    "code/ca_kernel_v2.py", "code/canonical_partition.py", "code/routing.py", "code/river_geometry.py",
    # inputs
    "data/grilla_con_W_hidrica.csv", "data/celdas_fluviales_1020.csv", "data/HydroRIVERS_canete_UTM18S.gpkg",
    "out/e4/scenario_fields.pkl", "out/b6/design_counts_H4m.csv", "out/b6/fingerprint_H4m.json",
    "out/e3/partitions.npz", "out/e3/regionalization_comparison.csv", "out/kmeans_labels_100.npy", "out/rowmajor_iy_ix.npy",
    "out/rowmajor_sites.npy", "out/state_v2_repro.pkl", "out/river_cell_mask.npy",
    "out/hydro/hydrorivers_model_lines_classified.gpkg", "out/b4/t1/hydrorivers_official_match.csv",
    "out/b5/t5/decision_metrics.csv", "out/b5/t5/attractors_changed.csv",
    *[f"out/b4/runs/{s}_s{k}.npz" for s in ("H4m", "H4m-main", "H1+H2m") for k in (42, 43, 44)],
    *[f"out/b5/runs/H4m-excl8_s{k}.npz" for k in (42, 43, 44)],
    # procedures, pre-tests and the S6 package
    "out/final/PROCEDIMIENTOS.md", "out/final/previa/pruebas_previas.csv", "out/final/previa/pruebas_controles_v2.csv",
    "out/final/previa/ks_procedencia.csv", "out/final/previa/ppm_v2_reproduccion_R_aqui.csv",
    "out/final/previa/prueba_OK_Ryzen_2026-11-06.json", "out/b6/ryzen/JCAA80_barrido_H4m_Ryzen.zip",
]
MAN = Path("out/final/MANIFEST_congelado.sha256")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--verificar", action="store_true"); a = ap.parse_args()
    if a.verificar:
        bad = []
        for line in MAN.read_text().splitlines():
            if line.startswith("#"):
                continue
            h, f = line.split("  ", 1)
            if not Path(f).exists() or sha(f) != h:
                bad.append(f)
        print("manifest OK" if not bad else f"CHANGED: {bad}")
        _sys.exit(1 if bad else 0)
    missing = [f for f in FILES if not Path(f).exists()]
    assert not missing, missing
    MAN.write_text(f"# frozen {time.strftime('%Y-%m-%d %H:%M:%S %Z')} - JCAA-80 final recompute (H4m)\n"
                   + "\n".join(f"{sha(f)}  {f}" for f in FILES) + "\n")
    print(f"{len(FILES)} files frozen -> {MAN} ({sha(MAN)[:12]})")
