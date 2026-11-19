"""JCAA-80 final recompute: PPM input (replaces export_fields_v2_for_ppm.py). Same format and row
order as fields_canonical_v2_for_ppm.csv (grid CSV order; X, Y, F = F_norm, L, E, is_site) with
F and L of the S1 canonical run (seed 42) and E = E' (final affinity field). Then run
code/final_ppm.R inside out/final/ppm (Rscript final_ppm.R).
Regression-tested: export(v2 fields) reproduces the published v2 PPM input to 1e-15."""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import shutil
import numpy as np
import pandas as pd


def export(F, L, E, st, path):
    df = pd.read_csv("data/grilla_con_W_hidrica.csv")
    lat = st.lattice
    iy = np.array([lat["y_to_iy"][y] for y in df.Y]); ix = np.array([lat["x_to_ix"][x] for x in df.X])
    out = pd.DataFrame({"X": df.X, "Y": df.Y, "F": F[iy, ix], "L": L[iy, ix], "E": E[iy, ix],
                        "is_site": df.presencia_sitio_nuevo.astype(int)})
    out.to_csv(path, index=False)
    return out


if __name__ == "__main__":
    import final_common as fc, final_runs as fr
    from final_evaluate import fget, check_block
    check_block("S1")
    st = fc.state_for(fr.cfg())
    F, L = fget("S1", "canon_s42")
    d = fc.OUT / "ppm"; d.mkdir(parents=True, exist_ok=True)
    o = export(F, L, st.E_grid, st, d / "fields_H4m_for_ppm.csv")
    shutil.copy2("code/final_ppm.R", d / "final_ppm.R")
    print(len(o), "rows,", int(o.is_site.sum()), "sites ->", d / "fields_H4m_for_ppm.csv", "| next: cd out/final/ppm && Rscript final_ppm.R")
