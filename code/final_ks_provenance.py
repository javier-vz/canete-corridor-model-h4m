"""KS provenance (validation point A): reproduce the historical implementation behind the published
Table 6 to document where its numbers came from. It does NOT choose the final method (fixed a priori
in out/final/PROCEDIMIENTOS.md, section 6). Fields: E (same in v1 and v2), L of v1 and v2 (seed 42).
Residual-set definitions: (a) historical analyze_ks_test.py = sites outside the largest connected
component of the sites dilated by 3 cells; (b) final = sites outside the P2 zone with most sites.
Methods: SciPy 'auto', 'asymp', 'exact'. Output: out/final/previa/ks_procedencia.csv"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import pickle, warnings
import numpy as np, pandas as pd
from scipy.stats import ks_2samp
from scipy.ndimage import binary_dilation, label as scilabel
from canonical_partition import canonical_partition

OUT = _Path("out/final/previa"); OUT.mkdir(parents=True, exist_ok=True)
st = pickle.load(open("out/state_v2_repro.pkl", "rb")); V = st.valid_mask; S = st.site_mask
bg = V & ~S
# (a) historical dominant cluster
dil = binary_dilation(S, structure=np.ones((7, 7)))
lab, n = scilabel(dil); sizes = np.bincount(lab.ravel())[1:]
dom_a = (lab == int(sizes.argmax()) + 1) & S
# (b) P2 zone with most sites
P = canonical_partition(st)
zt = P.groupby("zone").presencia_sitio_nuevo.sum(); zdom = zt.idxmax()
dom_b = np.zeros_like(S); m = (P.zone.values == zdom) & (P.presencia_sitio_nuevo.values == 1)
dom_b[P.iy.values[m], P.ix.values[m]] = True
print(f"dominant (a) dilation: {dom_a.sum()} sites; (b) P2 zone {zdom}: {dom_b.sum()} sites; identical sets: {np.array_equal(dom_a, dom_b)}")

PUB = {("all84", "E"): (0.5309, 2.317e-22), ("all84", "L"): (0.5430, 1.867e-23),
       ("residual", "E"): (0.4336, 6.171e-4), ("residual", "L"): (0.2934, 0.05084)}
fields = {"E": st.E_grid, "L_v1": np.load("out/L_v1_s42.npy"), "L_v2": np.load("out/L_v2_s42.npy")}
rows = []
for fname, f in fields.items():
    for sset, mask in [("all84", S), ("residual_a_dilation", S & ~dom_a), ("residual_b_P2", S & ~dom_b)]:
        for meth in ("auto", "asymp", "exact"):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                r = ks_2samp(f[mask], f[bg], alternative="two-sided", method=meth)
            pub = PUB.get(("all84" if sset == "all84" else "residual", fname[0]))
            rows.append(dict(field=fname, site_set=sset, n_sites=int(mask.sum()), n_background=int(bg.sum()), method_requested=meth,
                             fallback_warning=any("exact" in str(x.message).lower() for x in w), D=r.statistic, p=r.pvalue,
                             D_published=pub[0], p_published=pub[1],
                             reproduces_published=bool(abs(r.statistic - pub[0]) < 5e-5 and abs(r.pvalue - pub[1]) / pub[1] < 5e-4)))
D = pd.DataFrame(rows); D.to_csv(OUT / "ks_procedencia.csv", index=False)
pd.set_option("display.width", 200)
print(D[["field", "site_set", "n_sites", "method_requested", "fallback_warning", "D", "p", "reproduces_published"]].to_string(index=False))
