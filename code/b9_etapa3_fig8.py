"""JCAA-80, Stage 3: Figure 8 with the diagonal labelled as a uniform-density reference (adjustment 3 of the Stage 2
approval). Report layer: the frozen script code/final_figures.py and its outputs in out/final/fig/ are not modified.

The figure is drawn by the frozen function final_figures.fig8() from the same CSVs (out/final/eval/zonas_prefijo.csv
and zonal_P2.csv); the only change is the legend text of the diagonal, substituted while the function runs.
Output: out/final/fig_etapa3/fig8_captura_zonas_ref_uniforme.(png|pdf) (manuscript file name: vz2_multi_hitrate.png).

Checks:
 1. the frozen function, run unchanged into a scratch folder, reproduces the frozen PNG pixel for pixel;
 2. the new PNG differs from the frozen one only inside the legend boxes of the two versions (the curves, points,
    axes and labels are identical);
 3. the step curves read back from the CSVs end at 100 % of the sites for n = 20, 30 and 40.
A JSON summary is written next to the figure."""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import hashlib, json, tempfile
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from PIL import Image
import final_figures as ff

OLD_LABEL, NEW_LABEL = "Random expectation", "Uniform-density reference"
FROZEN = Path("out/final/fig/fig8_captura_zonas.png")
OUTD = Path("out/final/fig_etapa3"); OUTD.mkdir(parents=True, exist_ok=True)
OUT = OUTD / "fig8_captura_zonas_ref_uniforme"
zonal = pd.read_csv(ff.EV / "zonal_P2.csv").set_index("n")
prefix = pd.read_csv(ff.EV / "zonas_prefijo.csv")


def run(path, relabel):
    """Run the frozen fig8(); optionally swap the diagonal's legend text; return the legend box in PNG pixels."""
    boxes, _plot, _close = [], Axes.plot, plt.close

    def plot(self, *a, **k):
        if relabel and k.get("label") == OLD_LABEL:
            k["label"] = NEW_LABEL
        return _plot(self, *a, **k)

    def close(fig=None):
        if fig is not None:                                    # legend box at the savefig resolution (300 dpi)
            fig.set_dpi(300); fig.canvas.draw()
            bb = fig.axes[0].get_legend().get_window_extent(fig.canvas.get_renderer())
            H = fig.bbox.height
            boxes.append((bb.x0, H - bb.y1, bb.x1, H - bb.y0))
        return _close(fig)

    Axes.plot, plt.close = plot, close
    try:
        ff.fig8(prefix, zonal, path)
    finally:
        Axes.plot, plt.close = _plot, _close
    return boxes[0]


with tempfile.TemporaryDirectory() as td:
    box_old = run(Path(td) / "fig8_rerun", relabel=False)
    rerun = np.asarray(Image.open(Path(td) / "fig8_rerun.png").convert("RGBA"))
box_new = run(OUT, relabel=True)

frozen = np.asarray(Image.open(FROZEN).convert("RGBA"))
new = np.asarray(Image.open(OUT.with_suffix(".png")).convert("RGBA"))
assert rerun.shape == frozen.shape and np.array_equal(rerun, frozen), "the frozen function no longer reproduces the frozen PNG"
assert new.shape == frozen.shape, (new.shape, frozen.shape)
diff = np.any(new != frozen, axis=2)
ys, xs = np.nonzero(diff)
pad = 3                                                                # anti-aliasing margin, pixels
x0 = min(box_old[0], box_new[0]) - pad; y0 = min(box_old[1], box_new[1]) - pad
x1 = max(box_old[2], box_new[2]) + pad; y1 = max(box_old[3], box_new[3]) + pad
inside = (xs >= x0) & (xs <= x1) & (ys >= y0) & (ys <= y1)
assert diff.any() and inside.all(), f"{(~inside).sum()} differing pixels outside the legend"
for n in (20, 30, 40):
    z = prefix[prefix.n == n].sort_values("rank")
    assert abs(z.cum_sites_pct.iloc[-1] - 100) < 1e-9 and abs(z.cum_area_pct.iloc[-1] - 100) < 1e-6, n

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
summary = dict(script="code/b9_etapa3_fig8.py", frozen_function="code/final_figures.py:fig8", frozen_function_sha256=sha("code/final_figures.py"),
               sources={str(ff.EV / "zonas_prefijo.csv"): sha(ff.EV / "zonas_prefijo.csv"), str(ff.EV / "zonal_P2.csv"): sha(ff.EV / "zonal_P2.csv")},
               legend_change=[OLD_LABEL, NEW_LABEL], frozen_png=str(FROZEN), frozen_png_sha256=sha(FROZEN),
               frozen_png_reproduced_pixel_identical=True, size_px=list(new.shape[:2][::-1]), dpi=300,
               differing_pixels=int(diff.sum()), differing_pixels_outside_legend=int((~inside).sum()),
               legend_box_px={"frozen": [round(v, 1) for v in box_old], "new": [round(v, 1) for v in box_new]},
               rho_in_legend={int(n): round(float(zonal.loc[n, "rho"]), 2) for n in (20, 30, 40)},
               outputs={str(OUT.with_suffix(e)): sha(OUT.with_suffix(e)) for e in (".png", ".pdf")})
(OUTD / "fig8_comprobacion.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

# side-by-side view of the legend area for the visual check (frozen | new)
cx0, cy0, cx1, cy1 = int(x0) - 40, int(y0) - 40, int(x1) + 20, int(y1) + 20
crop = lambda a: a[max(cy0, 0):cy1, max(cx0, 0):cx1]
side = np.concatenate([crop(frozen), np.full((crop(frozen).shape[0], 12, 4), 255, np.uint8), crop(new)], axis=1)
Image.fromarray(side).save(OUTD / "fig8_leyenda_antes_despues.png")
print(f"Figure 8 -> {OUT.with_suffix('.png')} | frozen PNG reproduced; {diff.sum()} differing pixels, all inside the legend boxes "
      f"| rho {summary['rho_in_legend']}")
