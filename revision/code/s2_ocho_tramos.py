"""S2: audit of the eight unverified HydroRIVERS reaches (no model run). Read-only on the frozen outputs.

Inputs: the official HydroRIVERS v1.0 South America shapefile (with its .dbf, in the authors' folder), the record numbers
(FID) of the eight reaches used by final_evaluate.py (EXCL_FID), the line match out/b4/t1/hydrorivers_official_match.csv
(dataset C87), the classified model lines, Dina's review form (8 priority lines, received 5 Nov 2026), the canonical
traffic field (S1, seed 42), the H4m state and out/final/informe/conteo_atractores_8_lineas.csv.
Outputs (03_S2_HIDROGRAFIA/): S2_OCHO_TRAMOS.csv, S2_ocho_tramos.gpkg (reaches with the review, a separate uncertainty layer and
the attractors of both scenarios; no archaeological sites), S2_MAPA.pdf/png, S2_log.txt.
Status rule: 'no verificable' when the opinion is missing, says it cannot be determined or repeats the observation; 'verificado'
only with high confidence; otherwise 'interpretación tentativa' with the stated confidence."""
import os, sys
from pathlib import Path
W = Path(__file__).resolve().parents[2]; os.chdir(W); sys.path.insert(0, str(W / "code"))
sys.dont_write_bytecode = True
import numpy as np, pandas as pd, geopandas as gpd, pyogrio, shapely
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import final_common as fc, final_runs as fr

OUT = W / "revision/out/S4_reaches"; OUT.mkdir(parents=True, exist_ok=True)
SHP = Path(os.environ["HYDRORIVERS_SA_SHP"])   # HydroRIVERS v1.0, South America (HydroSHEDS; not redistributed)
FORM = Path(os.environ["REACH_REVIEW_FORM"])   # imagery review form (not distributed)
EXCL_FID = [1131506, 1130585, 1129684, 1130586, 1127717, 1130415, 1126308, 1126640]   # = final_evaluate.py / b5_exclusion_runs.py
log = []

# 1. official records of the eight reaches
off = pyogrio.read_dataframe(SHP, fids=EXCL_FID).to_crs(32718)
off["shp_fid"] = EXCL_FID
cl = gpd.read_file("out/hydro/hydrorivers_model_lines_classified.gpkg")
match = pd.read_csv("out/b4/t1/hydrorivers_official_match.csv").query("dataset == 'C87'").set_index("line")
cl["shp_fid"] = [int(match.shp_fid[i]) for i in range(len(cl))]
cl["line"] = range(len(cl))
sel = cl[cl.shp_fid.isin(EXCL_FID)].merge(off.drop(columns="geometry"), on="shp_fid")
for _, r in sel.iterrows():
    hd = r.geometry.hausdorff_distance(off.set_index("shp_fid").geometry[r.shp_fid])
    log.append(f"line {r.line}: FID {r.shp_fid} -> HYRIV_ID {int(r.HYRIV_ID)}; Hausdorff to the model line {hd:.2f} m")
    assert hd < 1.0, "the official record does not match the model line"

# 2. Dina's review (sheet '8 prioritarias')
raw = pd.read_excel(FORM, sheet_name="8 prioritarias", header=None)
h = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "N°"][0]
dina = raw.iloc[h + 1:].copy(); dina.columns = [str(c).strip() for c in raw.iloc[h]]
dina = dina[pd.to_numeric(dina["N°"], errors="coerce").notna()].copy()
d56 = pd.read_csv("out/b4/t1/hydrorivers_official_match.csv").query("dataset == 'D56'").set_index("line").shp_fid   # Dina's package numbering
dina["linea_Dina"] = dina["Línea"].astype(int)
dina["shp_fid"] = [int(d56[l]) for l in dina.linea_Dina]
assert set(dina.shp_fid) == set(EXCL_FID), "the eight lines of the form are not the eight excluded reaches"
dina["line"] = [int(sel.set_index("shp_fid").line[f]) for f in dina.shp_fid]
log.append("the 8 lines of the form (numbering of Dina's package, D56) are the 8 excluded reaches: " + ", ".join(f"{a}->FID {b}" for a, b in zip(dina.linea_Dina, dina.shp_fid)))

# 3. traffic attributed to each reach (cells whose footprint intersects it), canonical run S1 seed 42
st = fc.state_for(fr.cfg())
F = np.load("out/final/runs/S1/canon_s42.npz")["F"]
lat = st.lattice
df = st.df
cells = shapely.box(df.X.values - 125, df.Y.values - 125, df.X.values + 125, df.Y.values + 125)
iy = np.array([lat["y_to_iy"][y] for y in df.Y]); ix = np.array([lat["x_to_ix"][x] for x in df.X])
Fc = F[iy, ix]
mass, ncell = {}, {}
for _, r in sel.iterrows():
    hit = shapely.intersects(cells, r.geometry)
    mass[r.line] = 100 * Fc[hit].sum() / Fc.sum(); ncell[r.line] = int(hit.sum())
allhit = shapely.intersects(cells, shapely.union_all(sel.geometry.values))
log.append(f"traffic of seed 42 on the cells of the eight reaches: {100 * Fc[allhit].sum() / Fc.sum():.2f}% ({int(allhit.sum())} cells)")

# 4. table S2
t11 = pd.read_csv("out/final/eval/tabla11_hidrologia.csv").set_index("scenario")
pct_form = sum(float(str(x).split(";")[-1].split("%")[0]) for x in dina["Qué dice la base de datos"])
NOTA_PCT = (f"los % de tránsito del formulario suman {pct_form:.1f} %, el valor de la Tabla 11 para la penalización solo en el cauce principal "
            f"({t11.loc['H4m-main', 'F_8lines_cells_pct_mean (B5)']:.1f} %); no son los de la configuración final "
            f"({t11.loc['H4m', 'F_8lines_cells_pct_mean (B5)']:.1f} %, semillas 42-44)")
log.append("form percentages: " + NOTA_PCT)
d = dina.set_index("line")
rows = []
for _, r in sel.sort_values("line").iterrows():
    q = d.loc[r.line]
    seen = str(q.get("Qué se ve en la imagen", "")).strip()
    opinion = str(q.get("Qué es, en tu opinión", "")).strip()
    conf = str(q.get("Confianza", "")).strip()
    same = opinion.rstrip(". ").lower() == seen.rstrip(". ").lower()     # the opinion column repeats the observation
    if not opinion or opinion.lower().startswith("no se puede") or same:
        estado = "no verificable"
    elif conf.lower() == "alta":
        estado = "verificado"
    else:
        estado = f"interpretación tentativa (confianza {conf.lower()})"
    rows.append({"orden_Dina": int(q["N°"]), "HYRIV_ID": int(r.HYRIV_ID), "codigo_linea_formulario_Dina": int(q["linea_Dina"]), "linea_modelo_C87": int(r.line), "registro_shapefile_FID": int(r.shp_fid),
                 "orden_Strahler_oficial": int(r.ORD_STRA) if "ORD_STRA" in r else None, "longitud_km_oficial": float(r.LENGTH_KM),
                 "texto_base_datos_formulario": str(q.get("Qué dice la base de datos", "")).strip(),
                 "nota_porcentaje_formulario": NOTA_PCT,
                 "fuente_de_observacion": str(q.get("De dónde lo sabes", "")).strip(),
                 "fecha_de_observacion": "formulario devuelto por Dina el 05/11/2026 (fecha de la imagen no registrada)",
                 "observacion_visible": seen, "interpretacion_tentativa": opinion, "confianza": conf,
                 "nota_revision": ("la columna 'Qué es, en tu opinión' repite la observación: no hay interpretación de la naturaleza de la línea" if same else ""),
                 "uso_como_camino": str(q.get("¿Se pudo usar como camino?", "")).strip(), "comentario": str(q.get("Comentario", "")).strip(),
                 "celdas_atravesadas": ncell[r.line], "trafico_semilla42_pct": round(mass[r.line], 3),
                 "estado": estado,
                 "decision_heuristica": "incluido en la configuración final; retirado solo en la sensibilidad H4m-excl8 (Tabla 11)"})
T2 = pd.DataFrame(rows).sort_values("orden_Dina")
T2.to_csv(OUT / "S2_OCHO_TRAMOS.csv", index=False)

# 5. layers (no archaeological sites) and map
ca = pd.read_csv("out/final/informe/conteo_atractores_8_lineas.csv")
s0 = ca[ca.set == "S0 (H4m)"]; new = ca[ca.set != "S0 (H4m)"]
log.append(f"attractors: {len(s0)} in H4m; {int(s0.disappears_in_excl8.sum())} disappear and {len(new)} appear without the eight reaches; "
           f"{len(s0) - int(s0.disappears_in_excl8.sum()) + len(new)} remain; {int(s0.within_1km.sum())} of the 25 lie within 1 km of the reaches")
g8 = sel[["line", "HYRIV_ID", "geometry"]].merge(T2[["linea_modelo_C87", "orden_Dina", "confianza", "observacion_visible", "estado"]],
                                                left_on="line", right_on="linea_modelo_C87").drop(columns="linea_modelo_C87")
(OUT / "S2_ocho_tramos.gpkg").unlink(missing_ok=True)
g8.to_file(OUT / "S2_ocho_tramos.gpkg", layer="ocho_tramos_con_revision")
unc = g8[["HYRIV_ID", "orden_Dina", "estado", "confianza", "geometry"]].copy()      # separate uncertainty layer
unc["fecha_imagen"] = "no registrada"; unc["fuente"] = "Google Earth (revisión de Dina, formulario del 05/11/2026)"
unc["nota"] = "revisión remota; sin verificación de campo; no indica uso prehispánico"
unc.to_file(OUT / "S2_ocho_tramos.gpkg", layer="incertidumbre_revision")
gpd.GeoDataFrame(s0, geometry=gpd.points_from_xy(s0.X, s0.Y), crs=32718).to_file(OUT / "S2_ocho_tramos.gpkg", layer="atractores_H4m")
excl = pd.concat([s0[~s0.disappears_in_excl8.astype(bool)], new])
gpd.GeoDataFrame(excl, geometry=gpd.points_from_xy(excl.X, excl.Y), crs=32718).to_file(OUT / "S2_ocho_tramos.gpkg", layer="atractores_sin_ocho_tramos")

TXT, MUTED, BLUE, ORANGE, SURF = "#0b0b0b", "#898781", "#2a78d6", "#eb6834", "#fcfcfb"
V = st.valid_mask
xs, ys = sorted(lat["x_to_ix"]), sorted(lat["y_to_iy"])
ext = (min(xs) - 125, max(xs) + 125, min(ys) - 125, max(ys) + 125)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.9), facecolor=SURF)
for k, ax in enumerate(axes):
    ax.set_facecolor(SURF)
    ax.imshow(np.where(V, 1.0, np.nan), origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=6, interpolation="nearest")
    (cl if k == 0 else cl[~cl.shp_fid.isin(EXCL_FID)]).plot(ax=ax, color=MUTED, lw=0.6)   # right panel: the eight reaches removed
    if k == 0:
        nv = g8.estado.str.startswith("no verificable")
        g8[nv].plot(ax=ax, color=ORANGE, lw=2.2, linestyle=(0, (3, 1.5)))
        g8[~nv].plot(ax=ax, color=ORANGE, lw=2.2)
        for _, r in g8.iterrows():
            p = r.geometry.interpolate(0.5, normalized=True); ax.annotate(str(int(r.orden_Dina)), (p.x, p.y), fontsize=7, color=TXT, xytext=(3, 3), textcoords="offset points")
        keep = s0[~s0.disappears_in_excl8.astype(bool)]; gone = s0[s0.disappears_in_excl8.astype(bool)]
        ax.scatter(keep.X, keep.Y, s=26, c=BLUE, edgecolor=SURF, lw=0.8, zorder=3, label="attractor kept in both")
        ax.scatter(gone.X, gone.Y, s=26, facecolor=SURF, edgecolor=BLUE, lw=1.4, zorder=3, label="attractor that disappears without the eight reaches")
        ax.set_title("Final configuration: 25 attractors; the eight reaches in orange (numbers = review order)", fontsize=8.5, loc="left")
    else:
        keep = s0[~s0.disappears_in_excl8.astype(bool)]
        ax.scatter(keep.X, keep.Y, s=26, c=BLUE, edgecolor=SURF, lw=0.8, zorder=3, label="attractor kept in both")
        ax.scatter(new.X, new.Y, s=34, marker="^", c=BLUE, edgecolor=SURF, lw=0.8, zorder=3, label="attractor that appears without the eight reaches")
        assert len(excl) == int(t11.loc["H4m-excl8", "attractors_mean (B5)"])
        ax.set_title(f"Without the eight reaches: {len(excl)} attractors and {int(t11.loc['H4m-excl8', 'pairs_mean (B5)'])} pairs (Table 11)", fontsize=8.5, loc="left")
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")
    x0, y0 = ext[0], ext[2] - 1100                                     # scale bar below the map frame
    ax.plot([x0, x0 + 5000], [y0, y0], color=TXT, lw=2, clip_on=False); ax.text(x0 + 5400, y0, "5 km", va="center", fontsize=7)
H = [plt.Line2D([], [], color=ORANGE, lw=2.2, label="reach with a tentative interpretation"),
     plt.Line2D([], [], color=ORANGE, lw=2.2, ls=(0, (3, 1.5)), label="reach not verifiable"),
     plt.Line2D([], [], marker="o", ls="none", color=BLUE, mec=SURF, ms=6, label="attractor in both scenarios"),
     plt.Line2D([], [], marker="o", ls="none", mfc=SURF, mec=BLUE, mew=1.4, ms=6, label="attractor that disappears without the eight reaches"),
     plt.Line2D([], [], marker="^", ls="none", color=BLUE, mec=SURF, ms=7, label="attractor that appears without the eight reaches")]
fig.legend(handles=H, loc="lower center", ncol=3, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.04))
fig.text(0.01, 0.01, "Grey: the 12,826 analysis cells and the 87 HydroRIVERS reaches. EPSG:32718; north is up. Archaeological sites are not shown.", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.15, 1, 1))
fig.savefig(OUT / "S2_MAPA.pdf"); fig.savefig(OUT / "S2_MAPA.png", dpi=200)
(OUT / "S2_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
print("\n".join(log)); print(T2[["orden_Dina", "HYRIV_ID", "linea_modelo_C87", "observacion_visible", "interpretacion_tentativa", "confianza", "trafico_semilla42_pct", "estado"]].to_string(index=False, max_colwidth=45))
