#!/usr/bin/env python3
"""
final_sweep_eps_kappa.py  -  JCAA-80, recalculo final (ADR-002 rev. 5), bloque S6 de la matriz.

Barrido completo de epsilon y kappa con la configuracion final H4m:
    epsilon en {0, 0.04, ..., 0.40}  x  kappa en {0, 0.01, ..., 0.10}  x  semillas 1000-1009
    = 1,210 corridas (la misma malla y las mismas semillas que exp1_phase_scan.py).
Reemplaza a exp1_phase_scan.py + run_phase1_local.py, que usaban el kernel v1. Usa el kernel
ca_kernel_v2.py SIN CAMBIOS y el estado H4m de final_h4m_state.py (identico al de la Fase B4).

Integridad (validacion de ChatGPT, 05/11):
  - la prueba de entorno queda vinculada a la "firma" del paquete: sha256 del manifiesto, del
    script y de cada archivo de datos y codigo; --lanzar exige una prueba OK con la misma firma;
  - cada JSON lleva la firma y se valida por contenido (epsilon, kappa, semilla, campos
    obligatorios, valores finitos y en rango, 25 atractores y 90 pares);
  - al reanudar solo se corren los trabajos faltantes o invalidos; los invalidos se mueven a
    resultados_barrido_H4m/invalidos/ (no se borran);
  - --resumen solo acepta exactamente las 1,210 combinaciones esperadas, validas y sin duplicados.

Uso (Anaconda Prompt, dentro de esta carpeta):
    python final_sweep_eps_kappa.py --prueba              # unos minutos: verifica el entorno
    python final_sweep_eps_kappa.py --lanzar              # varias horas: SOLO cuando se avise
    python final_sweep_eps_kappa.py --lanzar --workers 8  # con otro numero de procesos
    python final_sweep_eps_kappa.py --verificar           # cuenta validos, invalidos y faltantes
    python final_sweep_eps_kappa.py --resumen             # al terminar: resumen y zip para enviar
"""
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import shutil
import sys
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

CSV = HERE / "datos" / "grilla_con_W_hidrica.csv"
RIV = HERE / "datos" / "celdas_fluviales_1020.csv"
ZON = HERE / "datos" / "zonas_P2_n20.csv"
REF = HERE / "referencia" / "H4m_s42_referencia.npz"
HUELLA = HERE / "referencia" / "huella_H4m.json"
MANIFEST = HERE / "MANIFEST.sha256"
PRUEBA_OK = HERE / "prueba_OK.json"
OUT = HERE / "resultados_barrido_H4m"
N_TOTAL = 1210
N_ATTR, N_PAIRS = 25, 90
REQUIRED = ("eps", "kappa", "seed", "lcc_size", "n_threshold", "m1", "m2", "m3", "auc", "n_clusters",
            "F_norm_sum_check", "rho20_P2", "n_attractors", "n_pairs", "runtime_s", "firma_paquete")


def eps_values():
    import numpy as np
    return [float(f"{v:.4f}") for v in np.linspace(0.0, 0.40, 11)]   # como phase1_jobs.txt


def kappa_values():
    import numpy as np
    return [float(f"{v:.4f}") for v in np.linspace(0.0, 0.10, 11)]


SEEDS = list(range(1000, 1010))


def job_name(e, k, s):
    return f"e{e:.3f}_k{k:.3f}_s{s}.json"


def all_jobs(out_dir):
    jobs = [(e, k, s, str(out_dir / job_name(e, k, s))) for e in eps_values() for k in kappa_values() for s in SEEDS]
    assert len(jobs) == N_TOTAL and len({j[3] for j in jobs}) == N_TOTAL
    return jobs


# ----------------------------------------------------------------------------- firma
def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check_manifest():
    """Los archivos de codigo y datos deben ser exactamente los del paquete."""
    bad = []
    for line in MANIFEST.read_text().splitlines():
        h, name = line.split(maxsplit=1)
        name = name.strip().lstrip("*").replace("./", "", 1)
        if name.startswith("LEEME"):
            continue
        p = HERE / name
        if not p.exists() or sha256(p) != h:
            bad.append(name)
    return bad


def firma():
    """Huella del paquete: manifiesto + script + cada archivo de codigo y datos."""
    parts = {"MANIFEST.sha256": sha256(MANIFEST), "final_sweep_eps_kappa.py": sha256(Path(__file__))}
    for p in (HERE / "ca_kernel_v2.py", HERE / "final_h4m_state.py", CSV, RIV, ZON, REF, HUELLA):
        parts[p.relative_to(HERE).as_posix()] = sha256(p)
    total = hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()
    return total, parts


# ----------------------------------------------------------------------------- utilidades
def versions():
    import numpy, scipy, pandas, networkx
    return {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
            "pandas": pandas.__version__, "networkx": networkx.__version__,
            "sistema": f"{platform.system()} {platform.release()}", "procesador": platform.processor(),
            "nucleos_logicos": os.cpu_count()}


def build_state():
    import final_h4m_state as fs
    return fs.build_state_h4m(str(CSV), str(RIV))


def zone_index(state):
    import numpy as np, pandas as pd
    z = pd.read_csv(ZON)
    lat = state.lattice
    iy = np.array([lat["y_to_iy"][y] for y in z.Y]); ix = np.array([lat["x_to_ix"][x] for x in z.X])
    return iy, ix, z.zone.to_numpy()


def rho20(L, state, zidx):
    import pandas as pd
    from scipy.stats import spearmanr
    iy, ix, zone = zidx
    d = pd.DataFrame({"z": zone, "L": L[iy, ix], "s": state.site_mask[iy, ix].astype(int)}).groupby("z")
    return float(spearmanr(d.L.mean(), d.s.sum())[0])


def fmt(sec):
    sec = int(sec); h, m = sec // 3600, (sec % 3600) // 60
    return f"{h}h{m:02d}m" if h else f"{m}m{sec % 60:02d}s"


def validate(path, e, k, s, sig):
    """Devuelve (True, '') si el JSON es una corrida valida de (e, k, s) con esta firma."""
    try:
        r = json.loads(Path(path).read_text())
    except Exception as exc:
        return False, f"no se puede leer ({exc.__class__.__name__})"
    miss = [c for c in REQUIRED if c not in r]
    if miss:
        return False, "faltan campos: " + ", ".join(miss)
    if r["firma_paquete"] != sig:
        return False, "firma del paquete distinta"
    if not all(isinstance(r[c], (int, float)) and not isinstance(r[c], bool) for c in ("eps", "kappa", "seed")):
        return False, "epsilon/kappa/semilla no son numeros"
    if abs(float(r["eps"]) - e) > 1e-12 or abs(float(r["kappa"]) - k) > 1e-12 or float(r["seed"]) != s:
        return False, f"epsilon/kappa/semilla no coinciden ({r['eps']}, {r['kappa']}, {r['seed']})"
    f = {}
    for c in REQUIRED[:-1]:
        if isinstance(r[c], bool) or not isinstance(r[c], (int, float)):
            return False, f"{c} no es numero"
        f[c] = float(r[c])
        if not math.isfinite(f[c]):
            return False, f"{c} no es finito"
    checks = [(0 <= f["auc"] <= 1, "auc fuera de [0,1]"), (0 <= f["m1"] <= 1, "m1 fuera de rango"),
              (0 <= f["m2"] <= 1, "m2 fuera de rango"), (0 < f["m3"] <= 1, "m3 fuera de rango"),
              (-1 <= f["rho20_P2"] <= 1, "rho20 fuera de rango"), (f["lcc_size"] >= 1, "lcc_size < 1"),
              (f["n_threshold"] >= 1, "n_threshold < 1"), (f["n_attractors"] == N_ATTR, "atractores != 25"),
              (f["n_pairs"] == N_PAIRS, "pares != 90"), (abs(f["F_norm_sum_check"] - 1) < 1e-9, "F no suma 1"),
              (f["runtime_s"] > 0, "runtime_s <= 0")]
    for ok, msg in checks:
        if not ok:
            return False, msg
    return True, ""


def classify(out_dir, sig):
    """Clasifica las 1,210 corridas esperadas y los archivos que sobran."""
    jobs = all_jobs(out_dir)
    valid, invalid, missing = [], [], []
    for e, k, s, p in jobs:
        if not Path(p).exists():
            missing.append((e, k, s, p))
        else:
            ok, why = validate(p, e, k, s, sig)
            (valid if ok else invalid).append((e, k, s, p, why) if not ok else (e, k, s, p))
    expected = {Path(j[3]).name for j in jobs}
    extra = sorted(f.name for f in out_dir.glob("*.json") if f.name not in expected) if out_dir.exists() else []
    return jobs, valid, invalid, missing, extra


# ----------------------------------------------------------------------------- prueba
def prueba():
    import numpy as np
    import ca_kernel_v2 as ca
    print("== PRUEBA DE ENTORNO (no produce resultados del articulo) ==")
    v = versions()
    for k, val in v.items():
        print(f"   {k:16s}: {val}")
    fallas = []
    bad = check_manifest()
    print(f"[1] Archivos del paquete: {'OK' if not bad else 'DISTINTOS: ' + ', '.join(bad)}")
    if bad:
        fallas.append("archivos distintos al paquete")
    sig, parts = firma()
    t0 = time.time(); st = build_state(); t_state = time.time() - t0
    hu = json.loads(HUELLA.read_text())
    same_attr = [list(map(int, a)) for a in st.attractors] == hu["attractors"]
    same_pairs = [list(map(int, p)) for p in st.pair_candidates] == hu["pairs"]
    ref = np.load(REF)
    V = st.valid_mask
    dE = float(np.abs(st.E_grid - ref["E"])[V].max()); dc = float(np.abs(st.c_grid - ref["c"])[V].max())
    dT = float(np.abs(st.T_grid - ref["T"])[V].max())
    ok2 = same_attr and same_pairs and max(dE, dc, dT) <= 1e-12
    print(f"[2] Estado H4m: {len(st.attractors)} atractores, {len(st.pair_candidates)} pares; "
          f"iguales a la referencia: atractores {same_attr}, pares {same_pairs}; "
          f"dif. max E' {dE:.1e}, c' {dc:.1e}, T {dT:.1e}  -> {'OK' if ok2 else 'FALLA'}  ({t_state:.1f} s)")
    if not ok2:
        fallas.append("estado H4m distinto")
    print("[3] Corrida de referencia: semilla 42, eps 0.15, kappa 0.03 (tarda unos minutos)...", flush=True)
    t0 = time.time(); r = ca.run_one(st, eps=0.15, kappa=0.03, seed=42, return_fields=True); t_run = time.time() - t0
    dF = float(np.abs(r["F_norm"] - ref["F"])[V].max()); dL = float(np.abs(r["L"] - ref["L"])[V].max())
    dA = abs(r["auc"] - float(ref["auc"]))
    rho = rho20(r["L"], st, zone_index(st)); dR = abs(rho - float(ref["rho20_P2"]))
    ok3 = max(dF, dL, dA, dR) <= 1e-12
    print(f"    AUC {r['auc']:.6f} (ref. {float(ref['auc']):.6f}); rho20 P2 {rho:.6f} (ref. {float(ref['rho20_P2']):.6f}); "
          f"dif. max F {dF:.1e}, L {dL:.1e}  -> {'OK' if ok3 else 'FALLA'}  ({fmt(t_run)})")
    if not ok3:
        fallas.append("la corrida de referencia no coincide")
    w = max(1, (os.cpu_count() or 2) - 1)
    print(f"[4] Tiempo estimado del barrido con {w} procesos: entre {fmt(N_TOTAL * t_run / w)} y "
          f"{fmt(N_TOTAL * t_run / max(1, w // 2))} (depende de la temperatura y de los nucleos fisicos).")
    res = {"fecha": time.strftime("%Y-%m-%d %H:%M:%S"), "versiones": v, "firma_paquete": sig, "archivos": parts,
           "segundos_por_corrida": round(t_run, 1),
           "dif_max": {"E": dE, "c": dc, "T": dT, "F": dF, "L": dL, "AUC": dA, "rho20": dR},
           "auc_s42": r["auc"], "rho20_s42": rho, "fallas": fallas}
    if fallas:
        print("\nPRUEBA FALLO: " + "; ".join(fallas))
        print("No lances el barrido. Envia a Claude todo lo que aparece en esta pantalla.")
        (HERE / "prueba_FALLO.json").write_text(json.dumps(res, indent=1))
        if PRUEBA_OK.exists():
            PRUEBA_OK.unlink()
        return 1
    PRUEBA_OK.write_text(json.dumps(res, indent=1))
    print(f"\nPRUEBA OK (firma del paquete {sig[:12]}). Se guardo prueba_OK.json. Envialo a Claude.")
    print("No lances el barrido (--lanzar) hasta que Claude te avise.")
    return 0


def prueba_vigente():
    """La prueba OK debe existir, no tener fallas y corresponder a la firma actual del paquete."""
    if not PRUEBA_OK.exists():
        return None, "no hay prueba_OK.json: corre primero --prueba"
    pr = json.loads(PRUEBA_OK.read_text())
    sig, _ = firma()
    if pr.get("fallas"):
        return None, "la prueba registrada tiene fallas"
    if pr.get("firma_paquete") != sig:
        return None, "el paquete cambio desde la prueba (firma distinta): vuelve a correr --prueba"
    bad = check_manifest()
    if bad:
        return None, "archivos distintos al manifiesto: " + ", ".join(bad)
    return sig, ""


# ----------------------------------------------------------------------------- barrido
_ST = None
_Z = None
_SIG = None


def _init_worker(sig):
    global _ST, _Z, _SIG
    sys.path.insert(0, str(HERE))
    _ST = build_state()
    _Z = zone_index(_ST)
    _SIG = sig


def _job(job):
    import ca_kernel_v2 as ca
    e, k, s, out = job
    out = Path(out)
    try:
        t0 = time.time()
        r = ca.run_one(_ST, eps=e, kappa=k, seed=s, return_fields=True)
        res = {key: r[key] for key in ("eps", "kappa", "seed", "lcc_size", "n_threshold", "m1", "m2", "m3",
                                         "auc", "n_clusters", "F_norm_sum_check")}
        res["rho20_P2"] = rho20(r["L"], _ST, _Z)
        res["n_attractors"] = len(_ST.attractors); res["n_pairs"] = len(_ST.pair_candidates)
        res["runtime_s"] = round(time.time() - t0, 2)
        res["firma_paquete"] = _SIG
        tmp = out.with_suffix(".tmp")
        tmp.write_text(json.dumps(res, indent=1))
        os.replace(tmp, out)
        ok, why = validate(out, e, k, s, _SIG)
        return ("ok" if ok else "invalido", out.name, res["auc"] if ok else why)
    except Exception as exc:  # se informa y se sigue con las demas
        return ("error", out.name, repr(exc))


def quarantine(items, out_dir):
    q = out_dir / "invalidos"
    q.mkdir(exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    for e, k, s, p, why in items:
        shutil.move(p, q / f"{Path(p).stem}__{stamp}.json")
        print(f"   invalido -> invalidos/: {Path(p).name} ({why})")


def lanzar(workers, out_dir, solo=None):
    sig, why = prueba_vigente()
    if sig is None:
        print("No se puede lanzar: " + why)
        return 1
    out_dir.mkdir(exist_ok=True)
    jobs, valid, invalid, missing, extra = classify(out_dir, sig)
    if extra:
        print(f"Aviso: {len(extra)} archivos que no corresponden a ninguna corrida esperada (no se usan): {extra[:5]}")
    if invalid:
        quarantine(invalid, out_dir)
    todo = [j[:4] for j in invalid] + missing
    todo.sort(key=lambda j: (j[0], j[1], j[2]))
    if solo:
        todo = todo[:solo]
    w = workers or max(1, (os.cpu_count() or 2) - 1)
    print(f"Corridas: {N_TOTAL} en total; validas {len(valid)}; invalidas {len(invalid)}; faltan {len(missing)}. "
          f"Por hacer en esta tanda: {len(todo)}. Procesos: {w}.")
    if not todo:
        print("No hay nada pendiente. Siguiente paso: --resumen")
        return 0
    t0 = time.time(); n = 0; err = 0
    with mp.Pool(processes=w, initializer=_init_worker, initargs=(sig,)) as pool:
        try:
            for st, name, val in pool.imap_unordered(_job, todo, chunksize=1):
                n += 1
                if st != "ok":
                    err += 1
                el = time.time() - t0
                eta = (len(todo) - n) * el / n
                det = f"AUC={val:.4f}" if st == "ok" else str(val)
                print(f"[{n:4d}/{len(todo)}] {st:8s} {name}  {det}  transcurrido {fmt(el)}  faltan ~{fmt(eta)}", flush=True)
        except KeyboardInterrupt:
            print("\nInterrumpido. Lo hecho queda guardado; vuelve a escribir el mismo comando para seguir.")
            pool.terminate(); pool.join()
            return 130
    _, valid, invalid, missing, _ = classify(out_dir, sig)
    print(f"\nTanda terminada en {fmt(time.time() - t0)}. Estado: validas {len(valid)}, invalidas {len(invalid)}, "
          f"faltan {len(missing)} (errores en esta tanda: {err}).")
    if len(valid) == N_TOTAL:
        print("Completo. Siguiente paso: python final_sweep_eps_kappa.py --resumen")
    else:
        print("Incompleto: vuelve a escribir --lanzar; solo se correran las faltantes o invalidas.")
    return 0 if err == 0 else 1


def verificar(out_dir):
    sig, why = prueba_vigente()
    if sig is None:
        print("Sin prueba vigente: " + why)
        return 1
    _, valid, invalid, missing, extra = classify(out_dir, sig)
    print(f"Validas {len(valid)} de {N_TOTAL}; invalidas {len(invalid)}; faltan {len(missing)}; archivos sobrantes {len(extra)}.")
    for e, k, s, p, why in invalid[:20]:
        print(f"   invalido: {Path(p).name} ({why})")
    return 0 if len(valid) == N_TOTAL and not extra else 1


def resumen(out_dir, parcial=False):
    import pandas as pd
    sig, why = prueba_vigente()
    if sig is None:
        print("Sin prueba vigente: " + why)
        return 1
    jobs, valid, invalid, missing, extra = classify(out_dir, sig)
    rows = [json.loads(Path(p).read_text()) for e, k, s, p in valid]
    keys = [(round(r["eps"], 4), round(r["kappa"], 4), int(r["seed"])) for r in rows]
    expected = {(round(e, 4), round(k, 4), s) for e, k, s, _ in jobs}
    complete = (len(valid) == N_TOTAL and not invalid and not missing and not extra
                and len(set(keys)) == len(keys) == N_TOTAL and set(keys) == expected)
    print(f"Validas {len(valid)} de {N_TOTAL}; invalidas {len(invalid)}; faltan {len(missing)}; sobrantes {len(extra)}; "
          f"combinaciones unicas {len(set(keys))}.")
    if not complete and not parcial:
        print("El resumen solo se acepta con las 1,210 combinaciones exactas y validas. Vuelve a --lanzar para completar.")
        return 1
    df = pd.DataFrame(rows)
    g = df.groupby(["eps", "kappa"]).agg(n=("auc", "size"), auc_media=("auc", "mean"), auc_de=("auc", "std"),
                                         rho20_media=("rho20_P2", "mean"), rho20_de=("rho20_P2", "std"),
                                         m1_media=("m1", "mean"), m2_media=("m2", "mean"), m3_media=("m3", "mean"),
                                         lcc_media=("lcc_size", "mean"), segundos_media=("runtime_s", "mean")).reset_index()
    if complete:
        assert len(g) == 121 and (g.n == 10).all()
    g.to_csv(HERE / "resumen_barrido_H4m.csv", index=False)
    print(f"Configuraciones: {len(g)}; AUC media por configuracion entre {g.auc_media.min():.4f} y {g.auc_media.max():.4f}.")
    zp = HERE / ("para_enviar_barrido_H4m.zip" if complete else "para_enviar_barrido_H4m_PARCIAL.zip")
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as zf:
        for e, k, s, p in valid:
            zf.write(p, f"resultados_barrido_H4m/{Path(p).name}")
        for extra_f in ("resumen_barrido_H4m.csv", "prueba_OK.json", "MANIFEST.sha256"):
            if (HERE / extra_f).exists():
                zf.write(HERE / extra_f, extra_f)
    print(f"Listo: {zp.name}. Ese es el archivo que se envia a Claude." + ("" if complete else " (PARCIAL: no es el resultado final)"))
    return 0


def main():
    ap = argparse.ArgumentParser(description="Barrido eps/kappa H4m (JCAA-80, bloque S6)")
    ap.add_argument("--prueba", action="store_true", help="verifica el entorno (unos minutos)")
    ap.add_argument("--lanzar", action="store_true", help="corre las corridas faltantes o invalidas (solo cuando se avise)")
    ap.add_argument("--verificar", action="store_true", help="cuenta corridas validas, invalidas y faltantes")
    ap.add_argument("--resumen", action="store_true", help="resume y crea el zip para enviar")
    ap.add_argument("--workers", type=int, default=0, help="procesos en paralelo (por defecto: nucleos logicos - 1)")
    ap.add_argument("--solo", type=int, default=0, help=argparse.SUPPRESS)          # pruebas del script
    ap.add_argument("--salida", default=str(OUT), help=argparse.SUPPRESS)          # pruebas del script
    ap.add_argument("--parcial", action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args()
    out_dir = Path(a.salida)
    if a.prueba:
        return prueba()
    if a.lanzar:
        return lanzar(a.workers, out_dir, a.solo or None)
    if a.verificar:
        return verificar(out_dir)
    if a.resumen:
        return resumen(out_dir, a.parcial)
    ap.print_help()
    return 0


if __name__ == "__main__":
    mp.freeze_support()          # necesario en Windows
    sys.exit(main())
