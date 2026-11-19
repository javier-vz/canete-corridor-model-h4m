"""JCAA-80 final recompute (matrix validated 05/11): runner for blocks S1-S5, S7-S9 and the v2
controls (CTRL). Frozen; the kernel is not modified. S6 (eps/kappa sweep) runs on the Ryzen with
final_sweep_eps_kappa.py.

Integrity (same rules as the sweep, validation point B):
  - each job writes <block>/<id>.json (+ .npz with F and L, + _paths.pkl for S1) atomically;
  - the JSON records the exact job specification, the signature of code and data (firma), the
    state fingerprint, the metrics, and the sha256 of its .npz / .pkl;
  - on resume every expected job is validated by content; only missing or invalid jobs are run;
    invalid outputs are moved to <block>/invalidos/ (never deleted).

Usage (from the project root):
  python code/final_runs.py --bloques S1 S2 S7 S8 S3 S4 S5 S9 --workers 2
  python code/final_runs.py --bloques CTRL --workers 2          # v2 controls (pre-tests)
  python code/final_runs.py --verificar --bloques S1 S2 ...
"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import argparse, json, math, multiprocessing as mp, pickle, shutil, time
from pathlib import Path
import numpy as np
import pandas as pd

import ca_kernel_v2 as ca
import final_common as fc
import final_h4m_state as fs

EPS, KAP = 0.15, 0.03
NOM = dict(rules="H4m", a=list(fs.A_NOMINAL), b=list(fs.B_NOMINAL), q=fs.Q_NOMINAL, T=list(fs.T_NOMINAL), slope=None)
SLOPE_VARIANTS = ("S2_clip_p1_p99", "S3_log1p_minmax", "S4_FAO_classes", "S1agg_cellmean_minmax")
ALL_BLOCKS = ("S1", "S2", "S3", "S4", "S5", "S7", "S8", "S9", "CTRL")
CA_KEYS = ("auc", "lcc_size", "n_threshold", "m1", "m2", "m3", "n_clusters", "F_norm_sum_check")
FETE_KEYS = ("auc", "lcc_size", "n_threshold", "m1", "m2", "m3", "n_successful_paths", "n_random_pairs")


def cfg(**kw):
    c = json.loads(json.dumps(NOM)); c.update(kw)
    return json.loads(json.dumps(c))        # canonical JSON types (lists, floats)


def safe(label):
    return label.replace("=", "").replace("+", "p").replace("-", "m")


def jobs():
    DC = pd.read_csv("out/b6/design_counts_H4m.csv")
    dc = {(r.family, r.config): (int(r.attractors), int(r.pairs)) for r in DC.itertuples()}
    J = []

    def add(block, jid, kind, conf, seed, eps=EPS, kappa=KAP, n_pairs=None, cost=None, fields=False, paths=False, expect=None):
        J.append(dict(block=block, id=jid, kind=kind, cfg=conf, seed=int(seed), eps=float(eps), kappa=float(kappa),
                      n_pairs=n_pairs, cost=cost, fields=fields, paths=paths, expect=list(expect) if expect else None))
    add("S1", "canon_s42", "ca", cfg(), 42, fields=True, paths=True, expect=(25, 90))
    add("S1", "det_s42", "ca", cfg(), 42, eps=0.0, kappa=0.0, fields=True, paths=True, expect=(25, 90))
    for s in range(2000, 2030):
        add("S2", f"canon_s{s}", "ca", cfg(), s, fields=True, expect=(25, 90))
    for lab, a, b in fs.weight_configs():
        for s in (3000, 3001, 3002):
            add("S3", f"{safe(lab)}_s{s}", "ca", cfg(a=list(a), b=list(b)), s, fields=True, expect=dc[("weights", lab)])
    for lab, r, sv, w in fs.T_configs():
        for s in (5000, 5001, 5002):
            add("S4", f"{safe(lab)}_s{s}", "ca", cfg(T=[r, sv, w]), s, fields=True, expect=dc[("T", lab)])
    for q in fs.QUANTILES:
        for s in (4000, 4001, 4002):
            add("S5", f"q{q}_s{s}", "ca", cfg(q=q), s, fields=True, expect=dc[("quantile", f"q={q}")])
    for n in (50, 90, 200, 500, 1000):
        for s in range(2000, 2030):
            add("S7", f"fete_c_n{n}_s{s}", "fete", cfg(), s, n_pairs=n, cost="c")
    for s in range(2000, 2030):
        add("S8", f"fete_1mE_n90_s{s}", "fete", cfg(), s, n_pairs=90, cost="1-E")
    for v in SLOPE_VARIANTS:
        add("S9", f"{v}_s42", "ca", cfg(slope=v), 42, fields=True)
    a_m20 = fs.perturb(fs.A_NOMINAL, 2, -1)                    # a_w -20 % (Table 3 row)
    for s in (3000, 3001, 3002):
        add("CTRL", f"v2_a_w_m20_s{s}", "ca", cfg(rules="v2", a=list(a_m20)), s)
    for s in (5000, 5001, 5002):
        add("CTRL", f"v2_r0.80_s{s}", "ca", cfg(rules="v2", T=[0.80, 0.20, 5]), s, expect=(24, None))
    for s in (4000, 4001, 4002):
        add("CTRL", f"v2_q0.975_s{s}", "ca", cfg(rules="v2", q=0.975), s, expect=(22, 78))
    assert len({(j["block"], j["id"]) for j in J}) == len(J)
    return J


# ----------------------------------------------------------------------------- outputs and validation
def paths_for(job, outdir):
    d = Path(outdir) / job["block"]
    return d / f"{job['id']}.json", d / f"{job['id']}.npz", d / f"{job['id']}_paths.pkl"


def validate(job, sig, outdir):
    pj, pz, pp = paths_for(job, outdir)
    if not pj.exists():
        return None, "falta"
    try:
        r = json.loads(pj.read_text())
    except Exception as exc:
        return False, f"JSON ilegible ({exc.__class__.__name__})"
    if r.get("spec") != json.loads(json.dumps(job)):
        return False, "especificación distinta"
    if r.get("firma") != sig:
        return False, "firma de código/datos distinta"
    m = r.get("metrics", {})
    keys = CA_KEYS if job["kind"] == "ca" else FETE_KEYS
    for k in keys:
        v = m.get(k)
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            return False, f"métrica {k} ausente o no finita"
    if not (0 <= m["auc"] <= 1 and 0 <= m["m1"] <= 1 and 0 <= m["m2"] <= 1 and 0 < m["m3"] <= 1 and m["lcc_size"] >= 1):
        return False, "métricas fuera de rango"
    if job["kind"] == "ca" and abs(m["F_norm_sum_check"] - 1) > 1e-9:
        return False, "F no suma 1"
    if job["kind"] == "fete" and m["n_random_pairs"] != job["n_pairs"]:
        return False, "n_random_pairs distinto"
    if job["expect"]:
        ea, ep = job["expect"]
        if r.get("n_attractors") != ea or (ep is not None and r.get("n_pairs") != ep):
            return False, f"atractores/pares {r.get('n_attractors')}/{r.get('n_pairs')} (esperado {ea}/{ep})"
    if job["fields"]:
        if not pz.exists() or fc.sha256(pz) != r.get("npz_sha256"):
            return False, "npz ausente o distinto"
    if job["paths"]:
        if not pp.exists() or fc.sha256(pp) != r.get("paths_sha256"):
            return False, "rutas ausentes o distintas"
        if not (r.get("F_rebuild_maxdiff", 1) < 1e-12 and r.get("n_path_records", 0) == r.get("n_pairs")):
            return False, "F reconstruido o registros de ruta no cuadran"
    return True, ""


# ----------------------------------------------------------------------------- worker
_ORIG = ca.k_simple_paths_weighted
_REC = []


def _observer(G, source, target, n_paths=4):
    p = _ORIG(G, source, target, n_paths=n_paths)
    _REC.append(dict(source=source, target=target, paths=[np.asarray(x, dtype=np.int32) for x in p],
                     costs=np.array([ca.path_cost(G, x) for x in p], dtype=float)))
    return p


def traffic_from_records(st, records, beta_path=1.5):
    T = np.zeros_like(st.E_grid, dtype=float)
    for r in records:
        if not r["paths"]:
            continue
        c0 = r["costs"].min(); b = np.exp(-beta_path * (r["costs"] - c0)); b /= b.sum() + 1e-12
        (i1, j1), (i2, j2) = r["source"], r["target"]
        pw = np.exp(2.0 * st.E_grid[i1, j1]) * np.exp(2.0 * st.E_grid[i2, j2]) * np.exp(-3.0 * c0)
        for p, pp in zip(r["paths"], b):
            T[p[:, 0], p[:, 1]] += pw * pp
    return T / T.sum()


_SIG = None


def _init(sig):
    global _SIG
    _os.chdir(_ROOT)
    _SIG = sig


def run_job(job, outdir):
    pj, pz, pp = paths_for(job, outdir)
    st = fc.state_for(job["cfg"])
    t0 = time.time()
    out = dict(spec=job, firma=_SIG, state_hash=fc.state_key_hash(st), n_attractors=len(st.attractors),
               n_pairs=len(st.pair_candidates), fecha=time.strftime("%Y-%m-%d %H:%M:%S"))
    if job["kind"] == "ca":
        if job["paths"]:
            _REC.clear(); ca.k_simple_paths_weighted = _observer
        try:
            r = ca.run_one(st, eps=job["eps"], kappa=job["kappa"], seed=job["seed"], return_fields=True)
        finally:
            ca.k_simple_paths_weighted = _ORIG
        out["metrics"] = {k: (float(r[k]) if isinstance(r[k], float) else int(r[k])) for k in CA_KEYS}
        if job["paths"]:
            F_rec = traffic_from_records(st, _REC)
            out["F_rebuild_maxdiff"] = float(np.abs(F_rec - r["F_norm"]).max())
            out["n_path_records"] = len(_REC)
            out["n_path_alternatives"] = int(sum(len(x["paths"]) for x in _REC))
            fc.atomic_bytes(pp, pickle.dumps(list(_REC)))
            out["paths_sha256"] = fc.sha256(pp)
        if job["fields"]:
            fc.atomic_npz(pz, F=r["F_norm"], L=r["L"])
            out["npz_sha256"] = fc.sha256(pz)
    else:
        cost = None if job["cost"] == "c" else 1.0 - st.E_grid
        r = ca.run_fete(st, n_random_pairs=job["n_pairs"], seed=job["seed"], cost_surface=cost)
        out["metrics"] = {k: (float(r[k]) if isinstance(r[k], float) else int(r[k])) for k in FETE_KEYS}
    out["runtime_s"] = round(time.time() - t0, 2)
    fc.atomic_json(pj, out)
    return out


def _work(args):
    job, outdir = args
    try:
        out = run_job(job, outdir)
        ok, why = validate(job, _SIG, outdir)
        return job["block"], job["id"], ("ok" if ok else "invalido"), (f"AUC={out['metrics']['auc']:.4f} {out['runtime_s']:.0f}s" if ok else why)
    except Exception as exc:
        return job["block"], job["id"], "error", repr(exc)


def classify(J, sig, outdir):
    valid, invalid, missing = [], [], []
    for j in J:
        ok, why = validate(j, sig, outdir)
        if ok is None:
            missing.append(j)
        elif ok:
            valid.append(j)
        else:
            invalid.append((j, why))
    return valid, invalid, missing


def quarantine(invalid, outdir):
    stamp = time.strftime("%Y%m%d-%H%M%S")
    for j, why in invalid:
        q = Path(outdir) / j["block"] / "invalidos"; q.mkdir(parents=True, exist_ok=True)
        for p in paths_for(j, outdir):
            if p.exists():
                shutil.move(str(p), str(q / f"{p.stem}__{stamp}{p.suffix}"))
        print(f"   inválido -> {j['block']}/invalidos/: {j['id']} ({why})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bloques", nargs="+", default=[], choices=ALL_BLOCKS)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--verificar", action="store_true")
    ap.add_argument("--salida", default=str(fc.OUT / "runs"))
    ap.add_argument("--solo", type=int, default=0)
    a = ap.parse_args()
    sig, parts = fc.firma_calculo()
    J = [j for j in jobs() if j["block"] in a.bloques]
    valid, invalid, missing = classify(J, sig, a.salida)
    print(f"firma {sig[:12]} | bloques {a.bloques} | trabajos {len(J)}: válidos {len(valid)}, inválidos {len(invalid)}, faltan {len(missing)}", flush=True)
    if a.verificar:
        for j, why in invalid[:30]:
            print(f"   inválido: {j['block']}/{j['id']} ({why})")
        return 0 if len(valid) == len(J) else 1
    if invalid:
        quarantine(invalid, a.salida)
    todo = [j for j, _ in invalid] + missing
    # long jobs first so the two workers finish together
    order = {"S1": 0, "S9": 1, "S3": 2, "S4": 3, "S5": 4, "S2": 5, "CTRL": 6, "S8": 7, "S7": 8}
    todo.sort(key=lambda j: (order[j["block"]], -(j["n_pairs"] or 0)))
    if a.solo:
        todo = todo[:a.solo]
    if not todo:
        print("nada pendiente"); return 0
    t0 = time.time(); n = 0; bad = 0
    with mp.Pool(processes=a.workers, initializer=_init, initargs=(sig,)) as pool:
        for blk, jid, st, det in pool.imap_unordered(_work, [(j, a.salida) for j in todo], chunksize=1):
            n += 1; bad += st != "ok"
            el = time.time() - t0
            print(f"[{n:4d}/{len(todo)}] {st:8s} {blk}/{jid}  {det}  {el/60:.1f} min, faltan ~{(len(todo)-n)*el/n/60:.0f} min", flush=True)
    valid, invalid, missing = classify(J, sig, a.salida)
    print(f"FIN: válidos {len(valid)}/{len(J)}, inválidos {len(invalid)}, faltan {len(missing)}, errores {bad} ({(time.time()-t0)/60:.1f} min)")
    return 0 if len(valid) == len(J) else 1


if __name__ == "__main__":
    _sys.exit(main())
