"""Restricted files of the repository: the analysis grid columns and files that give the location of archaeological sites.

The public repository keeps the analysis grid without its three site columns (data/grilla_con_W_hidrica_publica.csv)
and leaves out six other files (RESTRICTED_MANIFEST.sha256). With an authorised bundle of the restricted data,
tools/restore_restricted.py rebuilds the original files byte for byte; verify_repository.py then runs in full.
Standard library only."""
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID = "data/grilla_con_W_hidrica.csv"
GRID_PUBLIC = "data/grilla_con_W_hidrica_publica.csv"
GRID_COPIES = [GRID, "out/b6/ryzen/barrido_H4m/datos/grilla_con_W_hidrica.csv"]   # identical files
SITE_FILE = "grid_site_columns.csv"                                             # name inside the bundle
SITE_COLUMNS = [b"distancia_sitios_nuevos", b"presencia_sitio_nuevo", b"distancia_mincul"]
ORIGINAL_HEADER = (b"X,Y,distancia_sitios_nuevos,presencia_sitio_nuevo,distancia_mincul,altitud,pendiente,dist_river_m,"
                   b"dist_river_norm,p_norm,h_norm,W_river,W_corridor").split(b",")
WHOLE_FILES = ["data/ref/fields_canonical_v2_for_ppm.csv", "out/b6/ryzen/JCAA80_barrido_H4m_Ryzen.zip",
               "out/final/ppm/fields_H4m_for_ppm.csv", "out/final/qgis/sitios_84.gpkg", "out/rowmajor_sites.npy",
               "out/state_v2_repro.pkl"]


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256(path):
    return sha256_bytes(Path(path).read_bytes())


def restricted_manifest(root=ROOT):
    rows = {}
    for line in (Path(root) / "RESTRICTED_MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            h, p = line.split(maxsplit=1)
            rows[p.strip()] = h
    return rows


def _cells(line):
    body = line.rstrip(b"\r\n")
    return body.split(b","), line[len(body):]


def split_grid(data):
    """Full grid bytes -> (public grid bytes, site-column bytes). Cells are kept as text, so nothing is re-formatted."""
    lines = data.splitlines(keepends=True)
    head, _ = _cells(lines[0])
    assert head == ORIGINAL_HEADER, head
    idx = [head.index(c) for c in SITE_COLUMNS]
    keep = [i for i in range(len(head)) if i not in idx]
    pub, site = [], []
    for ln in lines:
        cells, end = _cells(ln)
        assert len(cells) == len(head)
        pub.append(b",".join(cells[i] for i in keep) + end)
        site.append(b",".join([cells[0], cells[1]] + [cells[i] for i in idx]) + end)
    return b"".join(pub), b"".join(site)


def merge_grid(pub, site):
    """Inverse of split_grid: rebuilds the full grid in its original column order."""
    P, S = pub.splitlines(keepends=True), site.splitlines(keepends=True)
    assert len(P) == len(S), "the public grid and the site columns have different numbers of rows"
    hp, _ = _cells(P[0]); hs, _ = _cells(S[0])
    out = []
    for lp, ls in zip(P, S):
        cp, end = _cells(lp); cs, _ = _cells(ls)
        assert cp[:2] == cs[:2], f"rows do not match: {cp[:2]} vs {cs[:2]}"
        row = dict(zip(hp, cp)); row.update(zip(hs, cs))
        out.append(b",".join(row[c] for c in ORIGINAL_HEADER) + end)
    return b"".join(out)
