"""Least-cost routing on a regular lattice with scipy.sparse (8 or 16 directions).
Edge weight = mean cost of the cells the move crosses x physical length (map units).
Orthogonal/diagonal: 0.5 (c_u + c_v) * {1, sqrt2} * res (the kernel rule x res).
Knight (1,2)/(2,1): mean of the 4 cells crossed * sqrt5 * res; forbidden if any crossed
cell is outside the domain. Shared by E2, the GRASS check and the route diagnostics."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from shapely.geometry import LineString

OFF8 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)]
OFFK = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def knight_cells(di, dj):
    if abs(di) == 1:
        return [(0, int(np.sign(dj))), (di, int(np.sign(dj)))]
    return [(int(np.sign(di)), 0), (int(np.sign(di)), dj)]


def build_graph(c, valid, res, knight):
    ny, nx = valid.shape; idx = -np.ones((ny, nx), dtype=np.int64); idx[valid] = np.arange(valid.sum())
    rows, cols, w = [], [], []
    for di, dj in OFF8 + (OFFK if knight else []):
        is_k = (abs(di), abs(dj)) in [(1, 2), (2, 1)]
        i0, i1 = max(0, -di), min(ny, ny - di); j0, j1 = max(0, -dj), min(nx, nx - dj)
        src = valid[i0:i1, j0:j1] & valid[i0 + di:i1 + di, j0 + dj:j1 + dj]
        cs, ct = c[i0:i1, j0:j1], c[i0 + di:i1 + di, j0 + dj:j1 + dj]
        if is_k:
            acc, ok = cs + ct, src.copy()
            for a, b in knight_cells(di, dj):
                ok &= valid[i0 + a:i1 + a, j0 + b:j1 + b]; acc = acc + c[i0 + a:i1 + a, j0 + b:j1 + b]
            src, wt = ok, acc / 4 * np.sqrt(5) * res
        else:
            wt = 0.5 * (cs + ct) * np.hypot(di, dj) * res
        ii, jj = np.nonzero(src)
        rows.append(idx[i0:i1, j0:j1][ii, jj]); cols.append(idx[i0 + di:i1 + di, j0 + dj:j1 + dj][ii, jj]); w.append(wt[ii, jj])
    n = int(valid.sum())
    G = coo_matrix((np.concatenate(w), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    return G, idx


def route(G, idx, s, t):
    d, pred = dijkstra(G, directed=True, indices=idx[s], return_predecessors=True)
    inv = np.argwhere(idx >= 0); path = [idx[t]]
    while path[-1] != idx[s]:
        path.append(pred[path[-1]])
    path = path[::-1]
    return [tuple(inv[p]) for p in path], float(d[idx[t]])


def polyline(cells, xs, ys):
    return LineString([(xs[j], ys[i]) for i, j in cells])


def sample_pts(line, step=25.0):
    n = max(2, int(line.length // step) + 1)
    return [line.interpolate(t, normalized=True) for t in np.linspace(0, 1, n)]


def mean_sep(a, b):
    pa, pb = sample_pts(a), sample_pts(b)
    return 0.5 * (np.mean([b.distance(p) for p in pa]) + np.mean([a.distance(p) for p in pb]))


def angular(cells):
    v = np.diff(np.array(cells), axis=0)
    ang = np.degrees(np.arctan2(v[:, 0], v[:, 1])) % 180
    seglen = np.hypot(v[:, 0], v[:, 1])
    on45 = np.isclose(ang % 45, 0) | np.isclose(ang % 45, 45)
    return float(seglen[on45].sum() / seglen.sum())


