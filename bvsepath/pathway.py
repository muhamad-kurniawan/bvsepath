import heapq
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .cube import read_cube

NEIGHBORS = [np.array(x, int) for x in itertools.product([-1, 0, 1], repeat=3) if x != (0, 0, 0)]


def _basins(E):
    shape = np.asarray(E.shape, int)
    best = E.copy()
    step = np.zeros(E.shape + (3,), np.int8)
    for off in NEIGHBORS:
        q = np.roll(E, shift=tuple(-off), axis=(0, 1, 2))
        mask = q < best
        best[mask] = q[mask]
        step[mask] = off
    coords = np.indices(E.shape, dtype=np.int32).transpose(1, 2, 3, 0)
    raw_next = coords + step
    wrapped = raw_next % shape
    nxt = np.ravel_multi_index((wrapped[..., 0], wrapped[..., 1], wrapped[..., 2]), E.shape).ravel()
    tile_step = ((raw_next - wrapped) // shape).reshape(-1, 3).astype(np.int16)
    n = E.size
    root = np.full(n, -1, np.int64)
    shift = np.zeros((n, 3), np.int16)
    for start in range(n):
        if root[start] >= 0:
            continue
        trail = []
        x = start
        while root[x] < 0:
            trail.append(x)
            y = int(nxt[x])
            if y == x:
                root[x] = x
                shift[x] = 0
                break
            x = y
        for node in reversed(trail):
            y = int(nxt[node])
            if y == node:
                root[node] = node
                shift[node] = 0
            else:
                root[node] = root[y]
                shift[node] = tile_step[node] + shift[y]
    minima_flat = np.where(root == np.arange(n))[0]
    remap = {int(r): i for i, r in enumerate(minima_flat)}
    labels = np.array([remap[int(r)] for r in root], np.int32).reshape(E.shape)
    min_grid = np.array(np.unravel_index(minima_flat, E.shape)).T
    min_energy = E.ravel()[minima_flat]
    return labels, shift.reshape(E.shape + (3,)), min_grid, min_energy


def _half_neighbors():
    out = []
    for off in NEIGHBORS:
        for x in off:
            if x != 0:
                if x > 0:
                    out.append(off)
                break
    return out


def _saddles(cube, labels, shifts, min_grid, min_energy, max_relative_eV):
    E = cube.values
    shape = cube.shape
    coords = np.indices(E.shape, dtype=np.int32).transpose(1, 2, 3, 0)
    global_min = float(E.min())
    best = {}
    for off in _half_neighbors():
        raw = coords + off
        wrapped = raw % shape
        edge_tile = ((raw - wrapped) // shape).astype(np.int16)
        other_label = np.roll(labels, shift=tuple(-off), axis=(0, 1, 2))
        other_E = np.roll(E, shift=tuple(-off), axis=(0, 1, 2))
        other_shift = np.roll(shifts, shift=tuple(-off), axis=(0, 1, 2))
        trans = other_shift + edge_tile - shifts
        mask = (labels != other_label) | ((labels == other_label) & np.any(trans != 0, axis=-1))
        for p in np.argwhere(mask):
            pt = tuple(p)
            i = int(labels[pt])
            j = int(other_label[pt])
            t = trans[pt].astype(int)
            if i == j and np.all(t == 0):
                continue
            saddle = max(float(E[pt]), float(other_E[pt]))
            if saddle - global_min > float(max_relative_eV):
                continue
            sg = p.astype(int) if float(E[pt]) >= float(other_E[pt]) else p.astype(int) + off
            if i > j:
                i, j = j, i
                t = -t
            elif i == j:
                nz = np.flatnonzero(t)
                if len(nz) and t[nz[0]] < 0:
                    t = -t
            key = (i, j, int(t[0]), int(t[1]), int(t[2]))
            bi = max(0.0, saddle - float(min_energy[i]))
            bj = max(0.0, saddle - float(min_energy[j]))
            row = {
                "edge_id": -1,
                "minimum_i": i,
                "minimum_j": j,
                "ta": int(t[0]),
                "tb": int(t[1]),
                "tc": int(t[2]),
                "saddle_energy_eV": saddle,
                "barrier_i_to_j_eV": bi,
                "barrier_j_to_i_eV": bj,
                "barrier_eV": max(bi, bj),
                "saddle_grid_i": int(sg[0]),
                "saddle_grid_j": int(sg[1]),
                "saddle_grid_k": int(sg[2]),
            }
            if key not in best or row["saddle_energy_eV"] < best[key]["saddle_energy_eV"]:
                best[key] = row
    rows = list(best.values())
    rows.sort(key=lambda x: (x["barrier_eV"], x["minimum_i"], x["minimum_j"], x["ta"], x["tb"], x["tc"]))
    for eid, row in enumerate(rows):
        row["edge_id"] = eid
    return rows


def _edge_length(cube, min_grid, edge):
    fi = min_grid[edge["minimum_i"]] / cube.shape
    fj = min_grid[edge["minimum_j"]] / cube.shape + np.array([edge["ta"], edge["tb"], edge["tc"]], float)
    return float(np.linalg.norm((fj - fi) @ cube.cell_A))


def _periodic_path(cube, min_grid, edges, target, tile_radius):
    adj = {}
    for e in edges:
        w = float(e["barrier_eV"])
        length = _edge_length(cube, min_grid, e)
        t = np.array([e["ta"], e["tb"], e["tc"]], int)
        for a in range(-tile_radius, tile_radius + 1):
            for b in range(-tile_radius, tile_radius + 1):
                for c in range(-tile_radius, tile_radius + 1):
                    u = (e["minimum_i"], a, b, c)
                    v = (e["minimum_j"], a + t[0], b + t[1], c + t[2])
                    if max(abs(v[1]), abs(v[2]), abs(v[3])) <= tile_radius:
                        adj.setdefault(u, []).append((v, w, length, e["edge_id"]))
                        adj.setdefault(v, []).append((u, w, length, e["edge_id"]))
    best = None
    target = np.asarray(target, int)
    for s in range(len(min_grid)):
        src = (s, 0, 0, 0)
        dst = (s, int(target[0]), int(target[1]), int(target[2]))
        dist = {src: (0.0, 0.0, 0)}
        prev = {}
        heap = [(0.0, 0.0, 0, src)]
        while heap:
            mb, length, hops, u = heapq.heappop(heap)
            cost = (mb, length, hops)
            if cost != dist.get(u):
                continue
            if u == dst:
                break
            for v, w, elen, eid in adj.get(u, []):
                nc = (max(mb, w), length + elen, hops + 1)
                if nc < dist.get(v, (np.inf, np.inf, 10**9)):
                    dist[v] = nc
                    prev[v] = (u, eid)
                    heapq.heappush(heap, (*nc, v))
        if dst not in dist:
            continue
        if best is None or dist[dst] < tuple(best["cost"]):
            nodes = [dst]
            edge_ids = []
            u = dst
            while u != src:
                pu, eid = prev[u]
                nodes.append(pu)
                edge_ids.append(eid)
                u = pu
            best = {
                "cost": [float(dist[dst][0]), float(dist[dst][1]), int(dist[dst][2])],
                "Ecrit_eV": float(dist[dst][0]),
                "path_length_A": float(dist[dst][1]),
                "n_hops": int(dist[dst][2]),
                "translation": [int(x) for x in target],
                "nodes": [[int(y) for y in x] for x in reversed(nodes)],
                "edge_ids": [int(x) for x in reversed(edge_ids)],
            }
    return best


class _PeriodicDSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n
        self.delta = [np.zeros(3, int) for _ in range(n)]

    def find(self, x):
        if self.parent[x] == x:
            return x, np.zeros(3, int)
        p = self.parent[x]
        root, d = self.find(p)
        total = self.delta[x] + d
        self.parent[x] = root
        self.delta[x] = total
        return root, total.copy()

    def add(self, i, j, t):
        ri, di = self.find(i)
        rj, dj = self.find(j)
        t = np.asarray(t, int)
        if ri == rj:
            return t - (dj - di)
        w = t + di - dj
        if self.rank[ri] < self.rank[rj]:
            self.parent[ri] = rj
            self.delta[ri] = -w
        else:
            self.parent[rj] = ri
            self.delta[rj] = w
            if self.rank[ri] == self.rank[rj]:
                self.rank[ri] += 1
        return np.zeros(3, int)


def _rank_thresholds(nmin, edges):
    dsu = _PeriodicDSU(nmin)
    cycles = []
    result = {1: None, 2: None, 3: None}
    for e in sorted(edges, key=lambda x: x["barrier_eV"]):
        v = dsu.add(e["minimum_i"], e["minimum_j"], [e["ta"], e["tb"], e["tc"]])
        if np.any(v):
            cycles.append(v.astype(float))
            rank = int(np.linalg.matrix_rank(np.vstack(cycles), tol=1e-10))
            for r in (1, 2, 3):
                if rank >= r and result[r] is None:
                    result[r] = float(e["barrier_eV"])
    return result


def _path_steps(cube, min_grid, edges, result):
    if result is None:
        return pd.DataFrame()
    by_id = {e["edge_id"]: e for e in edges}
    rows = []
    for order, eid in enumerate(result["edge_ids"]):
        e = by_id[eid]
        u = result["nodes"][order]
        v = result["nodes"][order + 1]
        sg = np.array([e["saddle_grid_i"], e["saddle_grid_j"], e["saddle_grid_k"]], float)
        sf = sg / cube.shape
        sc = sf @ cube.cell_A
        rows.append({
            "step": order,
            "edge_id": eid,
            "source_minimum": int(u[0]),
            "destination_minimum": int(v[0]),
            "source_tile_a": int(u[1]),
            "source_tile_b": int(u[2]),
            "source_tile_c": int(u[3]),
            "destination_tile_a": int(v[1]),
            "destination_tile_b": int(v[2]),
            "destination_tile_c": int(v[3]),
            "barrier_eV": float(e["barrier_eV"]),
            "saddle_energy_eV": float(e["saddle_energy_eV"]),
            "saddle_frac_a": float(sf[0]),
            "saddle_frac_b": float(sf[1]),
            "saddle_frac_c": float(sf[2]),
            "saddle_x_A": float(sc[0]),
            "saddle_y_A": float(sc[1]),
            "saddle_z_A": float(sc[2]),
        })
    return pd.DataFrame(rows)


def analyze(cube_path, outdir="pathway", max_relative_eV=3.0, tile_radius=2):
    cube = read_cube(cube_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    labels, shifts, min_grid, min_energy = _basins(cube.values)
    edges = _saddles(cube, labels, shifts, min_grid, min_energy, max_relative_eV)
    minima = pd.DataFrame({
        "minimum_id": np.arange(len(min_grid)),
        "grid_i": min_grid[:, 0],
        "grid_j": min_grid[:, 1],
        "grid_k": min_grid[:, 2],
        "frac_a": min_grid[:, 0] / cube.shape[0],
        "frac_b": min_grid[:, 1] / cube.shape[1],
        "frac_c": min_grid[:, 2] / cube.shape[2],
        "energy_eV": min_energy,
    })
    saddles = pd.DataFrame(edges)
    pa = _periodic_path(cube, min_grid, edges, [1, 0, 0], int(tile_radius))
    pb = _periodic_path(cube, min_grid, edges, [0, 1, 0], int(tile_radius))
    pc = _periodic_path(cube, min_grid, edges, [0, 0, 1], int(tile_radius))
    ranks = _rank_thresholds(len(min_grid), edges)
    minima.to_csv(outdir / "minima.csv", index=False)
    saddles.to_csv(outdir / "saddles.csv", index=False)
    _path_steps(cube, min_grid, edges, pa).to_csv(outdir / "path_a.csv", index=False)
    _path_steps(cube, min_grid, edges, pb).to_csv(outdir / "path_b.csv", index=False)
    _path_steps(cube, min_grid, edges, pc).to_csv(outdir / "path_c.csv", index=False)
    summary = {
        "cube": str(cube_path),
        "n_minima": int(len(min_grid)),
        "n_saddles": int(len(edges)),
        "Ecrit_a_eV": None if pa is None else pa["Ecrit_eV"],
        "Ecrit_b_eV": None if pb is None else pb["Ecrit_eV"],
        "Ecrit_c_eV": None if pc is None else pc["Ecrit_eV"],
        "Ecrit_1D_eV": ranks[1],
        "Ecrit_2D_eV": ranks[2],
        "Ecrit_3D_eV": ranks[3],
        "path_a": pa,
        "path_b": pb,
        "path_c": pc,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary
