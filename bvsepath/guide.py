import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .cube import read_cube
from .structure import read_structure, write_poscar


def _cart_delta(frac1, frac2, cell):
    d = (np.asarray(frac1) - np.asarray(frac2) + 0.5) % 1.0 - 0.5
    return d @ cell


def _distance_matrix(frac_a, frac_b, cell):
    out = np.zeros((len(frac_a), len(frac_b)), float)
    for i, a in enumerate(frac_a):
        for j, b in enumerate(frac_b):
            out[i, j] = np.linalg.norm(_cart_delta(a, b, cell))
    return out


def _greedy_assignment(frac_mobile, frac_minima, cell):
    if len(frac_mobile) == 0 or len(frac_minima) == 0:
        return {}, {}
    dist = _distance_matrix(frac_mobile, frac_minima, cell)
    pairs = [(dist[i, j], i, j) for i in range(dist.shape[0]) for j in range(dist.shape[1])]
    pairs.sort(key=lambda x: x[0])
    used_i = set()
    used_j = set()
    min_to_atom = {}
    atom_to_min = {}
    for d, i, j in pairs:
        if i in used_i or j in used_j:
            continue
        used_i.add(i)
        used_j.add(j)
        min_to_atom[j] = i
        atom_to_min[i] = j
        if len(used_i) == len(frac_mobile):
            break
    return min_to_atom, atom_to_min


def _unique_path_nodes(path_df):
    if len(path_df) == 0:
        return []
    out = [int(path_df.iloc[0]["source_minimum"])]
    for _, row in path_df.iterrows():
        out.append(int(row["destination_minimum"]))
    return out


def _candidate_from_path(path_df, occ_min_to_atom):
    if len(path_df) == 0:
        return None
    rows = []
    for idx, row in path_df.iterrows():
        s = int(row["source_minimum"])
        d = int(row["destination_minimum"])
        os = s in occ_min_to_atom
        od = d in occ_min_to_atom
        rows.append((float(row["barrier_eV"]), idx, s, d, os, od))
    one_occ = [x for x in rows if x[4] != x[5]]
    if one_occ:
        one_occ.sort(key=lambda x: (-x[0], x[1]))
        _, idx, s, d, os, od = one_occ[0]
        if os:
            source_min, dest_min = s, d
            atom_local = occ_min_to_atom[s]
        else:
            source_min, dest_min = d, s
            atom_local = occ_min_to_atom[d]
        return {"path_row": int(idx), "source_minimum": int(source_min), "destination_minimum": int(dest_min), "mobile_local_index": int(atom_local)}
    rows.sort(key=lambda x: (-x[0], x[1]))
    _, idx, s, d, os, od = rows[0]
    if os:
        source_min, dest_min, atom_local = s, d, occ_min_to_atom[s]
    elif od:
        source_min, dest_min, atom_local = d, s, occ_min_to_atom[d]
    else:
        source_min, dest_min, atom_local = s, d, None
    return {"path_row": int(idx), "source_minimum": int(source_min), "destination_minimum": int(dest_min), "mobile_local_index": None if atom_local is None else int(atom_local)}


def _path_polyline(path_df, minima_df):
    if len(path_df) == 0:
        return np.zeros((0, 3), float), np.zeros((0, 3), float)
    pts = []
    labels = []
    for k, row in enumerate(path_df.itertuples(index=False)):
        sf = minima_df.loc[int(row.source_minimum), ["frac_a", "frac_b", "frac_c"]].to_numpy(float) + np.array([row.source_tile_a, row.source_tile_b, row.source_tile_c], float)
        df = minima_df.loc[int(row.destination_minimum), ["frac_a", "frac_b", "frac_c"]].to_numpy(float) + np.array([row.destination_tile_a, row.destination_tile_b, row.destination_tile_c], float)
        saddle = np.array([row.saddle_frac_a, row.saddle_frac_b, row.saddle_frac_c], float)
        mid = 0.5 * (sf + df)
        saddle = saddle + np.round(mid - saddle)
        if k == 0:
            pts.append(sf)
            labels.append(("min", int(row.source_minimum)))
        pts.append(saddle)
        labels.append(("saddle", int(row.edge_id)))
        pts.append(df)
        labels.append(("min", int(row.destination_minimum)))
    return np.asarray(pts, float), labels


def _project_energy(E, mode):
    if mode == "ab":
        return E.min(axis=2)
    if mode == "ac":
        return E.min(axis=1)
    return E.min(axis=0)


def _project_points(frac_pts, mode):
    if mode == "ab":
        return frac_pts[:, 0], frac_pts[:, 1]
    if mode == "ac":
        return frac_pts[:, 0], frac_pts[:, 2]
    return frac_pts[:, 1], frac_pts[:, 2]


def _plot_projection(ax, E2, frac_poly, source_frac, dest_frac, saddle_frac, mode, title):
    tiled = np.tile(E2.T, (3, 3))
    ax.imshow(tiled, origin="lower", extent=[-1, 2, -1, 2], aspect="equal")
    x, y = _project_points(frac_poly, mode)
    ax.plot(x, y, linewidth=2.0)
    sx, sy = _project_points(source_frac[None, :], mode)
    dx, dy = _project_points(dest_frac[None, :], mode)
    bx, by = _project_points(saddle_frac[None, :], mode)
    ax.scatter(sx, sy, s=50, marker="o")
    ax.scatter(dx, dy, s=50, marker="s")
    ax.scatter(bx, by, s=90, marker="*")
    ax.text(float(sx[0]), float(sy[0]), " source", fontsize=8)
    ax.text(float(dx[0]), float(dy[0]), " target", fontsize=8)
    ax.text(float(bx[0]), float(by[0]), " bottleneck", fontsize=8)
    ax.set_xlim(float(np.floor(min(x.min(), sx.min(), dx.min(), bx.min()) - 0.2)), float(np.ceil(max(x.max(), sx.max(), dx.max(), bx.max()) + 0.2)))
    ax.set_ylim(float(np.floor(min(y.min(), sy.min(), dy.min(), by.min()) - 0.2)), float(np.ceil(max(y.max(), sy.max(), dy.max(), by.max()) + 0.2)))
    ax.set_xlabel(mode[0])
    ax.set_ylabel(mode[1])
    ax.set_title(title)


def _write_marker_xyz(path, structure, source_frac, dest_frac, saddle_frac, source_min, dest_min, barrier):
    cell = np.asarray(structure["cell"], float)
    pts = [source_frac @ cell, dest_frac @ cell, saddle_frac @ cell]
    labels = [f"SRC_min{source_min}", f"DST_min{dest_min}", f"BOT_{barrier:.3f}eV"]
    with Path(path).open("w") as f:
        f.write("3\n")
        f.write("neb guide markers\n")
        for lab, p in zip(labels, pts):
            f.write(f"H {p[0]:.8f} {p[1]:.8f} {p[2]:.8f} {lab}\n")


def _write_candidate_files(outdir, tag, structure, mobile_indices, mobile_local_to_global, candidate, minima_df, path_df):
    row = path_df.iloc[int(candidate["path_row"])]
    source_min = int(candidate["source_minimum"])
    dest_min = int(candidate["destination_minimum"])
    source_frac = minima_df.loc[source_min, ["frac_a", "frac_b", "frac_c"]].to_numpy(float) % 1.0
    dest_frac = minima_df.loc[dest_min, ["frac_a", "frac_b", "frac_c"]].to_numpy(float) % 1.0
    saddle_frac = np.array([row["saddle_frac_a"], row["saddle_frac_b"], row["saddle_frac_c"]], float) % 1.0
    global_atom = None if candidate["mobile_local_index"] is None else int(mobile_local_to_global[int(candidate["mobile_local_index"])])
    initial = {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else list(v) if isinstance(v, list) else v) for k, v in structure.items()}
    initial["frac"] = np.array(structure["frac"], float).copy()
    final = {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else list(v) if isinstance(v, list) else v) for k, v in structure.items()}
    final["frac"] = np.array(structure["frac"], float).copy()
    if global_atom is not None:
        final["frac"][global_atom] = dest_frac
    labels = [None] * len(structure["symbols"])
    if global_atom is not None:
        labels[global_atom] = f"MOVE_atom{global_atom+1}_to_min{dest_min}"
    order = np.arange(len(structure["symbols"]))
    if global_atom is not None:
        order = np.array([global_atom] + [i for i in order if i != global_atom], int)
    info = {
        "direction": tag,
        "source_minimum": source_min,
        "destination_minimum": dest_min,
        "barrier_eV": float(row["barrier_eV"]),
        "edge_id": int(row["edge_id"]),
        "source_frac": source_frac.tolist(),
        "destination_frac": dest_frac.tolist(),
        "saddle_frac": saddle_frac.tolist(),
        "mobile_global_atom_index_0based": global_atom,
        "mobile_global_atom_index_1based": None if global_atom is None else int(global_atom + 1),
    }
    comment = f"NEB candidate {tag}; source_min={source_min}; dest_min={dest_min}; barrier={float(row['barrier_eV']):.6f} eV; move_atom_1based={'' if global_atom is None else global_atom+1}"
    write_poscar(outdir / f"POSCAR_{tag}_initial.vasp", initial, comment=comment, order=order)
    write_poscar(outdir / f"POSCAR_{tag}_final.vasp", final, comment=comment, order=order)
    write_poscar(outdir / f"POSCAR_{tag}_initial_tagged.vasp", initial, comment=comment, labels=labels, order=order)
    write_poscar(outdir / f"POSCAR_{tag}_final_tagged.vasp", final, comment=comment, labels=labels, order=order)
    pd.DataFrame([info]).to_csv(outdir / f"candidate_{tag}.csv", index=False)
    (outdir / f"candidate_{tag}.json").write_text(json.dumps(info, indent=2))
    _write_marker_xyz(outdir / f"candidate_{tag}_markers.xyz", structure, source_frac, dest_frac, saddle_frac, source_min, dest_min, float(row["barrier_eV"]))
    return info


def make_guides(structure_path, pathway_dir, mobile, cube_path=None, directions=None, outdir="neb_guides"):
    pathway_dir = Path(pathway_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads((pathway_dir / "summary.json").read_text())
    cube_path = cube_path if cube_path is not None else summary["cube"]
    cube = read_cube(cube_path)
    structure = read_structure(structure_path)
    minima = pd.read_csv(pathway_dir / "minima.csv").set_index("minimum_id")
    mobile_global = [i for i, s in enumerate(structure["symbols"]) if s == mobile]
    if len(mobile_global) == 0:
        raise ValueError(f"No mobile-ion atoms with symbol {mobile} were found in the structure.")
    frac_mobile = np.asarray(structure["frac"], float)[mobile_global]
    frac_min = minima[["frac_a", "frac_b", "frac_c"]].to_numpy(float)
    min_to_local, atom_to_min = _greedy_assignment(frac_mobile, frac_min, np.asarray(structure["cell"], float))
    occ_min_to_atom = {int(k): int(v) for k, v in min_to_local.items()}
    if directions is None:
        directions = ["a", "b", "c"]
    outputs = {"structure": str(structure_path), "cube": str(cube_path), "mobile": mobile, "directions": {}}
    for tag in directions:
        path_file = pathway_dir / f"path_{tag}.csv"
        if not path_file.exists():
            continue
        path_df = pd.read_csv(path_file)
        if len(path_df) == 0:
            continue
        cand = _candidate_from_path(path_df, occ_min_to_atom)
        if cand is None:
            continue
        info = _write_candidate_files(outdir, tag, structure, mobile_global, mobile_global, cand, minima, path_df)
        poly, _ = _path_polyline(path_df, minima)
        row = path_df.iloc[int(cand["path_row"])]
        source_frac = minima.loc[int(info["source_minimum"]), ["frac_a", "frac_b", "frac_c"]].to_numpy(float) + 0.0
        dest_frac = minima.loc[int(info["destination_minimum"]), ["frac_a", "frac_b", "frac_c"]].to_numpy(float) + 0.0
        saddle_frac = np.array([row["saddle_frac_a"], row["saddle_frac_b"], row["saddle_frac_c"]], float)
        mid = 0.5 * (source_frac + dest_frac)
        saddle_frac = saddle_frac + np.round(mid - saddle_frac)
        fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
        for ax, mode in zip(axs, ["ab", "ac", "bc"]):
            _plot_projection(ax, _project_energy(cube.values - float(cube.values.min()), mode), poly, source_frac, dest_frac, saddle_frac, mode, f"{tag}-direction {mode}")
        fig.suptitle(f"NEB guide {tag}: edge {int(row['edge_id'])}, barrier {float(row['barrier_eV']):.3f} eV")
        fig.tight_layout()
        fig.savefig(outdir / f"candidate_{tag}_guide.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        outputs["directions"][tag] = info
    (outdir / "neb_guides_summary.json").write_text(json.dumps(outputs, indent=2))
    return outputs
