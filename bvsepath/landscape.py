from collections import Counter
from pathlib import Path
import json
import math
import numpy as np
from scipy.special import erfc
from scipy.spatial import cKDTree
from scipy.optimize import brentq, minimize_scalar
from .cube import Cube, write_cube
from .parameters import effective_charges, load_bvse, load_elements, pair_parameters, parse_oxidation
from .structure import read_structure

COULOMB_EV_A = 14.3996454784255
COULOMB_CUTOFF_A = 10.0


def _periodic_images(cell, positions, cutoff):
    inv = np.linalg.inv(cell)
    heights = [1.0 / np.linalg.norm(inv[:, i]) for i in range(3)]
    nr = [int(math.ceil(cutoff / h)) + 1 for h in heights]
    trans = np.array([(i, j, k) for i in range(-nr[0], nr[0] + 1) for j in range(-nr[1], nr[1] + 1) for k in range(-nr[2], nr[2] + 1)], int)
    shifts = trans @ cell
    return np.vstack([positions + s for s in shifts])


def _radial_sum(points, tree, cutoff, fn):
    neighborhoods = tree.query_ball_point(points, cutoff, workers=-1)
    counts = np.fromiter((len(x) for x in neighborhoods), dtype=np.int64, count=len(neighborhoods))
    if counts.sum() == 0:
        return np.zeros(len(points))
    pi = np.repeat(np.arange(len(points)), counts)
    ai = np.concatenate(neighborhoods)
    d = np.linalg.norm(points[pi] - tree.data[ai], axis=1)
    out = np.zeros(len(points))
    np.add.at(out, pi, fn(d))
    return out



def _scaled_structure(st, scale):
    return {
        "cell": np.asarray(st["cell"], float) * float(scale),
        "symbols": list(st["symbols"]),
        "numbers": np.asarray(st["numbers"], int).copy(),
        "frac": np.asarray(st["frac"], float).copy(),
    }


def _occupied_mobile_environment_energy(st, mobile, ox, elements, table, sf):
    symbols = list(st["symbols"])
    cell = np.asarray(st["cell"], float)
    positions = np.asarray(st["frac"], float) @ cell
    mobile_valence = int(ox[mobile])
    mobile_mask = np.array(symbols) == mobile
    points = positions[mobile_mask]
    if len(points) == 0:
        raise ValueError(f"sf auto requires at least one occupied {mobile} site in the input structure.")
    q = effective_charges(symbols, ox, elements)
    grouped = {}
    for s in sorted(set(symbols)):
        grouped[s] = positions[np.array(symbols) == s]
    trees = {s: cKDTree(_periodic_images(cell, grouped[s], COULOMB_CUTOFF_A)) for s in grouped}
    e = np.zeros(len(points), float)
    for s in grouped:
        if int(ox[s]) * mobile_valence >= 0:
            continue
        p = pair_parameters(mobile, mobile_valence, s, int(ox[s]), table)
        alpha = p["alpha_Ainv"]
        rmin = p["Rmin_A"]
        d0 = p["D0_eV"]
        cutoff = p["Rcut_A"]
        cutoff_value = (math.exp(alpha * (rmin - cutoff)) - 1.0) ** 2 - 1.0
        def morse(d, alpha=alpha, rmin=rmin, d0=d0, cutoff_value=cutoff_value):
            x = np.exp(alpha * (rmin - d))
            return 0.5 * d0 * (((x - 1.0) ** 2 - 1.0) - cutoff_value)
        e += _radial_sum(points, trees[s], cutoff, morse)
    mobile_radius = elements[(mobile, mobile_valence)]["radius_A"]
    for s in grouped:
        if int(ox[s]) * mobile_valence <= 0:
            continue
        v = int(ox[s])
        radius = elements[(s, v)]["radius_A"]
        rho = float(sf) * (mobile_radius + radius)
        qprod = q[mobile] * q[s]
        shift = erfc(COULOMB_CUTOFF_A / rho) / COULOMB_CUTOFF_A
        if s != mobile:
            def screened(d, rho=rho, qprod=qprod, shift=shift):
                return COULOMB_EV_A * qprod * (erfc(d / rho) / d - shift)
            e += _radial_sum(points, trees[s], COULOMB_CUTOFF_A, screened)
        else:
            neighborhoods = trees[s].query_ball_point(points, COULOMB_CUTOFF_A, workers=-1)
            for i, inds in enumerate(neighborhoods):
                if not inds:
                    continue
                d = np.linalg.norm(points[i] - trees[s].data[np.asarray(inds)], axis=1)
                d = d[d > 1e-8]
                if len(d):
                    e[i] += np.sum(COULOMB_EV_A * qprod * (erfc(d / rho) / d - shift))
    return float(np.sum(e))


def pressure_proxy(st, mobile, ox, elements, table, sf, strain=1.0e-3):
    sm = _scaled_structure(st, 1.0 - float(strain))
    sp = _scaled_structure(st, 1.0 + float(strain))
    em = _occupied_mobile_environment_energy(sm, mobile, ox, elements, table, sf)
    ep = _occupied_mobile_environment_energy(sp, mobile, ox, elements, table, sf)
    vm = abs(np.linalg.det(sm["cell"]))
    vp = abs(np.linalg.det(sp["cell"]))
    return -(ep - em) / (vp - vm)


def optimize_sf_auto(st, mobile, ox, elements, table, sf_min=0.50, sf_max=0.85, strain=1.0e-3):
    f = lambda x: pressure_proxy(st, mobile, ox, elements, table, x, strain=strain)
    p0 = float(f(sf_min))
    p1 = float(f(sf_max))
    if np.sign(p0) != np.sign(p1):
        sf = float(brentq(f, float(sf_min), float(sf_max), xtol=1e-5, rtol=1e-8))
        method = "pressure_proxy_root"
    else:
        res = minimize_scalar(lambda x: abs(f(x)), bounds=(float(sf_min), float(sf_max)), method="bounded", options={"xatol": 1e-5})
        sf = float(res.x)
        method = "pressure_proxy_minabs"
    meta = {
        "method": method,
        "sf_min": float(sf_min),
        "sf_max": float(sf_max),
        "strain": float(strain),
        "pressure_at_sf_eV_A3": float(f(sf)),
        "pressure_at_min_eV_A3": p0,
        "pressure_at_max_eV_A3": p1,
        "note": "Independent finite-difference mobile-environment pressure-balance proxy inspired by the softBV zero-pressure screening-factor concept; not the internal softBV sf_auto implementation.",
    }
    return sf, meta

def generate(structure, mobile, oxidation, output="bvse.cube", resolution_A=0.18, sf=0.74, elements_path=None, bvse_path=None, metadata=None, chunk_size=5000, sf_min=0.50, sf_max=0.85, sf_strain=1.0e-3):
    st = read_structure(structure)
    symbols = list(st["symbols"])
    ox = parse_oxidation(oxidation) if isinstance(oxidation, str) else {str(k): int(v) for k, v in oxidation.items()}
    missing = sorted(set(symbols) - set(ox))
    if missing:
        raise ValueError("Missing oxidation states for: " + ", ".join(missing))
    if mobile not in ox:
        raise ValueError(f"Missing oxidation state for mobile ion {mobile}.")
    mobile_valence = int(ox[mobile])
    elements = load_elements(elements_path)
    table = load_bvse(bvse_path)
    for s in set(symbols):
        if (s, int(ox[s])) not in elements:
            raise KeyError(f"No elements.dat entry for {(s, int(ox[s]))}.")
    q = effective_charges(symbols, ox, elements)
    sf_meta = {}
    if isinstance(sf, str) and sf.strip().lower() == "auto":
        sf, sf_meta = optimize_sf_auto(st, mobile, ox, elements, table, sf_min=sf_min, sf_max=sf_max, strain=sf_strain)
        sf_mode = "auto_pressure_proxy"
    else:
        sf = float(sf)
        sf_mode = "fixed"
    cell = np.asarray(st["cell"], float)
    lengths = np.linalg.norm(cell, axis=1)
    shape = np.maximum(2, np.ceil(lengths / float(resolution_A)).astype(int))
    grid = cell / shape[:, None]
    positions = np.asarray(st["frac"], float) @ cell
    numbers = np.asarray(st["numbers"], int)
    atom_rows = np.array([[int(z), q[s], *r] for z, s, r in zip(numbers, symbols, positions)], float)
    grouped = {}
    for s in sorted(set(symbols)):
        grouped[s] = positions[np.array(symbols) == s]
    trees = {s: cKDTree(_periodic_images(cell, grouped[s], COULOMB_CUTOFF_A)) for s in grouped}
    opposite = [s for s in grouped if ox[s] * mobile_valence < 0]
    same = [s for s in grouped if ox[s] * mobile_valence > 0 and s != mobile]
    pair = {s: pair_parameters(mobile, mobile_valence, s, int(ox[s]), table) for s in opposite}
    mobile_radius = elements[(mobile, mobile_valence)]["radius_A"]
    values = np.empty(int(np.prod(shape)), float)
    for start in range(0, values.size, int(chunk_size)):
        stop = min(start + int(chunk_size), values.size)
        flat = np.arange(start, stop)
        ijk = np.array(np.unravel_index(flat, tuple(shape))).T
        pts = ijk @ grid
        e = np.zeros(len(pts), float)
        for s in opposite:
            p = pair[s]
            alpha = p["alpha_Ainv"]
            rmin = p["Rmin_A"]
            d0 = p["D0_eV"]
            cutoff = p["Rcut_A"]
            cutoff_value = (math.exp(alpha * (rmin - cutoff)) - 1.0) ** 2 - 1.0
            def morse(d, alpha=alpha, rmin=rmin, d0=d0, cutoff_value=cutoff_value):
                x = np.exp(alpha * (rmin - d))
                return 0.5 * d0 * (((x - 1.0) ** 2 - 1.0) - cutoff_value)
            e += _radial_sum(pts, trees[s], cutoff, morse)
        for s in same:
            v = int(ox[s])
            radius = elements[(s, v)]["radius_A"]
            rho = float(sf) * (mobile_radius + radius)
            qprod = q[mobile] * q[s]
            shift = erfc(COULOMB_CUTOFF_A / rho) / COULOMB_CUTOFF_A
            def screened(d, rho=rho, qprod=qprod, shift=shift):
                return COULOMB_EV_A * qprod * (erfc(d / rho) / d - shift)
            e += _radial_sum(pts, trees[s], COULOMB_CUTOFF_A, screened)
        values[start:stop] = e
    values = values.reshape(tuple(shape), order="C")
    counts = Counter(symbols)
    desc = " ".join(f"{s}{abs(ox[s])}{'+' if ox[s] > 0 else '-'}:{counts[s]}" for s in sorted(counts))
    cube = Cube(desc, f"mobile={mobile}{mobile_valence:+d}; sf={float(sf):.6f}; sf_mode={sf_mode}; resolution_A={float(resolution_A):.6f}", np.zeros(3), shape, grid, cell, atom_rows, values)
    write_cube(cube, output)
    info = {
        "structure": str(structure),
        "mobile": mobile,
        "mobile_valence": mobile_valence,
        "oxidation_states": ox,
        "resolution_A": float(resolution_A),
        "sf": float(sf),
        "sf_mode": sf_mode,
        "sf_auto": sf_meta,
        "grid_shape": [int(x) for x in shape],
        "grid_step_A": [float(np.linalg.norm(v)) for v in grid],
        "output": str(output),
    }
    if metadata is not None:
        Path(metadata).write_text(json.dumps(info, indent=2))
    return info
