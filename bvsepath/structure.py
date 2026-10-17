from pathlib import Path
import math
import re
import numpy as np

ELEMENTS = [
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
    "Md","No","Lr"
]
Z = {s: i + 1 for i, s in enumerate(ELEMENTS)}


def _cell_from_lengths_angles(a, b, c, alpha, beta, gamma):
    ar = math.radians(alpha)
    br = math.radians(beta)
    gr = math.radians(gamma)
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(gr), b * math.sin(gr), 0.0])
    cx = c * math.cos(br)
    cy = c * (math.cos(ar) - math.cos(br) * math.cos(gr)) / math.sin(gr)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([va, vb, [cx, cy, cz]], float)


def read_poscar(path):
    lines = [x.rstrip() for x in Path(path).read_text().splitlines() if x.strip()]
    scale = float(lines[1].split()[0])
    raw_cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)], float)
    if scale < 0:
        target_volume = abs(scale)
        current = abs(np.linalg.det(raw_cell))
        scale = (target_volume / current) ** (1.0 / 3.0)
    cell = raw_cell * scale
    symbols = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].lower().startswith("s"):
        idx += 1
    direct = lines[idx].lower().startswith("d")
    idx += 1
    species = []
    for s, n in zip(symbols, counts):
        species.extend([s] * n)
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]] for i in range(sum(counts))], float)
    if direct:
        frac = coords
    else:
        cart = coords * scale
        frac = cart @ np.linalg.inv(cell)
    frac %= 1.0
    numbers = np.array([Z[s] for s in species], int)
    return {"cell": cell, "symbols": species, "numbers": numbers, "frac": frac}


def _clean_cif_token(x):
    x = x.strip().strip("'").strip('"')
    x = re.sub(r"\([^)]*\)$", "", x)
    return x


def read_cif(path):
    lines = Path(path).read_text().splitlines()
    values = {}
    for line in lines:
        q = line.strip()
        if q.startswith("_cell_"):
            p = q.split(None, 1)
            if len(p) == 2:
                values[p[0]] = float(_clean_cif_token(p[1]))
    keys = ["_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]
    if any(k not in values for k in keys):
        raise ValueError("CIF cell parameters are incomplete.")
    cell = _cell_from_lengths_angles(*(values[k] for k in keys))
    headers = []
    rows = []
    i = 0
    while i < len(lines):
        if lines[i].strip().lower() == "loop_":
            j = i + 1
            local = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                local.append(lines[j].strip().split()[0])
                j += 1
            needed = {"_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"}
            if needed.issubset(local):
                while j < len(lines):
                    q = lines[j].strip()
                    if not q or q.startswith("#"):
                        j += 1
                        continue
                    if q.startswith("_") or q.lower() == "loop_" or q.lower().startswith("data_"):
                        break
                    p = q.split()
                    if len(p) >= len(local):
                        rows.append(p[:len(local)])
                    j += 1
                headers = local
                break
            i = j
        else:
            i += 1
    if not rows:
        raise ValueError("No fractional atom-site loop was found in the CIF.")
    ix = headers.index("_atom_site_fract_x")
    iy = headers.index("_atom_site_fract_y")
    iz = headers.index("_atom_site_fract_z")
    if "_atom_site_type_symbol" in headers:
        isym = headers.index("_atom_site_type_symbol")
    elif "_atom_site_label" in headers:
        isym = headers.index("_atom_site_label")
    else:
        raise ValueError("No atom symbol column was found in the CIF.")
    species = []
    frac = []
    for row in rows:
        raw = _clean_cif_token(row[isym])
        m = re.match(r"([A-Z][a-z]?)", raw)
        if not m:
            raise ValueError(f"Cannot parse element symbol from {raw}.")
        s = m.group(1)
        species.append(s)
        frac.append([float(_clean_cif_token(row[ix])), float(_clean_cif_token(row[iy])), float(_clean_cif_token(row[iz]))])
    frac = np.asarray(frac, float) % 1.0
    numbers = np.array([Z[s] for s in species], int)
    return {"cell": cell, "symbols": species, "numbers": numbers, "frac": frac}


def read_structure(path):
    path = Path(path)
    if path.suffix.lower() == ".cif":
        return read_cif(path)
    return read_poscar(path)


def write_poscar(path, structure, comment="Generated by bvsepath", labels=None, order=None):
    cell = np.asarray(structure["cell"], float)
    frac = np.asarray(structure["frac"], float) % 1.0
    symbols = list(structure["symbols"])
    if order is None:
        order = np.arange(len(symbols))
    else:
        order = np.asarray(order, int)
    frac = frac[order]
    symbols = [symbols[i] for i in order]
    labels = None if labels is None else [labels[i] for i in order]
    species = []
    for s in symbols:
        if s not in species:
            species.append(s)
    grouped = []
    grouped_labels = [] if labels is not None else None
    counts = []
    for s in species:
        idx = [i for i, x in enumerate(symbols) if x == s]
        counts.append(len(idx))
        grouped.extend(frac[idx])
        if labels is not None:
            grouped_labels.extend([labels[i] for i in idx])
    with Path(path).open("w") as f:
        f.write(str(comment).rstrip() + "\n")
        f.write("1.0\n")
        for v in cell:
            f.write(f"{v[0]:16.10f} {v[1]:16.10f} {v[2]:16.10f}\n")
        f.write(" ".join(species) + "\n")
        f.write(" ".join(str(x) for x in counts) + "\n")
        if labels is not None:
            f.write("Selective dynamics\n")
        f.write("Direct\n")
        for i, q in enumerate(grouped):
            line = f"{q[0]:16.10f} {q[1]:16.10f} {q[2]:16.10f}"
            if labels is not None:
                lab = grouped_labels[i] if grouped_labels[i] is not None else f"A{i+1}"
                line += f" T T T {lab}"
            f.write(line + "\n")
