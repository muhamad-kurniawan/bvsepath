from importlib.resources import files
from pathlib import Path
import json
import math


def parse_oxidation(value):
    p = Path(value)
    if p.exists():
        data = json.loads(p.read_text())
        return {str(k): int(v) for k, v in data.items()}
    out = {}
    for item in value.split(","):
        key, val = item.split("=", 1)
        out[key.strip()] = int(val.strip())
    return out


def _data_path(name, override=None):
    if override is not None:
        return Path(override)
    return files("bvsepath.data").joinpath(name)


def load_elements(path=None):
    out = {}
    source = _data_path("elements.dat", path)
    with source.open() as f:
        for line in f:
            if not line.strip():
                continue
            p = line.split()
            if len(p) < 10:
                continue
            symbol = p[1]
            valence = int(float(p[2]))
            out[(symbol, valence)] = {
                "Z": int(float(p[0])),
                "n": float(p[5]),
                "c": float(p[7]),
                "radius_A": float(p[8]),
                "softness": float(p[9]),
            }
    return out


def load_bvse(path=None):
    out = {}
    source = _data_path("bvse.dat", path)
    with source.open() as f:
        header = f.readline().split()
        if len(header) < 10:
            raise ValueError("bvse.dat must contain the 10-column BVSE parameter schema.")
        for line in f:
            if not line.strip():
                continue
            p = line.split()
            if len(p) < 10:
                continue
            cation = p[0]
            cv = int(float(p[1]))
            anion = p[2]
            av = int(float(p[3]))
            key = (cation, cv, anion, av)
            out.setdefault(key, []).append({
                "Nc": float(p[4]),
                "R0": float(p[5]),
                "Rcut_A": float(p[6]),
                "D0_eV": float(p[7]),
                "Rmin_A": float(p[8]),
                "alpha_Ainv": float(p[9]),
            })
    return out


def effective_charges(symbols, oxidation, elements):
    counts = {}
    for s in symbols:
        counts[s] = counts.get(s, 0) + 1
    raw = {}
    pos = 0.0
    neg = 0.0
    for s, count in counts.items():
        v = oxidation[s]
        entry = elements[(s, v)]
        q = v / math.sqrt(entry["n"])
        raw[s] = q
        if q > 0:
            pos += count * q
        elif q < 0:
            neg += count * q
    if pos <= 0 or neg >= 0:
        raise ValueError("Both positive and negative framework charge are required.")
    scale = math.sqrt(-neg / pos)
    return {s: q * scale if q > 0 else q / scale for s, q in raw.items()}


def pair_parameters(mobile, mobile_valence, other, other_valence, table):
    if mobile_valence > 0 and other_valence < 0:
        key = (mobile, mobile_valence, other, other_valence)
    elif mobile_valence < 0 and other_valence > 0:
        key = (other, other_valence, mobile, mobile_valence)
    else:
        return None
    rows = table.get(key)
    if not rows:
        raise KeyError(f"No BVSE parameters for {key}.")
    return rows[0]
