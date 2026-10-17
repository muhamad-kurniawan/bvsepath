from importlib.resources import files
from pathlib import Path
import json
import math
import pandas as pd

COULOMB_EV_A = 14.3996454784255

ELEMENTS = [
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
    "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
]
Z = {s: i + 1 for i, s in enumerate(ELEMENTS)}
D_BLOCK = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Rf Db Sg Bh Hs Mt Ds Rg Cn".split())
F_BLOCK = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split())


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


def principal_n(element):
    z = Z[element]
    if z <= 2:
        return 1.0
    if z <= 10:
        return 2.0
    if z <= 18:
        return 3.0
    if z <= 36:
        return 4.0
    if z <= 54:
        return 5.0
    if z <= 86:
        return 6.0
    return 7.0


def block_factor(element):
    return 2.0 if element in D_BLOCK or element in F_BLOCK else 1.0


def load_ions(path=None):
    df = pd.read_csv(_data_path("softbv_ion_properties.csv", path))
    out = {}
    for row in df.itertuples(index=False):
        out[(str(row.element), int(row.oxidation))] = {
            "n": principal_n(str(row.element)),
            "c": block_factor(str(row.element)),
            "radius_A": float(row.screening_radius_A),
            "softness": float(row.softness_eV_inv),
        }
    return out


def load_pairs(path=None):
    df = pd.read_csv(_data_path("softbv_pairs_2017.csv", path))
    out = {}
    for row in df.itertuples(index=False):
        key = (str(row.cation), int(row.cation_valence), str(row.anion), int(row.anion_valence))
        out[key] = {
            "n_reference": int(row.n_reference),
            "b_A": float(row.b_A),
            "R0_A": float(row.R0_A),
            "Rcut_A": float(row.Rcut_A),
            "deltaV_softBV": float(row.deltaV_softBV),
            "Nc": float(row.Nc),
            "R0_softNC1_A": float(row.R0_softNC1_A),
            "Rcut_softNC1_A": float(row.Rcut_softNC1_A),
            "deltaV_softNC1": float(row.deltaV_softNC1),
            "source": str(row.source),
        }
    return out


def effective_charges(symbols, oxidation, ions):
    counts = {}
    for s in symbols:
        counts[s] = counts.get(s, 0) + 1
    raw = {}
    pos = 0.0
    neg = 0.0
    for s, count in counts.items():
        v = int(oxidation[s])
        if (s, v) not in ions:
            raise KeyError(f"No SoftBV ion-property entry for {(s, v)}.")
        q = v / math.sqrt(principal_n(s))
        raw[s] = q
        if q > 0:
            pos += count * q
        elif q < 0:
            neg += count * q
    if pos <= 0 or neg >= 0:
        raise ValueError("Both positive and negative framework charge are required.")
    scale = math.sqrt(-neg / pos)
    return {s: q * scale if q > 0 else q / scale for s, q in raw.items()}


def pair_key(mobile, mobile_valence, other, other_valence):
    if mobile_valence > 0 and other_valence < 0:
        return (mobile, int(mobile_valence), other, int(other_valence))
    if mobile_valence < 0 and other_valence > 0:
        return (other, int(other_valence), mobile, int(mobile_valence))
    return None


def derived_pair_parameters(mobile, mobile_valence, other, other_valence, pairs, ions):
    key = pair_key(mobile, mobile_valence, other, other_valence)
    if key is None or key not in pairs:
        raise KeyError(f"No published 2017 SoftBV pair entry for {key} in the bundled table.")
    p = pairs[key]
    cation, cv, anion, av = key
    if (cation, cv) not in ions or (anion, av) not in ions:
        raise KeyError(f"Missing ion properties required for {key}.")
    ic = ions[(cation, cv)]
    ia = ions[(anion, av)]
    b = p["b_A"]
    softness_factor = 0.9185 + 0.2285 * abs(ic["softness"] - ia["softness"])
    rmin = softness_factor * p["R0_A"] - b * math.log(abs(cv) / p["Nc"])
    c = block_factor(cation)
    d0 = (b * b / 2.0) * COULOMB_EV_A * c * (abs(cv * av) ** (1.0 / c)) / (rmin * math.sqrt(principal_n(cation) * principal_n(anion)))
    return {
        "Nc": p["Nc"],
        "R0_A": p["R0_A"],
        "Rcut_A": p["Rcut_A"],
        "b_A": b,
        "alpha_Ainv": 1.0 / b,
        "Rmin_A": rmin,
        "D0_eV": d0,
        "softness_factor": softness_factor,
        "pair_source": p["source"],
        "Nc_source": "softNC1_2017",
    }
