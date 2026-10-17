from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .cube import read_cube
from .landscape import generate
from .pathway import analyze
from .guide import make_guides, _path_polyline, _project_energy, _project_points


_PATH_FILES = {
    "a": "path_a.csv",
    "b": "path_b.csv",
    "c": "path_c.csv",
    "periodic_a": "path_periodic_a.csv",
    "periodic_b": "path_periodic_b.csv",
    "periodic_c": "path_periodic_c.csv",
    "periodic_1d": "path_periodic_1d.csv",
    "1d": "path_periodic_1d.csv",
}


class BVSEWorkflow:
    def __init__(self, structure=None, mobile=None, oxidation=None, workdir="bvse_run", cube=None):
        self.structure = None if structure is None else str(structure)
        self.mobile = mobile
        self.oxidation = oxidation
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.cube_path = Path(cube) if cube is not None else self.workdir / "bvse.cube"
        self.landscape_metadata_path = self.workdir / "landscape.json"
        self.pathway_dir = self.workdir / "pathway"
        self.neb_dir = self.workdir / "neb_guides"
        self.landscape_metadata = None
        self.summary = None
        self.minima = None
        self.saddles = None

    @classmethod
    def from_cube(cls, cube, workdir="bvse_run", structure=None, mobile=None):
        return cls(structure=structure, mobile=mobile, oxidation=None, workdir=workdir, cube=cube)

    def landscape(self, resolution=0.18, sf=0.74, ions=None, pairs=None, chunk_size=5000, sf_min=0.50, sf_max=0.85, sf_strain=1.0e-3, grid_shape=None):
        if self.structure is None:
            raise ValueError("A structure path is required to generate a BVSE landscape.")
        if self.mobile is None:
            raise ValueError("mobile must be specified.")
        if self.oxidation is None:
            raise ValueError("oxidation must be specified.")
        meta = generate(
            self.structure,
            self.mobile,
            self.oxidation,
            output=str(self.cube_path),
            resolution_A=resolution,
            sf=sf,
            ions_path=ions,
            pairs_path=pairs,
            metadata=str(self.landscape_metadata_path),
            chunk_size=chunk_size,
            sf_min=sf_min,
            sf_max=sf_max,
            sf_strain=sf_strain,
            grid_shape=grid_shape,
        )
        self.landscape_metadata = meta
        return meta

    def pathways(self, max_relative=3.0, tile_radius=3, max_repeat=3):
        if not self.cube_path.exists():
            raise FileNotFoundError(f"Cube file not found: {self.cube_path}")
        self.summary = analyze(
            str(self.cube_path),
            outdir=str(self.pathway_dir),
            max_relative_eV=max_relative,
            tile_radius=tile_radius,
            max_repeat=max_repeat,
        )
        self._load_tables()
        return self.summary

    def run(self, resolution=0.18, sf=0.74, max_relative=3.0, tile_radius=3, max_repeat=3, ions=None, pairs=None, chunk_size=5000, sf_min=0.50, sf_max=0.85, sf_strain=1.0e-3, grid_shape=None):
        self.landscape(
            resolution=resolution,
            sf=sf,
            ions=ions,
            pairs=pairs,
            chunk_size=chunk_size,
            sf_min=sf_min,
            sf_max=sf_max,
            sf_strain=sf_strain,
            grid_shape=grid_shape,
        )
        self.pathways(max_relative=max_relative, tile_radius=tile_radius, max_repeat=max_repeat)
        return self

    def load(self):
        if self.landscape_metadata_path.exists():
            self.landscape_metadata = json.loads(self.landscape_metadata_path.read_text())
        if (self.pathway_dir / "summary.json").exists():
            self.summary = json.loads((self.pathway_dir / "summary.json").read_text())
            self._load_tables()
        return self

    def _load_tables(self):
        minima = self.pathway_dir / "minima.csv"
        saddles = self.pathway_dir / "saddles.csv"
        self.minima = pd.read_csv(minima) if minima.exists() else None
        self.saddles = pd.read_csv(saddles) if saddles.exists() else None

    def path(self, name="periodic_1d"):
        key = str(name).lower().strip()
        if key not in _PATH_FILES:
            raise ValueError("Path name must be one of: " + ", ".join(sorted(_PATH_FILES)))
        p = self.pathway_dir / _PATH_FILES[key]
        if not p.exists():
            raise FileNotFoundError(f"Path result not found: {p}. Run pathways() first.")
        return pd.read_csv(p)

    def metrics(self):
        if self.summary is None:
            self.load()
        if self.summary is None:
            raise RuntimeError("No pathway summary is available. Run pathways() first.")
        keys = [
            "Ecrit_periodic_1D_path_eV",
            "periodic_1D_translation",
            "Ecrit_1D_eV",
            "Ecrit_2D_eV",
            "Ecrit_3D_eV",
            "Ecrit_periodic_a_eV",
            "Ecrit_periodic_b_eV",
            "Ecrit_periodic_c_eV",
            "periodic_a_translation",
            "periodic_b_translation",
            "periodic_c_translation",
        ]
        row = {k: self.summary.get(k) for k in keys}
        e1 = self.summary.get("Ecrit_1D_eV")
        e2 = self.summary.get("Ecrit_2D_eV")
        row["delta_E_2D_minus_1D_eV"] = None if e1 is None or e2 is None else float(e2) - float(e1)
        return pd.DataFrame([row])

    def periodic_candidates(self):
        p = self.pathway_dir / "periodic_axis_candidates.csv"
        if not p.exists():
            raise FileNotFoundError(f"Periodic candidate table not found: {p}")
        return pd.read_csv(p)

    def neb_guides(self, directions=None, outdir=None):
        if self.structure is None:
            raise ValueError("A structure path is required to prepare NEB guides.")
        if self.mobile is None:
            raise ValueError("mobile must be specified to prepare NEB guides.")
        if not (self.pathway_dir / "summary.json").exists():
            raise FileNotFoundError("Pathway results were not found. Run pathways() first.")
        target = self.neb_dir if outdir is None else Path(outdir)
        result = make_guides(
            self.structure,
            self.pathway_dir,
            self.mobile,
            cube_path=str(self.cube_path),
            directions=directions,
            outdir=str(target),
        )
        return result

    def plot_path(self, name="periodic_1d", projection="ab", energy_max=None, figsize=(6.5, 5.5), save=None):
        projection = str(projection).lower().strip()
        if projection not in {"ab", "ac", "bc"}:
            raise ValueError("projection must be 'ab', 'ac', or 'bc'.")
        path_df = self.path(name)
        if len(path_df) == 0:
            raise ValueError(f"Path {name} is empty.")
        if self.minima is None:
            self._load_tables()
        minima = self.minima.set_index("minimum_id")
        cube = read_cube(self.cube_path)
        e2 = _project_energy(cube.values - float(cube.values.min()), projection)
        poly, _ = _path_polyline(path_df, minima)
        x, y = _project_points(poly, projection)
        fig, ax = plt.subplots(figsize=figsize)
        tiled = np.tile(e2.T, (3, 3))
        im = ax.imshow(tiled, origin="lower", extent=[-1, 2, -1, 2], aspect="equal", vmin=0.0, vmax=energy_max)
        ax.plot(x, y, linewidth=2.0)
        ax.scatter(x[0], y[0], marker="o", s=55, label="start")
        ax.scatter(x[-1], y[-1], marker="s", s=55, label="end")
        idx = int(np.argmax(path_df["barrier_eV"].to_numpy(float)))
        row = path_df.iloc[idx]
        source_frac = minima.loc[int(row["source_minimum"]), ["frac_a", "frac_b", "frac_c"]].to_numpy(float) + np.array([row["source_tile_a"], row["source_tile_b"], row["source_tile_c"]], float)
        destination_frac = minima.loc[int(row["destination_minimum"]), ["frac_a", "frac_b", "frac_c"]].to_numpy(float) + np.array([row["destination_tile_a"], row["destination_tile_b"], row["destination_tile_c"]], float)
        bottleneck = np.array([row["saddle_frac_a"], row["saddle_frac_b"], row["saddle_frac_c"]], float)
        bottleneck = bottleneck + np.round(0.5 * (source_frac + destination_frac) - bottleneck)
        bx, by = _project_points(bottleneck[None, :], projection)
        ax.scatter(bx, by, marker="*", s=140, label="bottleneck")
        ax.set_xlabel(projection[0])
        ax.set_ylabel(projection[1])
        ax.set_title(f"{name}: Ecrit = {path_df['barrier_eV'].max():.3f} eV")
        ax.legend()
        fig.colorbar(im, ax=ax, label="relative BVSE (eV)")
        fig.tight_layout()
        if save is not None:
            fig.savefig(save, dpi=300, bbox_inches="tight")
        return fig, ax

    def files(self):
        return {
            "cube": self.cube_path,
            "landscape_metadata": self.landscape_metadata_path,
            "pathway_dir": self.pathway_dir,
            "neb_dir": self.neb_dir,
        }
