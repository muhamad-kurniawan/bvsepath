import argparse
import json
from .landscape import generate


def main():
    p = argparse.ArgumentParser()
    p.add_argument("structure")
    p.add_argument("--mobile", required=True)
    p.add_argument("--oxidation", required=True)
    p.add_argument("--resolution", type=float, default=0.18)
    p.add_argument("--grid-shape", default=None)
    p.add_argument("--sf", default="0.74")
    p.add_argument("--sf-min", type=float, default=0.50)
    p.add_argument("--sf-max", type=float, default=0.85)
    p.add_argument("--sf-strain", type=float, default=1.0e-3)
    p.add_argument("--output", default="bvse.cube")
    p.add_argument("--metadata", default="landscape.json")
    p.add_argument("--ions", default=None)
    p.add_argument("--pairs", default=None)
    a = p.parse_args()
    grid_shape = None if a.grid_shape is None else [int(x) for x in a.grid_shape.split(",")]
    result = generate(a.structure, a.mobile, a.oxidation, a.output, a.resolution, a.sf, a.ions, a.pairs, a.metadata, sf_min=a.sf_min, sf_max=a.sf_max, sf_strain=a.sf_strain, grid_shape=grid_shape)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
