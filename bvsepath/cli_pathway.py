import argparse
import json
from .pathway import analyze


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cube")
    p.add_argument("--outdir", default="pathway")
    p.add_argument("--max-relative", type=float, default=3.0)
    p.add_argument("--tile-radius", type=int, default=3)
    p.add_argument("--max-repeat", type=int, default=3)
    a = p.parse_args()
    result = analyze(a.cube, a.outdir, a.max_relative, a.tile_radius, a.max_repeat)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
