import argparse
import json
from .guide import make_guides


def main():
    p = argparse.ArgumentParser()
    p.add_argument("structure")
    p.add_argument("pathway")
    p.add_argument("--mobile", required=True)
    p.add_argument("--cube", default=None)
    p.add_argument("--directions", default="a,b,c")
    p.add_argument("--outdir", default="neb_guides")
    a = p.parse_args()
    directions = [x.strip() for x in str(a.directions).split(",") if x.strip()]
    result = make_guides(a.structure, a.pathway, a.mobile, cube_path=a.cube, directions=directions, outdir=a.outdir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
