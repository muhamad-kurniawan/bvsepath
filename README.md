# bvsepath

`bvsepath` is a Python implementation for generating bond-valence site-energy landscapes and extracting periodic ion-migration pathways.

## Installation

```bash
git clone https://github.com/muhamad-kurniawan/bvsepath/
cd bvsepath
pip install -e .
```

## Energy landscape

Only the structure, mobile ion, and oxidation states parameters are required in general use. The default grid resolution is 0.18 A. The default screening factor is fixed at 0.74 for backward-compatible reproducibility, while `--sf auto` enables the independent pressure-proxy optimizer.

```bash
bvse-landscape CONTCAR \
  --mobile Na \
  --oxidation "Na=1,Co=3,Fe=3,Mg=2,Mn=7,Ni=2,O=-2" \
  --output bvse.cube
```

Li example:

```bash
bvse-landscape POSCAR \
  --mobile Li \
  --oxidation "Li=1,La=3,Zr=4,O=-2" \
  --output li_bvse.cube
```

Oxidation states may also be stored in JSON:

```json
{
  "Na": 1,
  "Co": 3,
  "Fe": 3,
  "Mg": 2,
  "Mn": 7,
  "Ni": 2,
  "O": -2
}
```

```bash
bvse-landscape CONTCAR --mobile Na --oxidation oxidation.json --output bvse.cube
```

Optional numerical controls are:

```bash
--resolution 0.18
--sf 0.74
```

Automatic screening-factor optimization can be requested with:

```bash
bvse-landscape CONTCAR \
  --mobile Na \
  --oxidation oxidation.json \
  --sf auto \
  --output bvse.cube
```

### Automatic screening factor: basis and implementation

The `--sf auto` option is an independent implementation inspired by the screening-factor optimization described for softBV. Chen, Wong, and Adams reported that softBV estimates the screening factor by varying it iteratively so that a fast static-pressure relaxation gives a pressure close to zero. The physical idea is to tune the range of the screened same-sign Coulomb repulsion until the repulsive Coulomb contribution and the attractive/short-range Morse contribution are approximately balanced for the structure.

In this package, the same physical idea is used. It uses a finite-difference pressure proxy based only on the occupied mobile-ion environment:

```text
P_proxy(sf) = - [E_mobile(V+) - E_mobile(V-)] / [V+ - V-]
```

with

```text
V- = V(1 - epsilon)^3
V+ = V(1 + epsilon)^3
```

and a default isotropic strain

```text
epsilon = 0.001
```

For each trial screening factor, the code evaluates the occupied-mobile-ion BVSE environment in the slightly compressed and expanded cells while keeping fractional coordinates fixed. It then searches for

```text
P_proxy(sf) = 0
```

within the requested interval:

```bash
--sf auto
--sf-min 0.50
--sf-max 0.85
--sf-strain 0.001
```

If a sign change exists, the root is obtained numerically. If no root occurs inside the interval, the value that minimizes `abs(P_proxy)` is returned. The selected screening factor, search interval, residual pressure proxy, and optimization mode are written to the landscape metadata JSON.

This is therefore a **pressure-balance proxy**, not a reproduction of the internal softBV `sf_auto` routine. The original softBV procedure can include a more complete force-field pressure balance.

The conceptual basis is:

Chen, H., Wong, L. L., and Adams, S. (2019), *SoftBV – a software tool for screening the materials genome of inorganic fast ion conductors*, Acta Crystallographica Section B 75, 18-33. https://doi.org/10.1107/S2052520618015718

A useful comparison is also provided by:

He, B. et al. (2020), *High-throughput screening platform for solid electrolytes combining hierarchical ion-transport prediction algorithms*, Scientific Data 7, 151. https://doi.org/10.1038/s41597-020-0474-y

He et al. use a fixed factor of 0.74 for high-throughput BVSE screening and explicitly contrast it with softBV, where the screening factor is iteratively adapted according to the balance between Morse and Coulomb interactions in the individual structure.

Custom parameter files can be supplied with `--elements` and `--bvse`, but the package includes both files by default.

## Pathway analysis

```bash
bvse-pathway bvse.cube --outdir pathway
```

The output directory contains:

```text
minima.csv
saddles.csv
path_a.csv
path_b.csv
path_c.csv
summary.json
```

For a local connection between minima `i` and `j`, the conservative edge barrier is

```text
w_ij = max(E_s - E_i, E_s - E_j)
```

where `E_s` is the basin-boundary saddle energy.

The directional periodic barrier is

```text
Ecrit(t) = min_P max_(ij in P) w_ij
```

where `P` spans one lattice translation `t`.

The graph search first minimizes `Ecrit`, then total path length, then hop count. This avoids arbitrary detours when several paths share the same bottleneck barrier.

`summary.json` reports

```text
Ecrit_a_eV
Ecrit_b_eV
Ecrit_c_eV
Ecrit_1D_eV
Ecrit_2D_eV
Ecrit_3D_eV
```

The 1D, 2D, and 3D thresholds are determined from the rank of periodic translation cycles that become available as the allowed edge barrier increases.

For a wider three-dimensional search:

```bash
bvse-pathway bvse.cube --outdir pathway --max-relative 6.0 --tile-radius 3
```
## NEB candidate guide (beta version)

This command prepares a small guide package for NEB setup.

```bash
bvse-neb-guide CONTCAR pathway --mobile Na --outdir neb_guides
```

The candidate-selection rule is simple and designed to make manual NEB setup easier:

1. Find the path with the requested periodic translation.
2. Identify the highest-barrier edge along that path.
3. Prefer an edge where exactly one endpoint minimum is occupied by the mobile ion in the input structure.
4. Use the occupied endpoint as the initial site and the unoccupied endpoint as the target site.

This is still in beta version. And it only work on certain problems or structures..

## References

H. Chen, L. L. Wong, and S. Adams, “SoftBV – a software tool for screening the materials genome of inorganic fast ion conductors,” *Acta Crystallographica Section B* **75**, 18–33 (2019). https://doi.org/10.1107/S2052520618015718

H. Chen and S. Adams, “Bond softness sensitive bond-valence parameters for crystal structure plausibility tests,” *IUCrJ* **4**, 614–625 (2017). https://doi.org/10.1107/S2052252517010211

L. L. Wong, K. C. Phuah, R. Dai, H. Chen, W. S. Chew, and S. Adams, “Bond Valence Pathway Analyzer—An Automatic Rapid Screening Tool for Fast Ion Conductors within softBV,” *Chemistry of Materials* **33**, 625–641 (2021). https://doi.org/10.1021/acs.chemmater.0c03893

B. He et al., “High-throughput screening platform for solid electrolytes combining hierarchical ion-transport prediction algorithms,” *Scientific Data* **7**, 151 (2020). https://doi.org/10.1038/s41597-020-0474-y

## Third-party data notice

The software source code and the bundled parameter datasets have different provenance. See `THIRD_PARTY_PARAMETERS.md`.
