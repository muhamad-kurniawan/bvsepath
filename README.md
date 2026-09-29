# bvsepath

`bvsepath` generates bond-valence site-energy landscapes for a selected mobile ion and extracts periodic 1D, 2D, and 3D migration pathways from the resulting 3D energy landscape.

The workflow is separated into two commands:

```text
structure -> BVSE cube
BVSE cube -> minima, saddles, directional paths, 1D/2D/3D Ecrit
```

## Installation

```bash
git clone <repository-url>
cd bvsepath
pip install -e .
```

`elements.dat` and `bvse.dat` are bundled with the package and used automatically.

## Energy landscape

### Na+

```bash
bvse-landscape CONTCAR \
  --mobile Na \
  --mobile-valence 1 \
  --oxidation examples/oxidation_na.json \
  --resolution 0.18 \
  --sf 0.74 \
  --output na_bvse.cube
```

### Li+

```bash
bvse-landscape POSCAR \
  --mobile Li \
  --mobile-valence 1 \
  --oxidation examples/oxidation_li.json \
  --resolution 0.18 \
  --sf 0.74 \
  --output li_bvse.cube
```

### Mg2+

```bash
bvse-landscape POSCAR \
  --mobile Mg \
  --mobile-valence 2 \
  --oxidation examples/oxidation_mg.json \
  --resolution 0.18 \
  --sf 0.74 \
  --output mg_bvse.cube
```

The default pair source is `table`, which reads the bundled `bvse.dat`.

The same interface can be used for K+, Ca2+, Ag+, Zn2+, and other ions when the required elemental and mobile-ion/anion pair entries exist in the bundled parameter tables.

Oxidation states can be supplied as a JSON file or directly:

```bash
--oxidation "Li=1,La=3,Zr=4,O=-2"
```

The grid is generated from a real-space resolution:

```text
Ni = ceil(|ai| / resolution)
```

`--sf auto` enables the included pressure-based screening-factor estimator:

```bash
bvse-landscape CONTCAR \
  --mobile Na \
  --mobile-valence 1 \
  --oxidation examples/oxidation_na.json \
  --resolution 0.18 \
  --sf auto \
  --output na_bvse.cube
```

The automatic screening-factor routine is an independent implementation and is not claimed to reproduce the internal softBV `sf_auto` implementation exactly.

Advanced users can override the bundled parameter tables with `--elements` and `--bvse`.

An optional Na-specific mode is available with:

```bash
--pair-source na-published
```

The default is:

```bash
--pair-source table
```

## Pathway analysis

```bash
bvse-pathway na_bvse.cube --outdir pathway
```

Outputs:

```text
pathway/
  minima.csv
  saddles.csv
  path_a.csv
  path_b.csv
  path_c.csv
  summary.json
```

For each pair of adjacent local minima, the conservative local barrier is

```text
wij = max(Es - Ei, Es - Ej)
```

The periodic critical barrier for translation vector `t` is

```text
Ecrit(t) = min_P max_(ij in P) wij
```

If several paths have the same `Ecrit`, the representative path is selected by the shortest total path length and then the fewest hops.

The analyzer evaluates periodic translations in three dimensions. `summary.json` reports:

```text
Ecrit_a_eV
Ecrit_b_eV
Ecrit_c_eV
Ecrit_1D_eV
Ecrit_2D_eV
Ecrit_3D_eV
```

`Ecrit_1D`, `Ecrit_2D`, and `Ecrit_3D` are the minimum edge-barrier thresholds at which the periodic saddle graph contains translation cycles of rank 1, rank 2, and rank 3, respectively.

For difficult three-dimensional networks:

```bash
bvse-pathway bvse.cube \
  --outdir pathway \
  --tile-radius 3 \
  --max-relative 4.0
```

No fully sodiated parent structure is required. Local minima and saddle connections are obtained directly from the 3D energy landscape.

## Recommended convergence check

Evaluate representative structures at several grid resolutions, for example:

```text
0.20 A
0.18 A
0.15 A
```

Verify that both `Ecrit` and the critical saddle are stable.

## Parameter data

The bundled `elements.dat` and `bvse.dat` are taken from the CAVD release repository:

Shuhebing/CAVD, release branch  
https://gitee.com/shuhebing/cavd/tree/release

If these parameter tables are used in published work, cite the CAVD source together with the relevant bond-valence and BVSE literature below.

## Method scope

The reported `Ecrit` values are BVSE-derived screening descriptors and should not be described as DFT or NEB activation energies. Selected hops can be validated independently by MLIP-NEB or DFT-NEB.

The energy-landscape calculation follows the bond-valence site-energy framework. The pathway analyzer is an independent local-minimum, basin-saddle, and periodic minimax implementation and is not the official softBV/BVPA pathway algorithm.

## References

H. Chen, L. L. Wong, and S. Adams, “SoftBV – a software tool for screening the materials genome of inorganic fast ion conductors,” *Acta Crystallographica Section B* **75**, 18–33 (2019). https://doi.org/10.1107/S2052520618015718

H. Chen and S. Adams, “Bond softness sensitive bond-valence parameters for crystal structure plausibility tests,” *IUCrJ* **4**, 614–625 (2017). https://doi.org/10.1107/S2052252517010211

L. L. Wong, K. C. Phuah, R. Dai, H. Chen, W. S. Chew, and S. Adams, “Bond Valence Pathway Analyzer—An Automatic Rapid Screening Tool for Fast Ion Conductors within softBV,” *Chemistry of Materials* **33**, 625–641 (2021). https://doi.org/10.1021/acs.chemmater.0c03893

B. He et al., “High-throughput screening platform for solid electrolytes combining hierarchical ion-transport prediction algorithms,” *Scientific Data* **7**, 151 (2020). https://doi.org/10.1038/s41597-020-0474-y

CAVD release repository and parameter tables: https://gitee.com/shuhebing/cavd/tree/release
