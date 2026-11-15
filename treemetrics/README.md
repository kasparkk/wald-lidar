# treemetrics: tree models, volume, DBH and species for segmented trees

Adds per-tree values to trees that are already segmented (one point cloud per tree): a light branch model,
wood volume, DBH, height, and for leaf-off scans conifer vs broadleaf. Input is just the points of a single
tree, so it works behind any segmentation.

Kaspar von Klitzing, with Claude (Anthropic), November 2026. Summary of the results:
https://claude.ai/artifact/9YyyMLRpiD4kEqnBxbL2CU

## What it gives, by scan type

| Scan | Height, crown | DBH | Species (leaf-off) | Branch model and volume |
| --- | --- | --- | --- | --- |
| Ground / backpack (TLS, MLS) | yes | about 0.6 cm (0.44 cm with pine bark offset) | possible | yes, about 6 % per tree after correction (MLS not yet tested) |
| Drone (ULS) | yes, height to 0.3-0.7 m | larger trees only, about 2-2.5 cm | 92-97 % conifer vs broadleaf, transfers to new species (intensity) | no |
| Aircraft (ALS) | yes | no | 94-97 % on 2 of 4 test tiles, needs calibration per flight | no |

## Files

| File | Does |
| --- | --- |
| `run_trees.py` | folder of single-tree LAS/LAZ files -> `trees.csv` + one branch model per tree (ground scans) |
| `stem.py` | measures the stem (circles every 50 cm, smooth taper) and puts it into the RayExtract model |
| `treevol.py` | overlap-free volume at forks, ground height at the stem base, DBH from circle and girth at 1.3 m |
| `corrections.py` | twig correction for conifers, broadleaf factor, DBH plausibility check and pine bark offset; all constants with their source |
| `validate_felled.py` | the test against 65 felled and weighed trees; also refits the corrections on your own felled trees |
| `drone_dbh.py` | DBH from drone points: tilted tapering cylinder 0.5-4 m plus plausibility checks; stub trimming for spruce |
| `species.py` | conifer vs broadleaf from echo intensity (default; crown shape does not transfer to new species); per-strip intensity normalisation |

## Setup

Python 3.9+ with `pip install numpy scipy pandas scikit-learn laspy[lazrs]`.
The ground-scan models need RayCloudTools (CSIRO), built on Linux or WSL:

```
sudo apt install -y build-essential cmake git libeigen3-dev
git clone --depth 1 https://github.com/norlab-ulaval/libnabo
cd libnabo && mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DLIBNABO_BUILD_TESTS=OFF -DLIBNABO_BUILD_EXAMPLES=OFF -DLIBNABO_BUILD_PYTHON=OFF && make -j4 && sudo make install && cd ~
git clone --depth 1 https://github.com/csiro-robotics/raycloudtools
cd raycloudtools && mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j4
```

## Use

```
python3 run_trees.py TREES/ --out results --bin ~/raycloudtools/build/raycloudtools --species species.csv
python3 validate_felled.py DEMOL/ --bin ~/raycloudtools/build/raycloudtools
```

`species.csv` has columns `tree,species` (Latin names). Without it you get everything except the corrected volume.
A few seconds per tree. Models are written as `x y z radius parent` (metres, parent -1 = root), small enough
to show a whole stand in a browser.

## How it works

- **RayExtract** with `--crop_length 0.25` (default 1 m cuts off crown wood). Single-tree files have no scanner
  positions, so rays are faked from above and a flat ground mesh is added under the tree.
- **Measured stem:** RANSAC circles every 50 cm from 0.8 m above the uphill ground, a smooth taper curve through
  them replaces RayExtract's stem radii; a cone from the last circle to the top.
  In dense plots the stem is often visible only part of the way up; a cone alone then underestimates (hiding
  the upper stem of the felled trees: -5, -11, -18 % at 50, 30, 15 % of the height). Below half the height the
  volume is therefore the mean of the cone and RayExtract's own taper scaled to the last circle (-4, -5, -3 %).
  This costs a little on fully measured stems (6.3 instead of 6.1 %). On the Schorfheide plot 61 of 114 stems
  were measured to less than half the height, in Belgium 2 of 65. RayExtract alone makes the stem
  4.7 cm too thin at breast height, with the measured stem 2.4 cm.
- **Volume** without overlaps: side branches count from the surface of the branch they grow from.
- **Ground** at the stem base from the lowest stem points in twelve sectors, uphill side (upper quartile).
- **DBH** at 1.3 m: mean of RANSAC circle and girth (convex hull perimeter / pi, like a tape), when the circle is
  well supported and the stem is covered. Otherwise projected to 1.3 m from slices between 0.6 and 2.5 m
  (as accurate as the direct value on the Belgian trees; about 1.3 cm off for stems scanned from one side only);
  the column `dbh_method` says which. On the Schorfheide plot: 99 direct, 12 projected, 3 not measurable.
- **Corrections:** conifers miss fine twigs; by the pipe model the missing volume above each cut tip grows with
  r^3, so `added volume = 2656 x sum(r_tip^3)`. Broadleaf models are already slightly high: x 0.94.

## Results on 65 felled trees (Demol et al. 2021; leave-one-out)

| Volume against weighed wood | within 25 % | within 10 % | mean deviation | total |
| --- | --- | --- | --- | --- |
| RayExtract | 59 | 37 | 10.9 % | -9 % |
| RayExtract + measured stem | 64 | 32 | 10.7 % | -6 % |
| + corrections | 65 | 52 | 6.3 % | -1 % |
| classic form factor x DBH^2 x height, tape DBH | - | 35 | 9.8 % | - |

## Limits

- One data set with weighed wood (65 Belgian trees, one scanner). The RayExtract authors used the same trees,
  defaults may be tuned to them. Correction constants and the pine bark offset are not yet checked elsewhere;
  refit them with `validate_felled.py` when you have felled trees.
- The weighed reference itself has a few percent of error; about 5 % per tree is probably the floor.
- Species and drone DBH were tested on two drone plots (BioDiv-3DTrees SEW04_GD pine/beech, AEW33_GD spruce) and four LGB tiles, leaf-off only.
  Drone DBH rule 2 was set after looking at failures.
- The test harness used cached RayExtract models; the RayExtract call itself is the same code as in earlier runs.

## Data and methods

- Demol, Gielen & Verbeeck (2021): QSMs, point cloud and harvest data from a destructive forest biomass experiment in Belgium, Zenodo 4557401, CC BY 4.0.
- Griese, Ritzert & Nölke (2025): BioDiv-3DTrees, GRO.data, doi:10.25625/8PB1IF, CC BY 4.0.
- Devereux et al. (2026): RayExtract, Remote Sensing of Environment 334, 115162; RayCloudTools, CSIRO.
- LGB Laserscandaten Brandenburg, © GeoBasis-DE/LGB, dl-de/by-2-0. Koeltzow (FU Berlin, 2025): Tree species map of Germany, Zenodo 15655451.
