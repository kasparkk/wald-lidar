# wald-lidar

Open laser scan data → single trees with height, crown, DBH, wood volume and conifer/broadleaf.
Built on open Brandenburg state data (LGB) and the BioDiv-3DTrees and Demol et al. datasets,
tested against felled and weighed trees.

Kaspar von Klitzing, with Claude (Anthropic), 2026. Code MIT, texts and result tables CC BY 4.0.

## What's in here

| Folder | What it does |
| --- | --- |
| [`treemetrics/`](treemetrics) | Per-tree values for already segmented trees (one point cloud per tree): branch model via RayExtract, corrected wood volume, DBH, height, position, conifer/broadleaf from intensity. Works behind any segmentation, e.g. 3Dtrees.earth. |
| [`treeflow/`](treeflow) | Point cloud in, single trees out: ground model, tree detection (own or lidR), measurements, report, compact 3D web viewer. |
| [`brandenburg/`](brandenburg) | Statewide pipeline over all ~30,500 1 km tiles of the open LGB canopy height model (nDOM): trees, stand height, crown cover, rough growing stock; laser calibration of the tree count; conifer/broadleaf from colour-infrared orthophotos; the Hammer 3D viewer. |
| [`dbh-test/`](dbh-test) | DBH-from-height fits checked on 1,368 ground-scanned trees (BioDiv-3DTrees, Schorfheide). |
| [`downloads/`](downloads) | Windows scripts that fetch the third-party data used (nothing third-party is stored here). |
| [`docs/`](docs) | Summary of results and methods. |

## Main results

Full numbers in [docs/RESULTS.md](docs/RESULTS.md).

| | Ground scan (TLS) | Drone (ULS) | Aircraft (ALS, LGB) |
| --- | --- | --- | --- |
| Height | yes | 0.3–0.7 m vs ground scan | about 1 m |
| DBH | 0.6 cm vs tape (0.44 cm with pine bark offset) | about 2–2.5 cm, trees > 20 cm | from height only: 5–7 cm |
| Wood volume per tree | about 6 % vs weighed wood (65 trees) | no | stand level only, needs local calibration |
| Conifer vs broadleaf (leaf-off) | possible | 92–97 %, transfers between sites | 94–97 % on 2 of 4 tiles, calibrate per flight |

Live viewer: ~800,000 clickable trees around Hammer, Brandenburg: https://wald-hammer-3d.netlify.app

## Quick start

```
pip install numpy scipy pandas laspy[lazrs] scikit-learn
python -m treeflow run SCAN.laz --out results              # inside treeflow/
python treemetrics/run_trees.py TREES/ --out results --bin ~/raycloudtools/build/raycloudtools
python brandenburg/bb_wald.py --out ergebnisse --bbox 405 5774 412 5781
```

RayExtract (RayCloudTools, CSIRO) must be built separately for the branch models: https://github.com/csiro-robotics/raycloudtools

## Limits

- Volume is validated on one felled-tree dataset (65 Belgian trees, one scanner).
- Drone tests used point clouds already cut into trees by the ground-scan segmentation.
- Species reference for the LGB tiles is itself a map with 77 % accuracy.
- The Brandenburg growing stock is a rough value scaled to the national forest inventory (BWI 2022), not a measurement.

See [CREDITS.md](CREDITS.md) for data and tool sources and their licences.
