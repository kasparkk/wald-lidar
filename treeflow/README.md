# treeflow – point cloud in, single trees out

An automatic chain that turns any laser scan of a forest (airborne, drone, handheld or terrestrial;
LAS or LAZ) into a table of single trees with position, height, crown and – where the stem is
visible – measured diameter. Built and tested on LGB Brandenburg data; meant to also run on
3Dtrees.earth / BioDiv-3DTrees point clouds.

## Run it

```
pip install laspy[lazrs] numpy scipy scikit-image scikit-learn pyproj rasterio pyyaml
python -m treeflow run SCAN.laz --out results
python -m treeflow run SCAN.laz --out results --bbox 410215 5785815 410275 5785875 --ortho DOP.tif
python -m treeflow run FOLDER_WITH_LAZ --out results           # every file, plus overview.csv
python -m treeflow run SCAN.laz --out results --reference trees_ref.csv   # scores the result
```

## Show already segmented trees in 3D (the application layer)

Trees that someone has already found – 3Dtrees.earth, SegmentAnyTree, FSCT, lidR or treeflow – become
compact 3D tree models for a web viewer:

```
python -m treeflow visualise SEGMENTED.laz --out bundle --site mysite --title "My forest" [--labels values.csv]
python -m treeflow visualise FOLDER_OF_SINGLE_TREES --out bundle --site plots      # one file per tree
```

- Finds the tree id field by itself (`treeID`, `PredInstance`, `instance`, `tree_id`, …) or take `--tree-field`.
- Per tree: crown as rings (6 heights × 12 directions) from its own points, stem from circle fits
  every metre where the scan shows it, height, crown base and width, DBH at 1.3 m if measured.
- `--labels` adds any values per tree id (species, source, NDVI …); they appear in the viewer and as
  colour options.
- Output is tiled (`tile_i_j.json` models, `pts_i_j.json` a few points per tree, `terrain.json`).
  Copy the bundle to `viewer/data/` and open `viewer/index.html` through any web server.
- To get a segmented file from treeflow itself: `write_segmented: true` in the settings.

## Use lidR for the tree detection (recommended for airborne data)

treeflow can hand step 5 to **lidR**, the established R package for airborne tree segmentation,
and keep doing everything else (ground model, measurements, photo values, checks, outputs).

```
# once: install R (cran.r-project.org), then in R:  install.packages("lidR")
python -m treeflow run SCAN.laz --out results --config lidr.yaml      # lidr.yaml: detector: lidr-dalponte
python -m treeflow compare SCAN.laz --out cmp --bbox ... [--reference trees_ref.csv]
```

`compare` runs the same area through treeflow's own detector and three lidR methods
(Dalponte 2016, Silva 2016, Li 2012) and writes `compare.csv`.

## What happens, step by step

| Step | What it does | Decided automatically by |
|---|---|---|
| 1 Read | Loads the cloud (optionally a box plus a 10 m margin) | – |
| 2 Describe | Point density → airborne (< 60/m²), drone (< 1,500/m²) or ground scan | density |
| 3 Ground | Uses the file's ground class; if missing, runs its own ground filter | ground class present? |
| 4 Normalise | Terrain model (1 m) → every point gets height above ground | – |
| 5 Find trees | **From above:** canopy model → tree tops with a height-dependent search window → crowns by watershed. **From below:** slice at 1.0–1.6 m → stem clusters → circle fit → each stem seeds a tree | scan type and whether enough stems are found |
| 6 Measure | Height, crown base, crown area/diameter, crown volume, DBH (stem fit only), point counts | – |
| 7 Photo (optional) | NDVI per crown from an RGB+NIR image; flags crowns that differ from the stand | `--ortho` given |
| 8 Check (optional) | Matches trees to reference trees; detection rate, height and DBH error | `--reference` given |
| 9 Write | `trees.csv`, `trees.geojson` (WGS84), `summary.json`, `report.html`, `viewer.json` | – |

Every DBH says where it came from (`dbh_source`). Estimates are only added with `estimate_dbh: true`.

## Test results so far

| Data | Trees found | Wrong detections | Height error | DBH error |
|---|---|---|---|---|
| Synthetic ground scan (53 trees, 8 under the canopy) | 92 % | 0 % | −0.34 m bias, RMSE 0.37 m | RMSE 0.5 cm (49 stems) |
| Synthetic airborne scan (same stand) | 77 % (all misses under the canopy or merged) | 0 % | RMSE 0.13 m | not measurable |
| LGB tile 33410-5785, 60 × 60 m pine stand | 78 trees (217/ha) | unknown – no reference yet | unknown | not measurable |
| LGB tile, 300 × 300 m | 2,876 trees in 26 s | unknown | unknown | – |

**Own detector vs lidR**

| Data | treeflow | lidR Dalponte | lidR Silva | lidR Li 2012 |
|---|---|---|---|---|
| Synthetic airborne, trees found | 77 % | 68 % | 68 % | 74 % |
| LGB 60 × 60 m patch, trees | 78 | 65 | 65 | 80 |

On the real patch, 64 of lidR-Dalponte's 65 trees are also in treeflow's result; treeflow's 14 extra
trees are most likely wide pine crowns split in two (visible on the aerial photo). Which result is
right needs reference trees.

The synthetic stands are idealised (perfect cylinders, no wind, no occlusion), so real errors will be
larger. On the real patch some wide pine crowns are split into two tops; `window_a`/`window_b`
tune that, but the right setting needs reference trees.

## Using 3Dtrees data

- **BioDiv-3DTrees** (Univ. Göttingen, Scientific Data 2025): 4,952 single-tree TLS/ULS clouds of
  19 species with species labels and QSMs, on GROdata, DOI 10.25625/8PB1IF. Use: run the stem route
  on the trees and compare DBH and height with their labels; later, train a species model on them.
- **3Dtrees.earth** (Univ. Freiburg): upload a scan, export their segmentation, save it as a CSV with
  `x, y, height_m, dbh_cm`, then run treeflow on the same scan with `--reference` to compare both
  methods tree by tree.

## Files

`treeflow/io.py` read and describe · `ground.py` ground filter and terrain model · `detect.py`
canopy model, tree tops, stem fitting, segmentation · `metrics.py` per-tree values · `enrich.py`
aerial photo · `validate.py` comparison with reference trees · `export.py` outputs · `pipeline.py`
the chain · `cli.py` command line · `tests/make_synthetic.py` test stands with known answers.
