# Results: tree volume, DBH and species from ground, drone and aircraft scans

Kaspar von Klitzing with Claude (Anthropic), November 2026. CC BY 4.0.
All tests use cross-validation: each tree (or 100 m block) is predicted by a model fitted without it.

## Data

| Data set | What it gives | Size | Used for |
| --- | --- | --- | --- |
| Demol et al. 2021, Belgium | ground scans + weighed wood, tape DBH, felled height | 65 trees (beech, pine, ash, larch) | volume against real wood |
| BioDiv-3DTrees SEW04_GD, Schorfheide | ground + drone scans of the same trees, leaf-off | 114 trees (59 beech, 54 pine, 1 oak) | drone tests |
| BioDiv-3DTrees AEW33_GD, Swabian Alb | drone + ground scans | 221 trees (206 spruce, 15 broadleaf) | transfer to other site and species |
| LGB laser scans Brandenburg | state flights, leaf-off, 21–24 pts/m², intensity | 4–5 tiles of 1 km², 2019–2024 | species on real state data |
| Tree species map Germany (FU Berlin) | dominant species per 10 m, 77 % accurate | Brandenburg cut-out | labels for LGB tiles |

## Wood volume from ground scans

| Method (65 felled trees) | within ±25 % | within ±10 % | mean deviation per tree | total |
| --- | --- | --- | --- | --- |
| TreeQSM v2.3 (dataset authors) | 43 | – | 26 % | +22 % |
| Kreis-QSM (own) | 54 | – | 14 % | −16 % |
| RayExtract + measured stem | 64 | 33 | 10.6 % | −6 % |
| + simple factor (conifer ×1.15, broadleaf ×0.95) | 65 | 53 | 6.3 % | −1 % |
| **+ fine-twig correction (pipe model, conifers), broadleaf ×0.95** | **65** | **52** | **6.1 %** | **−1 %** |
| Forestry formula (form factor × DBH² × height) | – | 31 | 11.3 % | – |

- Conifer crowns lose fine twigs in the scan (pine −14 %, larch −18 %). Missing volume above each cut branch tip grows with the cube of its radius; on average +13 %, capped at +25 %, tips thicker than 5 cm ignored.
- Broadleaf is already slightly high; ×0.95.
- Dense stands: when only the lower stem is visible, the volume above the last stem circle is the mean of a cone and RayExtract's own taper. Error stays within −3 to −5 % with the stem visible to 15–50 % of height.
- About 5 % per tree is probably the floor: the weighed reference itself has a few percent error.

## DBH from ground scans

- 0.6 cm against tape (62 of 65 clean fits). Pine reads about 0.9 cm thin (tape spans the bark furrows); with that offset 0.44 cm.
- Circle and girth must agree, otherwise the value is rejected (caught a forked ash 13 cm off).
- Stem hidden at 1.3 m: DBH from 0.6–2.5 m extrapolated, equally accurate (0.58 cm). If that fails, estimated from height and crown and flagged.
- Spruce plot (Alb, 221 trees, not felled): 1.3 cm and 0.5 m from the dataset's own values.

## DBH from drone points

A tilted, tapering cylinder through all stem points from 0.5–4 m (about 100 points) instead of one slice (4–8 points).

| Drone DBH vs ground scan | Trees | Mean error | Bias |
| --- | --- | --- | --- |
| Schorfheide, no checks | 97 | 7.2 cm | +5.8 cm |
| Schorfheide, with plausibility checks | 82 | 2.4 cm | +0.8 cm |
| Alb spruce, branch-stub filter | 117 | 2.5 cm | 0.0 cm |

Rejected trees get DBH from height and crown, flagged as estimated (all trees together about 2.7 cm).

## Tree position

Each tree's position is the centre of its lowest clean stem circle (ground scan) or of a tilted cylinder
through the stem points 0.5–4 m (drone), converted to WGS84 when the scan carries a coordinate system.

Ground vs drone scan of the same 112 trees (SEW04_GD, `treemetrics/validate_position.py`):

| | Median | 90 % within | 95 % within |
| --- | --- | --- | --- |
| all trees | 6 cm | 28 cm | 47 cm |
| DBH > 20 cm (94 trees) | 6 cm | 13 cm | – |
| DBH ≤ 20 cm (18 trees) | 19 cm | 55 cm | – |

68 % of trees agree within 10 cm, 88 % within 25 cm; no systematic shift (under 3 cm). Pine 5 cm, beech 10 cm.
Worst cases are thin or leaning beeches with few drone stem points. 2 trees had too few drone points.

What this does and does not show: both scans are in the same plot frame, so this measures how consistently
the method finds the stem across scan types, not absolute map accuracy. The absolute position is only as good
as the scan's own georeferencing. For comparison, BioDiv-3DTrees matched its ground-scan trees to the
Exploratories field inventory at a median distance of 14 cm (all 4,952 trees) and 27 cm (this plot).

## Volume and DBH from height and crown (what an aircraft sees)

| Mean deviation per tree | Beech BE | Pine BE | Pine Schorfheide | Beech Schorfheide |
| --- | --- | --- | --- | --- |
| height only | 21 % | 58 % | 38 % | 93 % |
| height + crown width | 18 % | 20 % | 26 % | 56 % |
| DBH + height (field) | 9 % | 7 % | 18 % | 30 % |

Stand totals within ±6 %, but only with local calibration: Belgian fits overestimate the Schorfheide by +49 % (pine) and +63 % (beech).

## Single trees from above

Pine about 75 % found, leafless beech about 40 %, the same at 1, 20 and ~200 points/m². Crown shape is the limit, not point density.

## Conifer or broadleaf from laser intensity (leaf-off)

| Using | Schorfheide drone | trained Schorfheide → tested Alb | trained Alb → tested Schorfheide |
| --- | --- | --- | --- |
| intensity only | 95 % | 92 % | 97 % |
| crown shape only | 95 % | 47 % | 37 % |
| everything | 96 % | 83 % | 70 % |

Intensity transfers between plots and conifer species; crown shape does not (spruce narrow, pine flat).

LGB state flights (about 23,600 crown cells, balanced accuracy):

| Using | 33465-5780 | 33424-5870 | 33425-5871 | 33417-5805 |
| --- | --- | --- | --- | --- |
| intensity | 94 % | 95 % | 72 % | 67 % |
| everything | 95 % | 97 % | 86 % | 73 % |
| trained on the other three tiles | 89 % | 74 % | 61 % | 50 % |

Raw ground intensity differs by a factor of 8 between flights and one flight strip was 20 % brighter: normalise per strip and calibrate per campaign.

## Deep learning compared

FSCT stem points gave the same DBH as our fit (0.80 vs 0.82 cm). DetailView species: 81 % on summer ground scans, 31–48 % on winter and drone scans; adapted with 3–20 known trees per species it reached 85–98 %, only within the same plot.

## Brandenburg nDOM pipeline

~30,500 tiles, 8–11 s each. Tree count corrected with laser on a random sample of 50 tiles (`brandenburg/results/kalibrierung.json`); growing stock scaled to BWI 2022. Results per tile in `brandenburg/results/kacheln_korrigiert.csv`, map in `uebersicht.png`.

## Limits

- One felled-tree data set; drone points already segmented by the ground scan; species reference 77 % accurate; aerial volume calibration only within one site.

## Sources

- Demol, Gielen & Verbeeck (2021), Zenodo 4557401, CC BY 4.0.
- Griese, Ritzert & Nölke (2025), BioDiv-3DTrees, doi:10.25625/8PB1IF, CC BY 4.0.
- Koeltzow (2025), Tree species map of Germany, Zenodo 15655451.
- LGB laser scan data Brandenburg, © GeoBasis-DE/LGB, dl-de/by-2-0.
- Devereux et al. (2026), RayExtract, Remote Sensing of Environment 334.
- Reutebuch, Andersen & McGaughey (2005), Journal of Forestry 103(6): 286–292.
