# Data and tools

No third-party raw data is stored in this repository. Derived tables and viewer data carry the
attributions below.

## Data

| Source | Used for | Licence |
| --- | --- | --- |
| Landesvermessung und Geobasisinformation Brandenburg (LGB): [laser scans, nDOM, DOP20 RGBI](https://geobasis-bb.de/lgb/de/geodaten/3d-produkte/laserscandaten/) | Brandenburg pipeline, species tests, viewer | © GeoBasis-DE/LGB, [dl-de/by-2-0](https://www.govdata.de/dl-de/by-2-0) |
| Demol, Gielen & Verbeeck (2021): [QSMs, point cloud and harvest data, Belgium](https://zenodo.org/records/4557401), Zenodo 4557401 | Volume and DBH validation against weighed wood | CC BY 4.0 |
| Griese, Ritzert & Nölke (2025): [BioDiv-3DTrees](https://doi.org/10.25625/8PB1IF), GRO.data | Ground and drone scans of single trees, DBH test | CC BY 4.0 |
| Koeltzow (FU Berlin, 2025): [Tree species map of Germany](https://zenodo.org/records/15655451), Zenodo 15655451 | Species labels for LGB tiles | see record |
| Bundeswaldinventur 2022 (Thünen-Institut) | Scaling the rough growing stock | public statistics |

## Tools

- [RayCloudTools / RayExtract](https://github.com/csiro-robotics/raycloudtools), CSIRO: branch models. Devereux et al. (2026), Remote Sensing of Environment 334.
- [TreeQSM](https://github.com/InverseTampere/TreeQSM), [rTwig](https://aidanmorales.github.io/rTwig/), [FSCT](https://github.com/SKrisanski/FSCT), [DetailView](https://github.com/JulFrey/DetailView): compared against.
- [lidR](https://github.com/r-lidar/lidR): optional tree detection in treeflow.
- laspy, NumPy, SciPy, pandas, scikit-learn, three.js.

## Literature

- Reutebuch, Andersen & McGaughey (2005): Light Detection and Ranging (LIDAR): An Emerging Tool for Multiple Resource Inventory. Journal of Forestry 103(6): 286–292.
- Jucker et al. (2017): Allometric equations for integrating remote sensing imagery into forest monitoring programmes. Global Change Biology 23: 177–190.
