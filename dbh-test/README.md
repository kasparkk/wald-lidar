# DBH from tree height

Checks how well DBH follows from height alone (what an aircraft scan sees) on the ground-scanned
trees of BioDiv-3DTrees (`labels.csv`, columns `dTLS`, `hTLS`, `cpaHTLS`; download from
https://doi.org/10.25625/8PB1IF, CC BY 4.0). Leave-one-out power-law fits per species group.

    python dbhfit.py labels.csv

Schorfheide (plot group SEW), mean absolute error:

| Group | Trees | Old viewer rule (1.5·h − 6) | Height fit | Height + crown width |
| --- | --- | --- | --- | --- |
| Pine | 512 | 5.4 cm (bias −3 cm) | 4.9 cm | 3.9 cm |
| Broadleaf | 856 | 7.4 cm | 6.8 cm | – |

Fits now used in the Hammer viewer: pine DBH = 0.9796·h^1.105, broadleaf DBH = 0.2172·h^1.546
(cm, m). With crown diameter cd: pine DBH = 1.1277·h^0.829·cd^0.444.
