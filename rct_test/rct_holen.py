"""
Lädt drei Laser-Kacheln aus dem Raum Hammer/Löpten für den Test mit RayCloudTools auf 3Dtrees.earth.

Die Kacheln sind bewusst verschieden gewählt:
  33411-5775  dichter, geschlossener Kiefernbestand (Kronenschluss 0,90, kaum Laub)
  33408-5778  lockerer, mit etwas Laubholz (Kronenschluss 0,67, Laub ca. 14 %)
  33409-5779  hoher Bestand (Median 22,4 m), offen, etwas Laubholz

Ergebnis:
  upload/als_<kachel>.laz   -> diese Dateien bei 3Dtrees hochladen
  ndom/NDOM_<kachel>.laz    -> braucht später rct_vergleich.py

Daten: GeoBasis-DE / LGB, Datenlizenz Deutschland Namensnennung 2.0
"""
import argparse, os, sys, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "brandenburg"))
from bb_laser_kalib import ALS, NDOM, get  # noqa: E402

KACHELN = ["33411-5775", "33408-5778", "33409-5779"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kacheln", nargs="*", default=KACHELN)
    a = ap.parse_args()
    up = os.path.join(HERE, "upload"); nd = os.path.join(HERE, "ndom")
    os.makedirs(up, exist_ok=True); os.makedirs(nd, exist_ok=True)
    for key in a.kacheln:
        laz = os.path.join(up, f"als_{key}.laz")
        if not os.path.exists(laz):
            print(f"{key}: Laser laden (ca. 120 MB) ...", flush=True)
            z = get(ALS + f"als_{key}.zip", os.path.join(up, f"als_{key}.zip"))
            with zipfile.ZipFile(z) as zf:
                name = [x for x in zf.namelist() if x.endswith(".laz")][0]
                with zf.open(name) as src, open(laz, "wb") as dst:
                    while True:
                        b = src.read(1 << 20)
                        if not b:
                            break
                        dst.write(b)
            os.remove(z)
        ndp = os.path.join(nd, f"NDOM_{key}.laz")
        if not os.path.exists(ndp):
            print(f"{key}: nDOM laden (ca. 15 MB) ...", flush=True)
            get(NDOM + f"NDOM_{key}_LAS12.laz", ndp)
        print(f"{key}: fertig, {os.path.getsize(laz) / 1e6:.0f} MB", flush=True)
    print("\nZum Hochladen bei 3Dtrees:", up)


if __name__ == "__main__":
    main()
