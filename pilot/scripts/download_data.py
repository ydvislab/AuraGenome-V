"""Download the pilot's ground-truth data (not redistributed in this repository).

Kircher M, et al. (2019) Saturation mutagenesis of twenty disease-associated regulatory elements
at single base-pair resolution. Nat Commun 10:3583. doi:10.1038/s41467-019-11526-w
Source: the data of the official MPRA portal (github.com/kircherlab/MPRA_SaturationMutagenesis),
pinned to a commit and checked against a SHA-256 so every run uses identical input.
Any use of these data must cite the publication above.

usage: python download_data.py <pilot_dir>
"""
import hashlib, os, sys, urllib.request

COMMIT = "16a48622a20754b16ee356d2c4e214459ffe1476"
FILES = {
    "elements.tsv.gz": "fec2eed91fe27af3aae07ebce2eca65e9bad4bb6abba5d8c27f478887dd7b134",
}
BASE = f"https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/{COMMIT}/data/"

raw = os.path.join(sys.argv[1], "data", "raw")
os.makedirs(raw, exist_ok=True)
for name, sha in FILES.items():
    path = os.path.join(raw, name)
    if not os.path.exists(path):
        print(f"downloading {name}")
        urllib.request.urlretrieve(BASE + name, path)
    got = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if got != sha:
        sys.exit(f"checksum mismatch for {name}: {got}")
    print(f"{name} ok ({os.path.getsize(path)} bytes)")
