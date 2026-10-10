"""
Échantillon indépendant de galaxies à faible brillance de surface (LSB) : de Blok, McGaugh
& Rubin (2001, AJ 122, 2396) et McGaugh, Rubin & de Blok (2001, AJ 122, 2381).

Source : https://astroweb.case.edu/ssm/data/RCsmooth.0701.dat (courbes de rotation hybrides
Halpha-HI lissées + composantes baryoniques radiales). Colonnes : R["], Rkpc, Vgas (M_gas = 1,4 M_HI),
Vdisk (M/L = 1 en bande R), Vbul (M/L = 1 en bande R), V, Err. Les galaxies sans modèle de masse
(« échantillon II ») ont Vgas = Vdisk = Vbul = 0 et sont ignorées.

Attention : le rapport masse/luminosité est défini en bande R, pas à 3,6 micron : utiliser
`sparc.upsilon_prior(centre, sigma_dex, ...)` pour changer la prior sur Upsilon.
"""
from __future__ import annotations

import os
import urllib.request

import numpy as np

from .sparc import Galaxy

URL = "https://astroweb.case.edu/ssm/data/RCsmooth.0701.dat"
FILE = "RCsmooth.0701.dat"


def download(dest: str) -> None:
    os.makedirs(dest, exist_ok=True)
    path = os.path.join(dest, FILE)
    if not os.path.exists(path):
        with urllib.request.urlopen(URL, timeout=60) as resp, open(path, "wb") as f:
            f.write(resp.read())


def load(dest: str) -> list[Galaxy]:
    galaxies, name, rows = [], None, []

    def flush():
        if name and rows:
            d = np.array(rows)
            if np.any(d[:, 2] != 0) or np.any(d[:, 3] != 0):  # modèle de masse disponible
                galaxies.append(Galaxy(name=name, r=d[:, 1], vobs=d[:, 5], err=np.maximum(d[:, 6], 1.0),
                                       vgas=d[:, 2], vdisk=d[:, 3], vbul=d[:, 4]))

    with open(os.path.join(dest, FILE), encoding="latin-1") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            tok = line.split()
            if len(tok) == 1 or (len(tok) >= 1 and not _is_number(tok[0]) and len(tok) < 7 and tok[0] not in ('"',)):
                if tok[0] in ("Name", '"'):
                    continue
                flush()
                name, rows = tok[0], []
            elif len(tok) == 7 and all(_is_number(t) for t in tok):
                rows.append([float(t) for t in tok])
    flush()
    return galaxies


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False
