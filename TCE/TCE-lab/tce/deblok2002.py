"""
Échantillon de galaxies naines et LSB de de Blok & Bosma (2002, A&A 385, 816) : 24 galaxies,
courbes de rotation « traitées » (directement utilisables pour la modélisation de masse).

Source : CDS J/A+A/385/816, fichier procurve.dat. Colonnes : nom, r["], r[kpc], Vgas (M_gas = 1,4 M_HI),
Vdisk (M/L = 1 en bande R), Vrot, erreur (mesure + inclinaison + asymétrie). Les vitesses sont déjà
corrigées de l'inclinaison.

Attention : le rapport masse/luminosité est défini en bande R : utiliser
`sparc.upsilon_prior(centre, sigma_dex)` avec un centre adapté (vers 1 en R, et non 0,5 à 3,6 micron).
Les distances et inclinaisons ne sont pas fournies : pas de marginalisation possible.
"""
from __future__ import annotations

import gzip
import os
import re
import urllib.request

import numpy as np

from .sparc import Galaxy

URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/385/816/procurve.dat.gz"
FILE = "procurve.dat"


def _num0(field: str) -> float:
    """Champ « ?=0.00 » du CDS : vide signifie 0."""
    field = field.strip()
    return float(field) if field else 0.0


def clean_name(raw: str) -> str:
    """« LSBCF563-1 » -> « F563-1 » (préfixe de catalogue LSBC retiré)."""
    return re.sub(r"^LSBC", "", raw.strip())


def download(dest: str) -> None:
    os.makedirs(dest, exist_ok=True)
    path = os.path.join(dest, FILE)
    if not os.path.exists(path):
        with urllib.request.urlopen(URL, timeout=60) as resp:
            data = gzip.decompress(resp.read())
        with open(path, "wb") as f:
            f.write(data)


def load(dest: str) -> list[Galaxy]:
    rows: dict[str, list[list[float]]] = {}
    with open(os.path.join(dest, FILE), encoding="latin-1") as f:
        for line in f:
            if len(line.rstrip("\n")) < 44:
                continue
            name = clean_name(line[0:11])
            try:
                vals = [float(line[11:16]), float(line[16:22]), _num0(line[22:29]),
                        _num0(line[29:36]), float(line[36:44]), float(line[44:51])]
            except ValueError:
                continue
            rows.setdefault(name, []).append(vals)
    out = []
    for name, v in rows.items():
        d = np.array(v)
        out.append(Galaxy(name=name, r=d[:, 1], vobs=d[:, 4], err=np.maximum(d[:, 5], 1.0),
                          vgas=d[:, 2], vdisk=d[:, 3], vbul=np.zeros(len(d))))
    return out
