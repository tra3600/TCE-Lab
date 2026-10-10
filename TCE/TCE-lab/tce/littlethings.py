"""
Échantillon indépendant LITTLE THINGS (Oh et al. 2015, AJ 149, 180 ; 26 galaxies naines HI).

Source : CDS J/AJ/149/180 (table1.dat, table2.dat). Seules les quantités intégrées sont
disponibles ici (masse de gaz, masse stellaire, vitesse au rayon maximal) : le CDS ne
fournit pas les courbes baryoniques radiales (Vgas, Vstar). Ce module permet donc un test
de la relation de Tully-Fisher baryonique (éq. 7 du manuscrit), pas un ajustement radial.

Unités : masses en M_sun (le fichier les donne en 1e7 M_sun), vitesses en km/s, rayons en kpc.
Colonnes de table2.dat (séparées par « | ») : 0 Nom, 1 Rmax, 2 R0.3, 3 V(Rmax), 4 Viso(Rmax),
..., 24 Mgas, 25 MstarK, 26 MstarSED (1e7 M_sun), 27 log Mdyn, 28 log M200.
"""
from __future__ import annotations

import gzip
import os
import re
import urllib.request
from dataclasses import dataclass

import numpy as np

from . import sparc

BASE_URL = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/149/180/"
FILES = ("table1.dat", "table2.dat")


@dataclass
class LTGalaxy:
    name: str
    distance: float      # Mpc
    rmax: float          # kpc
    v_rmax: float        # km/s (dérive asymétrique corrigée)
    v_iso: float         # km/s
    m_gas: float         # M_sun
    m_star_sed: float    # M_sun (nan si absent)
    m_star_kin: float    # M_sun (nan si absent)

    @property
    def m_star(self) -> float:
        """Masse stellaire : SED (Zhang+2012) si disponible, sinon cinématique."""
        return self.m_star_sed if np.isfinite(self.m_star_sed) else self.m_star_kin

    @property
    def m_bar(self) -> float:
        return self.m_gas + self.m_star


def normalize_name(name: str) -> str:
    """Normalise un nom de galaxie pour comparer SPARC et LITTLE THINGS
    (« DDO_154 » = « DDO154 », « UGC08508 » = « UGC8508 », « U5750 » = « UGC5750 »)."""
    s = re.sub(r"[^A-Z0-9]", "", name.upper())
    s = re.sub(r"^U(\d)", r"UGC\1", s)          # « U5750 » = « UGC 5750 »
    return re.sub(r"^([A-Z]+)0+(\d)", r"\1\2", s)


def download(dest: str) -> None:
    os.makedirs(dest, exist_ok=True)
    for name in FILES:
        path = os.path.join(dest, name)
        if not os.path.exists(path):
            with urllib.request.urlopen(BASE_URL + name + ".gz", timeout=60) as resp:
                data = gzip.decompress(resp.read())
            with open(path, "wb") as f:
                f.write(data)


def _num(s: str) -> float:
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return float("nan")


def load(dest: str) -> list[LTGalaxy]:
    dist = {}
    with open(os.path.join(dest, "table1.dat"), encoding="latin-1") as f:
        for line in f:
            p = line.split("|")
            if len(p) > 2:
                dist[p[0].strip()] = _num(p[2])
    out = []
    with open(os.path.join(dest, "table2.dat"), encoding="latin-1") as f:
        for line in f:
            p = line.split("|")
            if len(p) < 27:
                continue
            name = p[0].strip()
            out.append(LTGalaxy(
                name=name, distance=dist.get(name, float("nan")), rmax=_num(p[1]),
                v_rmax=_num(p[3]), v_iso=_num(p[4]),
                m_gas=_num(p[24]) * 1e7, m_star_kin=_num(p[25]) * 1e7, m_star_sed=_num(p[26]) * 1e7,
            ))
    return out


# --------------------------------------------------------------------------
# Relation de Tully-Fisher baryonique, éq. (7) : v^4 = G M_bar a0
# --------------------------------------------------------------------------
def btfr_velocity(m_bar, a0: float = sparc.A0_KMS2_KPC):
    """v = (G M_bar a0)^(1/4), en km/s (G en kpc (km/s)^2/M_sun, a0 en (km/s)^2/kpc)."""
    return (sparc.G_KPC * np.asarray(m_bar, dtype=float) * a0) ** 0.25


def btfr_residuals(v_obs, m_bar, a0: float = sparc.A0_KMS2_KPC):
    """log10(v_obs / v_pred) pour chaque galaxie."""
    return np.log10(np.asarray(v_obs, dtype=float) / btfr_velocity(m_bar, a0))


def fit_btfr_slope(v_obs, m_bar):
    """Pente et ordonnée de log10 M_bar = s log10 v + b (régression de M_bar sur v,
    adaptée à une erreur dominante sur M_bar). Éq. (7) prédit s = 4."""
    x, y = np.log10(np.asarray(v_obs, dtype=float)), np.log10(np.asarray(m_bar, dtype=float))
    s, b = np.polyfit(x, y, 1)
    resid = y - (s * x + b)
    return float(s), float(b), float(resid.std(ddof=2))


def sparc_btfr_sample(galaxies, quality_max: int = 1, stellar_upsilon: float = 0.5):
    """Échantillon SPARC pour comparaison : Vflat de la table et M_bar = Upsilon L[3.6] + 1,33 M_HI."""
    v, m, names = [], [], []
    for g in galaxies:
        if g.quality and g.quality > quality_max:
            continue
        if not np.all(np.isfinite([g.vflat, g.lum, g.mhi])) or g.vflat <= 0:
            continue
        v.append(g.vflat)
        m.append((stellar_upsilon * g.lum + sparc.HELIUM_FACTOR * g.mhi) * 1e9)
        names.append(g.name)
    return np.array(v), np.array(m), names
