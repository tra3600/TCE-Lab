"""
Ajustement de Chrono-Energetics sur la base SPARC (Lelli, McGaugh & Schombert 2016,
175 galaxies à disque, courbes de rotation HI/Hα + photométrie Spitzer 3,6 µm).

Données : http://astroweb.cwru.edu/SPARC/ (téléchargées par `download_sparc`,
non versionnées). Fichiers utilisés : `Rotmod_LTG.zip` (courbes de rotation et
composantes baryoniques à Upsilon = 1) et `SPARC_Lelli2016c.mrt` (Table 1).

Modèles comparés (accélération radiale g = v^2 / r, en (km/s)^2/kpc)
--------------------------------------------------------------------
newton   : g = g_bar
mond     : g = nu(g_bar/a0) g_bar, interpolation « simple »
mond_a0  : contrôle, MOND avec a0 libre par galaxie (« rc » = a0 dans la sortie)
tce      : g = g_bar sqrt((1 + x^2 a0/g_bar) / (1 + x^2)), x = r/rc, rc libre par galaxie
tce_k    : même forme avec rc = k * Rdisk, k unique pour toutes les galaxies
tce_v1   : variante de l'éq. (5), g = g_bar sqrt(1 + w a0/g_bar), w = x^2/(1+x^2)
tce_v2   : variante cohérente RG (G_eff = 1/R_v), g = g_bar + w sqrt(g_bar a0)
(suffixe _k : rc = k * Rdisk)

Paramètres libres par galaxie : Upsilon_disque (prior log-normale, 0,1 dex autour de
0,5 ; Upsilon_bulbe = 1,4 Upsilon_disque) et, pour `tce`, rc. Les incertitudes de
distance et d'inclinaison ne sont pas modélisées (limite connue).
"""
from __future__ import annotations

import io
import os
import urllib.request
import zipfile
from dataclasses import dataclass

import numpy as np

BASE_URL = "https://astroweb.cwru.edu/SPARC/"
ZIP_NAME = "Rotmod_LTG.zip"
TABLE_NAME = "SPARC_Lelli2016c.mrt"
KPC_M = 3.0856775814913673e19
#: a0 = 1.2e-10 m/s^2 converti en (km/s)^2/kpc
A0_KMS2_KPC = 1.2e-10 * KPC_M / 1e6
UPSILON_DISK_PRIOR = 0.5
UPSILON_BULGE_RATIO = 1.4
UPSILON_SIGMA_DEX = 0.1


@dataclass
class Galaxy:
    name: str
    r: np.ndarray        # kpc
    vobs: np.ndarray     # km/s
    err: np.ndarray      # km/s
    vgas: np.ndarray
    vdisk: np.ndarray    # à Upsilon = 1
    vbul: np.ndarray     # à Upsilon = 1
    distance: float = np.nan
    rdisk: float = np.nan   # kpc
    quality: int = 0
    vflat: float = np.nan


# --------------------------------------------------------------------------
# Données
# --------------------------------------------------------------------------
def download_sparc(dest: str) -> None:
    """Télécharge les deux fichiers SPARC dans `dest` (HTTPS)."""
    os.makedirs(dest, exist_ok=True)
    for name in (ZIP_NAME, TABLE_NAME):
        path = os.path.join(dest, name)
        if not os.path.exists(path):
            with urllib.request.urlopen(BASE_URL + name, timeout=60) as resp, open(path, "wb") as f:
                f.write(resp.read())


def _parse_table1(path: str) -> dict:
    rows = {}
    with open(path, encoding="latin-1") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 18:
                try:
                    float(parts[2])
                    int(parts[1])
                except ValueError:
                    continue
                rows[parts[0]] = {
                    "D": float(parts[2]), "rdisk": float(parts[11]),
                    "vflat": float(parts[15]), "Q": int(parts[17]),
                }
    return rows


def load_sparc(directory: str) -> list[Galaxy]:
    table = _parse_table1(os.path.join(directory, TABLE_NAME))
    galaxies = []
    with zipfile.ZipFile(os.path.join(directory, ZIP_NAME)) as z:
        for member in sorted(z.namelist()):
            if not member.endswith("_rotmod.dat"):
                continue
            name = os.path.basename(member).replace("_rotmod.dat", "")
            data = np.array(
                [[float(v) for v in ln.split()] for ln in
                 io.TextIOWrapper(z.open(member), encoding="latin-1")
                 if ln.strip() and not ln.startswith("#")]
            )
            meta = table.get(name, {})
            galaxies.append(Galaxy(
                name=name, r=data[:, 0], vobs=data[:, 1], err=data[:, 2],
                vgas=data[:, 3], vdisk=data[:, 4], vbul=data[:, 5],
                distance=meta.get("D", np.nan), rdisk=meta.get("rdisk", np.nan),
                quality=meta.get("Q", 0), vflat=meta.get("vflat", np.nan),
            ))
    return galaxies


# --------------------------------------------------------------------------
# Modèles
# --------------------------------------------------------------------------
def g_baryon(gal: Galaxy, upsilon_disk):
    """g_bar (nU, n) pour un tableau de Upsilon_disque."""
    u = np.atleast_1d(upsilon_disk)[:, None]
    v2 = (gal.vgas * np.abs(gal.vgas))[None, :] + u * (gal.vdisk * np.abs(gal.vdisk))[None, :] \
        + UPSILON_BULGE_RATIO * u * (gal.vbul * np.abs(gal.vbul))[None, :]
    return v2 / gal.r[None, :]


def g_newton(gbar, r, **_):
    return gbar


def g_mond_simple(gbar, r, a0=A0_KMS2_KPC, **_):
    y = np.maximum(gbar, 1e-12) / a0
    return gbar * (0.5 + np.sqrt(0.25 + 1.0 / y))


def g_tce(gbar, r, rc, a0=A0_KMS2_KPC, **_):
    """Lecture B de l'éq. (5) + g_eff = g_N sqrt(R0 / R_v)."""
    x2 = (r / rc) ** 2
    gb = np.maximum(gbar, 1e-12)
    return gb * np.sqrt((1.0 + x2 * a0 / gb) / (1.0 + x2))


def _w(r, rc):
    x2 = (r / rc) ** 2
    return x2 / (1.0 + x2)


def g_tce_v1(gbar, r, rc, a0=A0_KMS2_KPC, **_):
    """Variante V1 (couplage racine carrée) : g = g_bar sqrt(1 + w a0/g_bar)."""
    gb = np.maximum(gbar, 1e-12)
    return gb * np.sqrt(1.0 + _w(r, rc) * a0 / gb)


def g_tce_v2(gbar, r, rc, a0=A0_KMS2_KPC, **_):
    """Variante V2 (couplage G_eff = 1/R_v, cohérent RG) : g = g_bar + w sqrt(g_bar a0)."""
    gb = np.maximum(gbar, 1e-12)
    return gb + _w(r, rc) * np.sqrt(gb * a0)


#: modèles dépendant de rc
RC_MODELS = {"tce": g_tce, "tce_v1": g_tce_v1, "tce_v2": g_tce_v2}


def velocity(g, r):
    return np.sqrt(np.maximum(g, 0.0) * r)


# --------------------------------------------------------------------------
# Ajustement par grille
# --------------------------------------------------------------------------
UPSILON_GRID = np.logspace(np.log10(0.2), np.log10(1.2), 41)
RC_GRID = np.logspace(-1.0, 2.5, 141)  # 0,1 à ~300 kpc


def _prior_chi2(u):
    return ((np.log10(u) - np.log10(UPSILON_DISK_PRIOR)) / UPSILON_SIGMA_DEX) ** 2


def _good(gal: Galaxy):
    """Points exploitables : données valides ET g_bar > 0 pour le Upsilon du prior
    (la composante gazeuse signée peut rendre g_bar <= 0 au centre ; ces points sont
    exclus pour TOUS les modèles afin de garder la comparaison équitable)."""
    ok = (gal.r > 0) & (gal.vobs > 0) & (gal.err > 0)
    gb = (g_baryon(gal, UPSILON_DISK_PRIOR)[0]) if ok.any() else np.zeros_like(gal.r)
    return ok & (gb > 0)


def fit_galaxy(gal: Galaxy, model: str, a0: float = A0_KMS2_KPC, rc_fixed: float | None = None):
    """Meilleur ajustement sur la grille ; renvoie un dict (chi2, dof, Upsilon, rc...)."""
    m = _good(gal)
    r, vo, er = gal.r[m], gal.vobs[m], gal.err[m]
    g = Galaxy(gal.name, r, vo, er, gal.vgas[m], gal.vdisk[m], gal.vbul[m])
    gbar = g_baryon(g, UPSILON_GRID)  # (nU, n)
    prior = _prior_chi2(UPSILON_GRID)
    if model == "newton":
        gm = gbar[:, None, :]
        rcs = np.array([np.nan])
    elif model == "mond":
        gm = g_mond_simple(gbar, r, a0)[:, None, :]
        rcs = np.array([np.nan])
    elif model == "mond_a0":
        # contrôle : MOND avec a0 libre par galaxie (un paramètre libre de plus, comme rc)
        rcs = A0_KMS2_KPC * np.logspace(-0.7, 0.7, 57)
        gm = g_mond_simple(gbar[:, None, :], r[None, None, :], rcs[None, :, None])
    elif model.removesuffix("_k") in RC_MODELS:
        rcs = np.array([rc_fixed]) if rc_fixed is not None else RC_GRID
        gm = RC_MODELS[model.removesuffix("_k")](gbar[:, None, :], r[None, None, :], rcs[None, :, None], a0)
    else:
        raise ValueError(model)
    v = velocity(gm, r)
    chi2_data = (((v - vo) / er) ** 2).sum(axis=-1)          # (nU, nR)
    chi2 = chi2_data + prior[:, None]
    iu, ir = np.unravel_index(np.argmin(chi2), chi2.shape)
    k = 1 + (model in RC_MODELS or model == "mond_a0")
    return {
        "name": gal.name, "model": model, "chi2": float(chi2[iu, ir]),
        "chi2_data": float(chi2_data[iu, ir]), "n": int(len(r)), "dof": int(len(r) - k),
        "upsilon": float(UPSILON_GRID[iu]), "rc": float(rcs[ir]),
        "g_obs": vo**2 / r, "g_bar": gbar[iu], "g_mod": gm[iu, ir],
    }


def fit_all(galaxies, model, a0=A0_KMS2_KPC, k_rdisk: float | None = None, quality_max: int = 2):
    out = []
    for gal in galaxies:
        if gal.quality and gal.quality > quality_max:
            continue
        if model.endswith("_k"):
            if not np.isfinite(gal.rdisk) or gal.rdisk <= 0:
                continue
            out.append(fit_galaxy(gal, model, a0, rc_fixed=k_rdisk * gal.rdisk))
        else:
            out.append(fit_galaxy(gal, model, a0))
    return out


def fit_k_global(galaxies, a0=A0_KMS2_KPC, ks=np.logspace(-0.5, 1.5, 21), quality_max=2, model="tce_k"):
    """Cherche le k unique (rc = k Rdisk) qui minimise le chi2 total."""
    totals = []
    for k in ks:
        res = fit_all(galaxies, model, a0, k_rdisk=k, quality_max=quality_max)
        totals.append(sum(r["chi2"] for r in res))
    i = int(np.argmin(totals))
    return float(ks[i]), np.array(totals)


def summarize(results, n_global_params: int = 0):
    chi2 = np.array([r["chi2"] for r in results])
    dof = np.array([r["dof"] for r in results])
    n = np.array([r["n"] for r in results])
    k_free = (n - dof).sum() + n_global_params
    bic = chi2.sum() + k_free * np.log(n.sum())
    resid = np.concatenate([np.log10(np.maximum(r["g_obs"], 1e-12) / np.maximum(r["g_mod"], 1e-12)) for r in results])
    return {
        "galaxies": len(results), "points": int(n.sum()), "chi2_total": float(chi2.sum()),
        "chi2_red_median": float(np.median(chi2 / dof)), "chi2_red_total": float(chi2.sum() / dof.sum()),
        "BIC": float(bic), "rms_dex": float(np.sqrt(np.mean(resid**2))),
    }
