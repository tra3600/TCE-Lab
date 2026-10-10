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
mond_rs  : contrôle équitable, MOND (a0 fixé) avec rayon de transition libre rs (« portail »)
mond_n   : contrôle, MOND (a0 fixé) avec indice d'interpolation libre n
tce      : g = g_bar sqrt((1 + x^2 a0/g_bar) / (1 + x^2)), x = r/rc, rc libre par galaxie
tce_k    : même forme avec rc = k * Rdisk, k unique pour toutes les galaxies
tce_v1   : variante de l'éq. (5), g = g_bar sqrt(1 + w a0/g_bar), w = x^2/(1+x^2)
tce_v2   : variante cohérente RG (G_eff = 1/R_v), g = g_bar + w sqrt(g_bar a0)
(suffixe _k : rc = k * Rdisk)

Paramètres libres par galaxie : Upsilon_disque (prior log-normale, 0,1 dex autour de
0,5 ; Upsilon_bulbe = 1,4 Upsilon_disque) et, pour `tce`, rc. Les incertitudes de
distance et d'inclinaison sont modélisées avec `marginalize=True` (priors gaussiens
 sur D et i, grille 7 x 7 sur +-2 sigma).
"""
from __future__ import annotations

import contextlib
import io
import os
import re
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
    e_dist: float = np.nan   # Mpc
    inc: float = np.nan      # degrés
    e_inc: float = np.nan    # degrés
    lum: float = np.nan      # L[3.6], 1e9 L_sun
    mhi: float = np.nan      # masse HI, 1e9 M_sun
    reff: float = np.nan     # kpc
    sbdisk: float = np.nan   # brillance centrale du disque, L_sun/pc^2


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
                    "D": float(parts[2]), "e_D": float(parts[3]),
                    "inc": float(parts[5]), "e_inc": float(parts[6]),
                    "L": float(parts[7]), "reff": float(parts[9]),
                    "rdisk": float(parts[11]), "sbdisk": float(parts[12]),
                    "mhi": float(parts[13]),
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
                e_dist=meta.get("e_D", np.nan), inc=meta.get("inc", np.nan),
                e_inc=meta.get("e_inc", np.nan), lum=meta.get("L", np.nan),
                mhi=meta.get("mhi", np.nan), reff=meta.get("reff", np.nan),
                sbdisk=meta.get("sbdisk", np.nan),
            ))
    return galaxies


def load_rotmod_directory(directory: str) -> list[Galaxy]:
    """Charge un dossier de fichiers `*_rotmod.dat` au format SPARC (colonnes :
    Rad[kpc] Vobs errV Vgas Vdisk Vbul SBdisk SBbul ; en-tête « # Distance = X Mpc »).
    Sert à tester le modèle sur un échantillon indépendant (ex. THINGS, de Blok+2008) dont
    les courbes baryoniques ont été converties à ce format. Les méta-données (inclinaison,
    erreur de distance...) sont inconnues : la marginalisation est alors sans effet."""
    out = []
    for fn in sorted(os.listdir(directory)):
        if not fn.endswith("_rotmod.dat"):
            continue
        dist = np.nan
        rows = []
        with open(os.path.join(directory, fn), encoding="latin-1") as f:
            for ln in f:
                if ln.startswith("#"):
                    mt = re.search(r"Distance\s*=\s*([0-9.]+)", ln)
                    if mt:
                        dist = float(mt.group(1))
                elif ln.strip():
                    rows.append([float(v) for v in ln.split()])
        if not rows:
            continue
        d = np.array(rows)
        out.append(Galaxy(name=fn.replace("_rotmod.dat", ""), r=d[:, 0], vobs=d[:, 1], err=d[:, 2],
                          vgas=d[:, 3], vdisk=d[:, 4], vbul=d[:, 5] if d.shape[1] > 5 else np.zeros(len(d)),
                          distance=dist))
    return out


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


def g_mond_n(gbar, r=None, n=1.0, a0=A0_KMS2_KPC, **_):
    """MOND avec indice d'interpolation n : g = g_bar nu_n(g_bar/a0),
    nu_n = [(1 + sqrt(1 + 4 y^-n)) / 2]^(1/n) ; n = 1 : « simple », n = 2 : « standard »."""
    gb = np.maximum(gbar, 1e-12)
    y = gb / a0
    return gb * ((1.0 + np.sqrt(1.0 + 4.0 * y ** (-n))) / 2.0) ** (1.0 / n)


def g_mond_gated(gbar, r, rs, a0=A0_KMS2_KPC, **_):
    """MOND « à portail » : g = g_bar + w (g_MOND - g_bar), w = x^2/(1+x^2), x = r/rs.
    Même structure que TCE V2 (g_bar + w sqrt(g_bar a0)) mais avec l'interpolation MOND."""
    return gbar + _w(r, rs) * (g_mond_simple(gbar, r, a0) - gbar)


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


@contextlib.contextmanager
def upsilon_prior(center: float, sigma_dex: float, lo: float | None = None, hi: float | None = None):
    """Change temporairement la prior log-normale sur le rapport masse/luminosité du disque
    (par défaut : 0,5 à 3,6 micron). La grille couvre [lo, hi] (défaut : center/2.5 à center*2.4)."""
    global UPSILON_DISK_PRIOR, UPSILON_SIGMA_DEX, UPSILON_GRID
    saved = (UPSILON_DISK_PRIOR, UPSILON_SIGMA_DEX, UPSILON_GRID)
    lo = center / 2.5 if lo is None else lo
    hi = center * 2.4 if hi is None else hi
    UPSILON_DISK_PRIOR, UPSILON_SIGMA_DEX = center, sigma_dex
    UPSILON_GRID = np.logspace(np.log10(lo), np.log10(hi), 41)
    try:
        yield
    finally:
        UPSILON_DISK_PRIOR, UPSILON_SIGMA_DEX, UPSILON_GRID = saved


def _prior_chi2(u):
    return ((np.log10(u) - np.log10(UPSILON_DISK_PRIOR)) / UPSILON_SIGMA_DEX) ** 2


def _good(gal: Galaxy):
    """Points exploitables : données valides ET g_bar > 0 pour le Upsilon du prior
    (la composante gazeuse signée peut rendre g_bar <= 0 au centre ; ces points sont
    exclus pour TOUS les modèles afin de garder la comparaison équitable)."""
    ok = (gal.r > 0) & (gal.vobs > 0) & (gal.err > 0)
    gb = (g_baryon(gal, UPSILON_DISK_PRIOR)[0]) if ok.any() else np.zeros_like(gal.r)
    return ok & (gb > 0)


#: nombre de sigmas explorés pour la distance et l'inclinaison (grille de 7 points sur +-2 sigma)
NUISANCE_Z = np.linspace(-2.0, 2.0, 7)


#: contrôles MOND à un paramètre libre par galaxie (la sortie « rc » contient ce paramètre)
CONTROL_MODELS = ("mond_a0", "mond_rs", "mond_n")


def _base_model(model: str) -> str:
    return model.removesuffix("_k").removesuffix("_law")


def _nuisance_grid(gal: Galaxy):
    """Liste (facteur de distance, rapport sin i0/sin i, chi2 de prior) ; un seul
    point neutre si distance ou inclinaison sont inconnues."""
    ok = np.isfinite([gal.distance, gal.e_dist, gal.inc, gal.e_inc]).all() and gal.distance > 0
    if not ok:
        return [(1.0, 1.0, 0.0)]
    e_inc = max(gal.e_inc, 1.0)
    out = []
    for zd in NUISANCE_Z:
        fD = 1.0 + zd * gal.e_dist / gal.distance
        if fD < 0.3:
            continue
        for zi in NUISANCE_Z:
            inc = np.clip(gal.inc + zi * e_inc, 5.0, 89.0)
            ratio = np.sin(np.radians(gal.inc)) / np.sin(np.radians(inc))
            out.append((fD, ratio, zd**2 + zi**2))
    return out


def fit_galaxy(gal: Galaxy, model: str, a0: float = A0_KMS2_KPC, rc_fixed: float | None = None,
               marginalize: bool = False):
    """Meilleur ajustement sur la grille ; renvoie un dict (chi2, dof, Upsilon, rc...).

    marginalize=True ajoute deux paramètres de nuisance avec prior gaussienne :
      * la distance D' = D (1 + z e_D/D) : rayons r -> r D'/D et vitesses baryoniques
        V_bar -> V_bar sqrt(D'/D) (M ~ D^2, V^2 ~ M/r) ;
      * l'inclinaison i' = i + z e_i : V_obs et erreurs -> x sin(i)/sin(i').
    Le chi2 renvoyé inclut les priors (maximum a posteriori)."""
    m = _good(gal)
    base = _base_model(model)
    combos = _nuisance_grid(gal) if marginalize else [(1.0, 1.0, 0.0)]
    n_nuis = 2 if (marginalize and len(combos) > 1) else 0
    prior_u = _prior_chi2(UPSILON_GRID)
    best = None
    for fD, ratio, prior_nuis in combos:
        sq = np.sqrt(fD)
        r = gal.r[m] * fD
        vo, er = gal.vobs[m] * ratio, gal.err[m] * ratio
        g = Galaxy(gal.name, r, vo, er, gal.vgas[m] * sq, gal.vdisk[m] * sq, gal.vbul[m] * sq)
        gbar = g_baryon(g, UPSILON_GRID)  # (nU, n)
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
        elif model == "mond_rs":
            # contrôle équitable : MOND (a0 fixé) avec un rayon de transition libre
            rcs = RC_GRID
            gm = g_mond_gated(gbar[:, None, :], r[None, None, :], rcs[None, :, None], a0)
        elif model == "mond_n":
            # contrôle : MOND (a0 fixé) avec indice d'interpolation libre
            rcs = np.logspace(np.log10(0.3), np.log10(8.0), 57)
            gm = g_mond_n(gbar[:, None, :], n=rcs[None, :, None], a0=a0)
        elif base in RC_MODELS:
            rcs = np.array([rc_fixed * fD]) if rc_fixed is not None else RC_GRID
            gm = RC_MODELS[base](gbar[:, None, :], r[None, None, :], rcs[None, :, None], a0)
        else:
            raise ValueError(model)
        v = velocity(gm, r)
        chi2_data = (((v - vo) / er) ** 2).sum(axis=-1)          # (nU, nR)
        chi2 = chi2_data + prior_u[:, None] + prior_nuis
        iu, ir = np.unravel_index(np.argmin(chi2), chi2.shape)
        if best is None or chi2[iu, ir] < best["chi2"]:
            best = {
                "name": gal.name, "model": model, "chi2": float(chi2[iu, ir]),
                "chi2_data": float(chi2_data[iu, ir]), "n": int(len(r)),
                "upsilon": float(UPSILON_GRID[iu]), "rc": float(rcs[ir]) / (fD if rc_fixed is not None else 1.0),
                "dist_factor": float(fD), "inc_ratio": float(ratio),
                "g_obs": vo**2 / r, "g_bar": gbar[iu], "g_mod": gm[iu, ir],
            }
    free_model = model in RC_MODELS or model in CONTROL_MODELS
    best["dof"] = best["n"] - (1 + n_nuis + free_model)
    return best


MIN_POINTS = 8  # galaxies avec moins de points valides exclues (dof <= 0 avec 4 paramètres libres)


def fit_all(galaxies, model, a0=A0_KMS2_KPC, k_rdisk: float | None = None, quality_max: int = 2,
            rc_fn=None, marginalize: bool = False):
    """rc fixé par galaxie : soit rc = k_rdisk * Rdisk (modèles « _k »), soit rc = rc_fn(galaxie)
    (modèles « _law »). Sinon rc est libre pour les modèles TCE."""
    out = []
    for gal in galaxies:
        if gal.quality and gal.quality > quality_max:
            continue
        if _good(gal).sum() < MIN_POINTS:
            continue
        if rc_fn is not None:
            rc = rc_fn(gal)
            if not np.isfinite(rc) or rc <= 0:
                continue
            out.append(fit_galaxy(gal, model, a0, rc_fixed=rc, marginalize=marginalize))
        elif model.endswith("_k"):
            if not np.isfinite(gal.rdisk) or gal.rdisk <= 0:
                continue
            out.append(fit_galaxy(gal, model, a0, rc_fixed=k_rdisk * gal.rdisk, marginalize=marginalize))
        else:
            out.append(fit_galaxy(gal, model, a0, marginalize=marginalize))
    return out


def fit_k_global(galaxies, a0=A0_KMS2_KPC, ks=np.logspace(-0.5, 1.5, 21), quality_max=2,
                 model="tce_k", marginalize: bool = False):
    """Cherche le k unique (rc = k Rdisk) qui minimise le chi2 total."""
    totals = []
    for k in ks:
        res = fit_all(galaxies, model, a0, k_rdisk=k, quality_max=quality_max, marginalize=marginalize)
        totals.append(sum(r["chi2"] for r in res))
    i = int(np.argmin(totals))
    return float(ks[i]), np.array(totals)


# --------------------------------------------------------------------------
# Lois candidates pour rc (surface densité baryonique)
# --------------------------------------------------------------------------
G_KPC = 4.30091e-6        # kpc (km/s)^2 / M_sun
UPSILON_STAR = UPSILON_DISK_PRIOR  # pour estimer M_bar et Sigma_b
HELIUM_FACTOR = 1.33
#: Sigma_dagger = a0 / G, surface densité critique de MOND (M_sun/pc^2)
SIGMA_DAGGER = A0_KMS2_KPC / G_KPC / 1e6


def baryonic_mass(gal: Galaxy) -> float:
    """M_bar = Upsilon* L[3.6] + 1.33 M_HI, en M_sun (Upsilon* = 0,5 à 3,6 micron)."""
    return (UPSILON_STAR * gal.lum + HELIUM_FACTOR * gal.mhi) * 1e9


def central_surface_density(gal: Galaxy) -> float:
    """Sigma_b = Upsilon* x brillance centrale du disque, en M_sun/pc^2."""
    return UPSILON_STAR * gal.sbdisk


def rc_law_rdisk(gal, k):
    """L1 : rc = k Rdisk."""
    return k * gal.rdisk


def rc_law_mass(gal, k, a0=A0_KMS2_KPC):
    """L2 : rc = k sqrt(G M_bar / a0), rayon où g_bar = a0."""
    return k * np.sqrt(G_KPC * baryonic_mass(gal) / a0)


def rc_law_surface(gal, k, alpha):
    """L3 : rc = k Rdisk (Sigma_b / Sigma_dagger)^alpha ; alpha = 0 redonne L1."""
    return k * gal.rdisk * (central_surface_density(gal) / SIGMA_DAGGER) ** alpha


LAWS = {
    "L1_rdisk": (rc_law_rdisk, {"k": np.logspace(-1.5, 1.5, 25)}),
    "L2_masse": (rc_law_mass, {"k": np.logspace(-1.5, 1.5, 25)}),
    "L3_surface": (rc_law_surface, {"k": np.logspace(-1.5, 1.5, 25), "alpha": np.linspace(-1.5, 1.5, 13)}),
}


def law_chi2_table(galaxies, base_model: str, law: str, a0=A0_KMS2_KPC, quality_max: int = 2,
                   marginalize: bool = False):
    """chi2 par galaxie sur toute la grille de paramètres globaux de la loi.
    Renvoie (grille de paramètres [liste de dicts], noms de galaxies, matrice chi2 [n_params, n_gal])."""
    fn, grids = LAWS[law]
    keys = list(grids)
    mesh = np.meshgrid(*[grids[k] for k in keys], indexing="ij")
    params = [dict(zip(keys, vals)) for vals in zip(*[m.ravel() for m in mesh])]
    sel = [g for g in galaxies if not (g.quality and g.quality > quality_max) and _good(g).sum() >= MIN_POINTS
           and np.isfinite(g.rdisk) and g.rdisk > 0 and np.isfinite(g.lum) and np.isfinite(g.sbdisk) and g.sbdisk > 0]
    table = np.empty((len(params), len(sel)))
    for ip, prm in enumerate(params):
        for ig, gal in enumerate(sel):
            rc = fn(gal, **prm)
            table[ip, ig] = fit_galaxy(gal, base_model + "_law", a0, rc_fixed=rc, marginalize=marginalize)["chi2"] \
                if np.isfinite(rc) and rc > 0 else np.inf
    return params, [g.name for g in sel], table


def law_best(params, table, mask=None):
    """Paramètres globaux minimisant le chi2 total (éventuellement sur un sous-ensemble)."""
    tot = table.sum(axis=1) if mask is None else table[:, mask].sum(axis=1)
    i = int(np.argmin(tot))
    return (params[i] if params is not None else None), i, float(tot[i])


def cross_validate_law(table, seed: int = 0, folds: int = 2):
    """Validation croisée : paramètres ajustés sur les autres plis, chi2 évalué sur le pli retenu.
    Renvoie le chi2 total hors échantillon (somme sur tous les plis)."""
    n = table.shape[1]
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    held = 0.0
    for f in range(folds):
        test = np.zeros(n, dtype=bool)
        test[idx[f::folds]] = True
        _, i, _ = law_best(None, table, mask=~test)
        held += table[i, test].sum()
    return float(held)


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
