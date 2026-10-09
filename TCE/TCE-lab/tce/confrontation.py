"""
Confrontation de Chrono-Energetics avec la relativité générale (RG), MOND et
la matière noire (RG + halo NFW), et liens avec les cadres « émergents ».

Ce module ne modifie pas `chrono_energetics.py` ; il en réutilise les fonctions.

Contenu
-------
1. Horloges : taux d'écoulement du temps prédit par la RG (Schwarzschild et
   champ faible) contre celui qu'implique l'éq. (1), dt ∝ 1/R_v.
2. Couplage gravitationnel : si le vide de rigidité R_v est lu comme
   kappa_eff = 1/R_v (lecture cohérente avec la RG, où R_0 = 1/kappa),
   g_eff = g_N * R_0/R_v. Comparaison avec la racine carrée utilisée par le
   manuscrit (éq. 6).
3. Courbes de rotation : Newton, TCE, MOND (interpolations simple/standard),
   RG + halo NFW.
4. Échelle a0 et cosmologie : a0 ~ c H0 / (2 pi) (Milgrom), c H0 / 6 (Verlinde).
5. Pont avec la thermodynamique de l'espace-temps (Jacobson) :
   R_0 = hbar c^2 eta / (2 pi), eta = 1 / (4 l_P^2).
"""
from __future__ import annotations

import numpy as np

from . import chrono_energetics as ce
from .constants import G, c

HBAR = 1.054571817e-34  # J s
MPC = 3.0856775814913673e22  # m
KPC = MPC / 1e3
MSUN = 1.98847e30  # kg


# --------------------------------------------------------------------------
# 1. Horloges : RG contre éq. (1)
# --------------------------------------------------------------------------
def schwarzschild_radius(M):
    return 2.0 * G * M / c**2


def gr_clock_rate_schwarzschild(r, M):
    """d tau / d t_infini = sqrt(1 - r_s / r) (horloge statique, RG exacte)."""
    r = np.asarray(r, dtype=float)
    return np.sqrt(1.0 - schwarzschild_radius(M) / r)


def gr_clock_rate_weak_field(r, M):
    """Champ faible : d tau / d t = 1 + Phi/c^2, Phi = -G M / r."""
    return 1.0 - G * M / (np.asarray(r, dtype=float) * c**2)


def tce_clock_rate(R_v, R0: float | None = None):
    """Taux d'horloge implicite de l'éq. (1) pour une même énergie dissipée :
    dt(r) / dt(infini) = R_0 / R_v(r)  (ici « infini » = vide non perturbé)."""
    R0 = ce.vacuum_stiffness_R0() if R0 is None else R0
    return R0 / np.asarray(R_v, dtype=float)


def tce_gr_profile_for_black_hole(r, M, R0: float | None = None):
    """Profil de rigidité qui reproduit la dilatation du temps de Schwarzschild
    avec t = E/R : R_v(r) = R_0 / sqrt(1 - r_s/r) (dérivation du PDF, §5).
    C'est un choix qui retombe sur la RG par construction, pas une prédiction."""
    R0 = ce.vacuum_stiffness_R0() if R0 is None else R0
    return R0 / gr_clock_rate_schwarzschild(r, M)


def clock_test_galaxy(r, M_bar, rc, a0: float = ce.A0_MOND):
    """Compare, dans une galaxie, le taux d'horloge de la RG (champ faible) à
    celui que l'éq. (1) implique avec le profil R_v(r) de l'éq. (5)."""
    gN = ce.g_newton(M_bar, r)
    Rv = ce.stiffness_profile(r, rc, gN, a0)
    tce = tce_clock_rate(Rv)
    gr = gr_clock_rate_weak_field(r, M_bar)
    return {
        "r": np.asarray(r, dtype=float),
        "g_N": gN,
        "R_v_over_R0": Rv / ce.vacuum_stiffness_R0(),
        "tce_rate": tce,
        "gr_rate": gr,
        "tce_deviation": tce - 1.0,
        "gr_deviation": gr - 1.0,
        "ratio_of_deviations": (tce - 1.0) / (gr - 1.0),
    }


# --------------------------------------------------------------------------
# 2. Couplage gravitationnel : racine carrée (manuscrit) contre 1/R_v (RG)
# --------------------------------------------------------------------------
def g_eff_sqrt(gN, Rv, R0: float | None = None):
    """Éq. (6) du manuscrit, lue comme g_N * sqrt(R_0 / R_v)."""
    return ce.g_eff_from_stiffness(gN, Rv, R0)


def g_eff_varying_G(gN, Rv, R0: float | None = None):
    """Lecture cohérente avec la RG : G_Einstein = 1/R_v donne G_eff = G R_0/R_v,
    donc g_eff = g_N * R_0 / R_v. À faible accélération, R_v ~ R_0 g_N/a0 donne
    g_eff -> a0 (constante) et v^2 = a0 r : la courbe de rotation MONTE, elle
    n'est pas plate. La racine carrée du manuscrit n'est pas dérivée de 1/R_v."""
    R0 = ce.vacuum_stiffness_R0() if R0 is None else R0
    return np.asarray(gN, dtype=float) * R0 / np.asarray(Rv, dtype=float)


# --------------------------------------------------------------------------
# 3. Courbes de rotation : Newton, TCE, MOND, RG + NFW
# --------------------------------------------------------------------------
def mond_g(gN, a0: float = ce.A0_MOND, interpolation: str = "simple"):
    """MOND : g = nu(gN/a0) gN, interpolations « simple » et « standard »."""
    y = np.asarray(gN, dtype=float) / a0
    if interpolation == "simple":
        nu = 0.5 + np.sqrt(0.25 + 1.0 / y)
    elif interpolation == "standard":
        nu = np.sqrt(0.5 + np.sqrt(0.25 + 1.0 / y**2))
    else:
        raise ValueError("interpolation inconnue : 'simple' ou 'standard'")
    return nu * np.asarray(gN, dtype=float)


def nfw_velocity(r, M200, c200: float = 10.0, H0: float = 70.0):
    """Vitesse circulaire d'un halo NFW (RG/Newton + matière noire).
    M200 en kg, H0 en km/s/Mpc."""
    H = H0 * 1e3 / MPC
    rho_c = 3.0 * H**2 / (8.0 * np.pi * G)
    r200 = (3.0 * M200 / (4.0 * np.pi * 200.0 * rho_c)) ** (1.0 / 3.0)
    rs = r200 / c200
    m = lambda x: np.log1p(x) - x / (1.0 + x)  # noqa: E731
    Ms = M200 / m(c200)
    x = np.asarray(r, dtype=float) / rs
    return np.sqrt(G * Ms * m(x) / np.asarray(r, dtype=float))


def compare_rotation_curves(r, M_bar, rc, M200=None, a0: float = ce.A0_MOND):
    """Quatre modèles sur le même profil de rayons (masse baryonique ponctuelle).
    M200 : masse du halo NFW (défaut 10 x M_bar)."""
    M200 = 10.0 * M_bar if M200 is None else M200
    r = np.asarray(r, dtype=float)
    gN = ce.g_newton(M_bar, r)
    Rv = ce.stiffness_profile(r, rc, gN, a0)
    v_n = ce.circular_velocity(r, gN)
    v_tce = ce.circular_velocity(r, g_eff_sqrt(gN, Rv))
    v_tceG = ce.circular_velocity(r, g_eff_varying_G(gN, Rv))
    v_mond = ce.circular_velocity(r, mond_g(gN, a0))
    v_nfw = np.sqrt(v_n**2 + nfw_velocity(r, M200) ** 2)
    return {
        "newton": v_n,
        "tce_sqrt": v_tce,
        "tce_varying_G": v_tceG,
        "mond": v_mond,
        "gr_nfw": v_nfw,
    }


# --------------------------------------------------------------------------
# 4. Échelle d'accélération et cosmologie
# --------------------------------------------------------------------------
def a0_cosmic_candidates(H0: float = 70.0):
    """Coïncidence a0 ~ c H0 : Milgrom c H0/(2 pi), Verlinde c H0/6. H0 en km/s/Mpc."""
    H = H0 * 1e3 / MPC
    return {
        "c*H0": c * H,
        "Milgrom c*H0/(2 pi)": c * H / (2.0 * np.pi),
        "Verlinde c*H0/6": c * H / 6.0,
        "a0 mesuré (MOND)": ce.A0_MOND,
    }


# --------------------------------------------------------------------------
# 5. Pont avec la thermodynamique de l'espace-temps (Jacobson)
# --------------------------------------------------------------------------
def planck_length_squared():
    return HBAR * G / c**3


def bekenstein_hawking_eta():
    """Densité d'entropie par unité d'aire de l'horizon : eta = 1 / (4 l_P^2) (k_B = 1)."""
    return 1.0 / (4.0 * planck_length_squared())


def jacobson_R0(eta: float | None = None):
    """R_0 = hbar c^2 eta / (2 pi)  (= c^5 / (8 pi G) pour eta = 1/(4 l_P^2),
    en puissance ; le manuscrit exprime R_0 en newtons : c^4/(8 pi G) = R_0/c)."""
    eta = bekenstein_hawking_eta() if eta is None else eta
    return HBAR * c**2 * eta / (2.0 * np.pi)
