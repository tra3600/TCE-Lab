"""
Chrono-Energetics (G. Lombardo) : implémentation pas à pas des équations (1) à (8)
du manuscrit « Chrono-Energetics: An Emergent Approach to Spacetime Metric via
Local Dissipation and Vacuum Stiffness ».

Équations du manuscrit
----------------------
(1)  dt = alpha * dE_diss / R_v                        [section 2]
(2)  R_0 = c^4 / (8 pi G)                              [section 2]
(3)  [dE_diss / R_v] = J / N = m   (le « Chronom » chi)
(4)  dt_seconds = (8 pi G / c^5) * dE_diss             [alpha = 1/c]
(5)  R_v(r) = R_0 * (1 + (r/rc)^2) / (1 + (r/rc)^2 * a0/g_N(r))   [lecture B]
(6)  g_eff(r) = sqrt(g_N(r) * a0) = sqrt(G M_bar a0) / r
(7)  v^4 = G M_bar a0                                  (Tully-Fisher baryonique)
(8)  dR_v/dt = -(1/tau_v) * (R_v - R_v,static(rho))    (hystérésis du vide)

Choix d'interprétation (non précisés dans le manuscrit, à valider par l'auteur)
-------------------------------------------------------------------------------
* Éq. (5) : le texte extrait du PDF est ambigu (le rapport g_N / a0 n'est pas
  lisible). Deux lectures ont été testées :
    A) g_N/a0 au dénominateur : R_v *augmente* quand g_N << a0, ce qui contredit
       « la résistance de l'espace s'affaisse » (section 3) et ne redonne pas
       Newton quand g_N >> a0 ;
    B) a0/g_N au dénominateur : R_v *diminue* quand g_N << a0 (R_v ~ R_0 g_N/a0),
       et R_v -> R_0 pour r << rc. C'est la lecture retenue.
* Le lien entre R_v et g_eff n'est pas écrit ; avec la lecture B on utilise
  g_eff = g_N * sqrt(R_0 / R_v) (le couplage gravitationnel varie comme 1/R_v),
  qui redonne Newton pour r << rc et l'éq. (6) pour r >> rc avec g_N << a0
  (voir g_eff_from_stiffness et les tests).
* a0 = 1.2e-10 m/s^2 (valeur MOND usuelle) et rc sont des paramètres libres.
* R_v,static(rho) de l'éq. (8) n'est pas donné ; une forme illustrative
  (vide « ramolli » par la densité) est fournie et doit être remplacée.
* Erreur numérique du manuscrit : c^4/(8 pi G) vaut 4.82e42 N et non
  1.21e43 N (voir R0_PAPER_QUOTED et check_r0_value).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .constants import c, G

PI = np.pi

#: Accélération caractéristique de MOND (m/s^2), valeur usuelle.
A0_MOND = 1.2e-10
#: Valeur de R_0 imprimée dans le manuscrit (N) -- voir check_r0_value.
R0_PAPER_QUOTED = 1.21e43
#: Couplage alpha de l'éq. (1) en unités SI, choisi pour obtenir des secondes.
ALPHA_SI = 1.0 / c


# --------------------------------------------------------------------------
# Étape 1 -- Éq. (2) : rigidité du vide R_0
# --------------------------------------------------------------------------
def vacuum_stiffness_R0() -> float:
    """R_0 = c^4 / (8 pi G), l'inverse de la constante d'Einstein kappa (en N)."""
    return c**4 / (8.0 * PI * G)


def check_r0_value() -> dict:
    """Compare la valeur calculée de R_0 à celle écrite dans le manuscrit."""
    r0 = vacuum_stiffness_R0()
    return {
        "R0_calcule_N": r0,
        "R0_manuscrit_N": R0_PAPER_QUOTED,
        "rapport_manuscrit_sur_calcule": R0_PAPER_QUOTED / r0,
        "c4_sur_G_N": c**4 / G,  # force de Planck, 1.21e44 N
    }


# --------------------------------------------------------------------------
# Étape 2 -- Éq. (1) et (3) : le temps émergent et le Chronom
# --------------------------------------------------------------------------
def chronom_from_energy(dE_diss, R_v=None):
    """Éq. (3) : dE_diss / R_v, en joules / newtons = mètres (1 chronom = 1 m)."""
    R = vacuum_stiffness_R0() if R_v is None else R_v
    return np.asarray(dE_diss, dtype=float) / R


def emergent_time_dt(dE_diss, R_v=None, alpha: float = ALPHA_SI):
    """Éq. (1) : dt = alpha * dE_diss / R_v (secondes si alpha = 1/c)."""
    return alpha * chronom_from_energy(dE_diss, R_v)


def chronom_to_seconds(chi):
    """Convertit des chronoms (mètres) en secondes : t = chi / c."""
    return np.asarray(chi, dtype=float) / c


# --------------------------------------------------------------------------
# Étape 3 -- Éq. (4) : forme finale en secondes
# --------------------------------------------------------------------------
def dt_seconds(dE_diss):
    """Éq. (4) : dt = (8 pi G / c^5) * dE_diss (secondes, avec dE_diss en J)."""
    return (8.0 * PI * G / c**5) * np.asarray(dE_diss, dtype=float)


# --------------------------------------------------------------------------
# Étape 4 -- Éq. (5) à (7) : courbes de rotation et équivalence MOND
# --------------------------------------------------------------------------
def g_newton(M_bar, r):
    """Accélération newtonienne g_N(r) = G M_bar / r^2 (masse ponctuelle)."""
    r = np.asarray(r, dtype=float)
    return G * M_bar / r**2


def stiffness_profile(r, rc, g_N, a0: float = A0_MOND, R0: float | None = None):
    """Éq. (5), lecture B : R_v(r) = R0 (1+x^2) / (1 + x^2 a0/g_N), x = r/rc."""
    R0 = vacuum_stiffness_R0() if R0 is None else R0
    x2 = (np.asarray(r, dtype=float) / rc) ** 2
    return R0 * (1.0 + x2) / (1.0 + x2 * a0 / np.asarray(g_N, dtype=float))


def _w(r, rc):
    x2 = (np.asarray(r, dtype=float) / rc) ** 2
    return x2 / (1.0 + x2)


def stiffness_profile_v1(r, rc, g_N, a0: float = A0_MOND, R0: float | None = None):
    """Variante V1 de l'éq. (5) : R_v = R_0 / (1 + w a0/g_N), w = x^2/(1+x^2).
    Avec g_eff = g_N sqrt(R_0/R_v) (couplage en racine carrée du manuscrit) :
    g_eff = g_N sqrt(1 + w a0/g_N) >= g_N, Newton si r << rc ou g_N >> a0,
    sqrt(g_N a0) si r >> rc et g_N << a0. R_v <= R_0 partout (le vide ne fait que s'affaisser)."""
    R0 = vacuum_stiffness_R0() if R0 is None else R0
    return R0 / (1.0 + _w(r, rc) * a0 / np.asarray(g_N, dtype=float))


def stiffness_profile_v2(r, rc, g_N, a0: float = A0_MOND, R0: float | None = None):
    """Variante V2 : R_v = R_0 / (1 + w sqrt(a0/g_N)).
    Avec le couplage cohérent avec la RG, G_eff = 1/R_v, c.-à-d. g_eff = g_N R_0/R_v :
    g_eff = g_N + w sqrt(g_N a0) >= g_N, mêmes limites que V1, sans la racine
    carrée ajoutée à la main : le plateau vient directement de R_v ~ R_0 sqrt(g_N/a0)."""
    R0 = vacuum_stiffness_R0() if R0 is None else R0
    return R0 / (1.0 + _w(r, rc) * np.sqrt(a0 / np.asarray(g_N, dtype=float)))


def g_eff_deep_mond(g_N, a0: float = A0_MOND):
    """Éq. (6), régime g_N << a0 : g_eff = sqrt(g_N a0)."""
    return np.sqrt(np.asarray(g_N, dtype=float) * a0)


def g_eff_from_stiffness(g_N, R_v, R0: float | None = None):
    """Lecture proposée : g_eff = g_N * sqrt(R_0 / R_v).
    Redonne g_N pour r << rc et sqrt(g_N a0) pour r >> rc avec g_N << a0.
    Limite connue : pour r >> rc ET g_N >> a0, g_eff ~ g_N / (r/rc) (gravité
    affaiblie), conséquence de la forme (5) telle que lue."""
    R0 = vacuum_stiffness_R0() if R0 is None else R0
    return np.asarray(g_N, dtype=float) * np.sqrt(R0 / np.asarray(R_v, dtype=float))


def circular_velocity(r, g_eff):
    """Vitesse circulaire : v^2 / r = g_eff  ->  v = sqrt(g_eff r)."""
    return np.sqrt(np.asarray(g_eff, dtype=float) * np.asarray(r, dtype=float))


def tully_fisher_velocity(M_bar, a0: float = A0_MOND):
    """Éq. (7) : v^4 = G M_bar a0  ->  v = (G M_bar a0)^(1/4)."""
    return (G * M_bar * a0) ** 0.25


def rotation_curve(r, M_bar, rc, a0: float = A0_MOND):
    """Chaîne complète éq. (5) -> (6) -> v(r) pour un profil de rayons r."""
    gN = g_newton(M_bar, r)
    Rv = stiffness_profile(r, rc, gN, a0)
    geff = g_eff_from_stiffness(gN, Rv)
    return {"g_N": gN, "R_v": Rv, "g_eff": geff, "v": circular_velocity(r, geff)}


# --------------------------------------------------------------------------
# Étape 5 -- Éq. (8) : hystérésis du vide (amas de la Balle)
# --------------------------------------------------------------------------
def relax_stiffness(R_v, R_static, dt: float, tau_v: float):
    """Pas exact de l'éq. (8) avec R_static constant sur [t, t+dt] :
    R(t+dt) = R_static + (R(t) - R_static) * exp(-dt / tau_v)."""
    R_v = np.asarray(R_v, dtype=float)
    R_static = np.asarray(R_static, dtype=float)
    return R_static + (R_v - R_static) * np.exp(-dt / tau_v)


def default_static_stiffness(rho, rho_ref: float, R0: float | None = None):
    """Forme ILLUSTRATIVE de R_v,static(rho) : R0 / (1 + rho/rho_ref)
    (le vide « ramollit » là où la densité est élevée). À remplacer."""
    R0 = vacuum_stiffness_R0() if R0 is None else R0
    return R0 / (1.0 + np.asarray(rho, dtype=float) / rho_ref)


@dataclass
class BulletResult:
    x: np.ndarray
    t: np.ndarray
    rho: np.ndarray          # (nt, nx) densité baryonique
    R_v: np.ndarray          # (nt, nx) rigidité dynamique
    R_static: np.ndarray     # (nt, nx) rigidité statique
    x_gas: np.ndarray        # position du pic de gaz
    x_vacuum: np.ndarray     # position du minimum de R_v (pic « de lentille »)

    @property
    def offset(self) -> np.ndarray:
        """Décalage entre le minimum de rigidité (lentille) et le gaz."""
        return self.x_vacuum - self.x_gas


def simulate_bullet_cluster(
    tau_v: float,
    nx: int = 400,
    nt: int = 400,
    L: float = 1.0,
    t_end: float = 1.0,
    v0: float = 1.0,
    decel: float = 0.8,
    width: float = 0.04,
    rho_peak: float = 1.0,
    rho_ref: float = 0.5,
    static_fn: Callable | None = None,
    advect_with_inertia: bool = False,
) -> BulletResult:
    """Jouet 1D de l'éq. (8). Un nuage de gaz (gaussien) est freiné par friction
    hydrodynamique, x_gas(t) = x0 + v0 t - 0.5 decel v0 t^2 / t_end ; la rigidité
    R_v relaxe vers R_static(rho) avec le temps caractéristique tau_v.
    Unités arbitraires (L, t_end, v0). R_v est normalisée par R0 = 1.

    L'éq. (8) telle qu'écrite est une relaxation purement locale (aucun terme de
    transport) : le minimum de R_v retarde alors sur le gaz au lieu de le dépasser.
    Avec advect_with_inertia=True (EXTENSION, absente du manuscrit) le champ R_v
    est en plus transporté à la vitesse inertielle initiale v0, ce qui traduit
    l'affirmation du texte (« le profil continue sur les trajectoires inertielles »)."""
    x = np.linspace(0.0, L, nx)
    t = np.linspace(0.0, t_end, nt)
    dt = t[1] - t[0]
    x0 = 0.15 * L
    static_fn = static_fn or (lambda rho: default_static_stiffness(rho, rho_ref, R0=1.0))

    rho_all = np.empty((nt, nx))
    Rs_all = np.empty((nt, nx))
    Rv_all = np.empty((nt, nx))
    x_gas = np.empty(nt)
    Rv = static_fn(np.zeros(nx))  # vide au repos
    for k, tk in enumerate(t):
        xg = x0 + v0 * tk - 0.5 * decel * v0 * tk**2 / t_end
        rho = rho_peak * np.exp(-0.5 * ((x - xg) / width) ** 2)
        Rs = static_fn(rho)
        if k > 0:
            if advect_with_inertia:
                Rv = np.interp(x - v0 * dt, x, Rv, left=static_fn(0.0), right=Rv[-1])
            Rv = relax_stiffness(Rv, Rs, dt, tau_v)
        rho_all[k], Rs_all[k], Rv_all[k], x_gas[k] = rho, Rs, Rv, xg
    x_vac = x[np.argmin(Rv_all, axis=1)]
    x_vac[0] = x_gas[0]  # à t=0 le champ est plat : argmin non défini
    return BulletResult(x, t, rho_all, Rv_all, Rs_all, x_gas, x_vac)
