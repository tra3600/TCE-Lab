"""Démonstration pas à pas des équations (1)-(8) de Chrono-Energetics.

Lancer depuis TCE-lab :  python examples/chrono_energetics_demo.py [--plots]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import chrono_energetics as ce  # noqa: E402

KPC = 3.0857e19  # m
MSUN = 1.989e30  # kg


def titre(txt):
    print("\n" + "=" * 72 + f"\n{txt}\n" + "=" * 72)


titre("Étape 1 - Éq. (2) : rigidité du vide R_0 = c^4 / (8 pi G)")
chk = ce.check_r0_value()
print(f"R_0 calculé          = {chk['R0_calcule_N']:.4e} N")
print(f"R_0 écrit manuscrit  = {chk['R0_manuscrit_N']:.4e} N  (rapport {chk['rapport_manuscrit_sur_calcule']:.2f})")
print(f"c^4/G (force Planck) = {chk['c4_sur_G_N']:.4e} N  -> le manuscrit semble confondre ces valeurs")

titre("Étape 2 - Éq. (1) et (3) : temps émergent et Chronom")
for dE in (1e-19, 1.0, 1e9):
    chi = ce.chronom_from_energy(dE)
    print(f"dE = {dE:9.2e} J -> {float(chi):.3e} chronom (m) -> dt = {float(ce.chronom_to_seconds(chi)):.3e} s")

titre("Étape 3 - Éq. (4) : dt = (8 pi G / c^5) dE")
print(f"coefficient 8 pi G / c^5 = {8*np.pi*ce.G/ce.c**5:.4e} s/J")
print(f"éq.(1) avec alpha=1/c == éq.(4) : {np.isclose(ce.emergent_time_dt(1.0), ce.dt_seconds(1.0))}")

titre("Étape 4 - Éq. (5)-(7) : courbe de rotation (galaxie 1e11 Msun, rc = 3 kpc)")
M = 1e11 * MSUN
r = np.array([0.5, 1, 3, 10, 30, 100]) * KPC
out = ce.rotation_curve(r, M, rc=3 * KPC)
print(f"{'r (kpc)':>8} {'g_N':>10} {'R_v/R_0':>10} {'g_eff/g_N':>10} {'v (km/s)':>9}")
for i, ri in enumerate(r):
    print(f"{ri/KPC:8.1f} {out['g_N'][i]:10.2e} {out['R_v'][i]/ce.vacuum_stiffness_R0():10.3e} "
          f"{out['g_eff'][i]/out['g_N'][i]:10.3e} {out['v'][i]/1e3:9.1f}")
print(f"Éq.(7) v = (G M a0)^(1/4) = {ce.tully_fisher_velocity(M)/1e3:.1f} km/s (plateau attendu)")

titre("Étape 5 - Éq. (8) : hystérésis du vide, jouet 1D de l'amas de la Balle")
print("Décalage final (minimum de R_v) - (pic de gaz), unités arbitraires :")
for adv in (False, True):
    for tau in (0.05, 0.3, 1.0):
        res = ce.simulate_bullet_cluster(tau_v=tau, advect_with_inertia=adv)
        print(f"  transport inertiel={adv!s:5}  tau_v={tau:4.2f}  décalage={res.offset[-1]:+.3f}")
print("-> L'éq.(8) seule donne un décalage NÉGATIF (la lentille retarde sur le gaz).")
print("   Le décalage positif de l'amas de la Balle exige un terme de transport (extension).")

if "--plots" in sys.argv:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.join(os.path.dirname(__file__), "figures"), exist_ok=True)
    rr = np.logspace(-0.5, 2.5, 200) * KPC
    o = ce.rotation_curve(rr, M, rc=3 * KPC)
    newton = ce.circular_velocity(rr, o["g_N"])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].semilogx(rr / KPC, newton / 1e3, "--", label="newtonien")
    ax[0].semilogx(rr / KPC, o["v"] / 1e3, label="Chrono-Energetics")
    ax[0].axhline(ce.tully_fisher_velocity(M) / 1e3, color="gray", ls=":", label="éq. (7)")
    ax[0].set(xlabel="r (kpc)", ylabel="v (km/s)", title="Courbe de rotation"); ax[0].legend()
    for adv, ls in ((False, "--"), (True, "-")):
        res = ce.simulate_bullet_cluster(tau_v=0.3, advect_with_inertia=adv)
        ax[1].plot(res.t, res.offset, ls, label=f"transport inertiel = {adv}")
    ax[1].axhline(0, color="gray", lw=0.8)
    ax[1].set(xlabel="t (u.a.)", ylabel="x_vide - x_gaz", title="Éq. (8) : hystérésis"); ax[1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(os.path.dirname(__file__), "figures", "chrono_energetics.png"), dpi=130)
    print("\nFigure : examples/figures/chrono_energetics.png")
