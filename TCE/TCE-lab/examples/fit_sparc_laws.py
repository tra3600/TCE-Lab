"""Cherche une loi qui fixe rc à partir de la densité de surface baryonique (SPARC).

Lois testées (paramètres globaux, valables pour toutes les galaxies) :
  L1  rc = k Rdisk
  L2  rc = k sqrt(G M_bar / a0)             (rayon où g_bar = a0)
  L3  rc = k Rdisk (Sigma_b / Sigma_dagger)^alpha, Sigma_dagger = a0/G
Chaque loi est ajustée sur la moitié des galaxies et évaluée sur l'autre moitié
(validation croisée à 2 plis, 20 tirages), puis comparée à MOND (a0 fixé).

Usage :  python examples/fit_sparc_laws.py [--base tce_v2] [--marginalize]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import sparc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="data/sparc")
ap.add_argument("--base", default="tce_v2", choices=["tce", "tce_v1", "tce_v2"])
ap.add_argument("--marginalize", action="store_true")
ap.add_argument("--seeds", type=int, default=20)
args = ap.parse_args()

sparc.download_sparc(args.data)
gals = sparc.load_sparc(args.data)
names_all = None
tables, params_all = {}, {}
for law in sparc.LAWS:
    params, names, table = sparc.law_chi2_table(gals, args.base, law, marginalize=args.marginalize)
    tables[law], params_all[law] = table, params
    names_all = names
    print(f"{law}: {len(params)} jeux de paramètres x {len(names)} galaxies")

by_name = {g.name: g for g in gals}
sel = [by_name[n] for n in names_all]
mond = {g.name: sparc.fit_galaxy(g, "mond", marginalize=args.marginalize) for g in sel}
mond_a0 = {g.name: sparc.fit_galaxy(g, "mond_a0", marginalize=args.marginalize) for g in sel}
free = {g.name: sparc.fit_galaxy(g, args.base, marginalize=args.marginalize) for g in sel}
tot = lambda d: sum(v["chi2"] for v in d.values())  # noqa: E731
dof = lambda d: sum(v["dof"] for v in d.values())   # noqa: E731
npts = sum(v["n"] for v in mond.values())

print(f"\nÉchantillon : {len(sel)} galaxies, {npts} points ; base = {args.base} ; "
      f"distance/inclinaison marginalisées : {args.marginalize}\n")
print(f"{'modèle':34} {'params glob.':>12} {'chi2 total':>11} {'chi2/dof':>9} {'hors éch.':>10}")
print(f"{'MOND (a0 fixé)':34} {0:12d} {tot(mond):11.0f} {tot(mond)/dof(mond):9.2f} {'-':>10}")
print(f"{'MOND, a0 libre / galaxie':34} {'0 (+1/gal.)':>12} {tot(mond_a0):11.0f} {tot(mond_a0)/dof(mond_a0):9.2f} {'-':>10}")
print(f"{args.base + ', rc libre / galaxie':34} {'0 (+1/gal.)':>12} {tot(free):11.0f} {tot(free)/dof(free):9.2f} {'-':>10}")

summary = {}
for law, table in tables.items():
    prm, i, chi2_in = sparc.law_best(params_all[law], table)
    cv = np.mean([sparc.cross_validate_law(table, seed=s) for s in range(args.seeds)])
    nparam = len(params_all[law][0])
    summary[law] = (prm, chi2_in, cv, nparam)
    lab = f"{args.base} + {law}"
    prm_s = ", ".join(f"{k}={v:.2f}" for k, v in prm.items())
    print(f"{lab:34} {nparam:12d} {chi2_in:11.0f} {chi2_in/dof(mond):9.2f} {cv:10.0f}   [{prm_s}]")

cv_mond = np.mean([  # MOND n'a aucun paramètre global : son chi2 hors échantillon = son chi2
    tot(mond)])
print(f"\nChi2 hors échantillon de MOND (0 paramètre global) : {cv_mond:.0f}")
best = min(summary, key=lambda k: summary[k][2])
print(f"Meilleure loi hors échantillon : {best}  (chi2 hors éch. = {summary[best][2]:.0f})")
print("Une loi n'est intéressante que si son chi2 hors échantillon bat celui de MOND (a0 fixé) "
      "et se rapproche du rc libre.")

# L3 : alpha significativement différent de zéro ?
p3, t3 = params_all["L3_surface"], tables["L3_surface"]
alphas = np.array([p["alpha"] for p in p3]); ks = np.array([p["k"] for p in p3])
best_by_alpha = {a: t3[alphas == a].sum(axis=1).min() for a in np.unique(alphas)}
a_best = min(best_by_alpha, key=best_by_alpha.get)
print(f"\nL3 : alpha optimal = {a_best:.2f} ; delta chi2 (alpha=0 -> optimal) = "
      f"{best_by_alpha.get(0.0, np.nan) - best_by_alpha[a_best]:.0f}")

# ---- À quoi rc libre est-il corrélé ? (galaxies dont rc n'est pas au bord de la grille)
rows = []
for g in sel:
    rc_ = free[g.name]["rc"]
    if not (0.15 < rc_ < 250) or not all(np.isfinite([g.rdisk, g.reff, g.sbdisk, g.lum, g.mhi, g.vflat])) \
            or g.vflat <= 0 or g.sbdisk <= 0 or g.rdisk <= 0:
        continue
    mb = sparc.baryonic_mass(g)
    rows.append([rc_, g.rdisk, g.reff, sparc.central_surface_density(g), mb, g.vflat,
                 sparc.HELIUM_FACTOR * g.mhi * 1e9 / mb, free[g.name]["upsilon"]])
if len(rows) > 10:
    A = np.log10(np.array(rows))
    labels = ["Rdisk", "Reff", "Sigma_b", "M_bar", "Vflat", "f_gaz", "Upsilon ajusté"]
    print(f"\nCorrélation de log rc avec les propriétés ({len(rows)} galaxies, rc hors bords) :")
    for j, lab in enumerate(labels, start=1):
        print(f"  {lab:16} r = {np.corrcoef(A[:, 0], A[:, j])[0, 1]:+.2f}")
    X = np.column_stack([np.ones(len(A)), A[:, [1, 3, 4]]])  # Rdisk, Sigma_b, M_bar
    coef, *_ = np.linalg.lstsq(X, A[:, 0], rcond=None)
    resid = A[:, 0] - X @ coef
    print(f"  Régression log rc ~ Rdisk + Sigma_b + M_bar : R^2 = {1 - resid.var() / A[:, 0].var():.2f}, "
          f"dispersion résiduelle = {resid.std():.2f} dex (rc couvre {np.ptp(A[:, 0]):.1f} dex)")
