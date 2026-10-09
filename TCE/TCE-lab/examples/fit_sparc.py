"""Ajuste Newton, MOND et Chrono-Energetics sur SPARC (175 galaxies).

Usage (depuis TCE-lab) :
    python examples/fit_sparc.py [--data data/sparc] [--quality 2] [--plots]
Télécharge SPARC au premier lancement (réseau requis), puis écrit
results/sparc_fit.csv et, avec --plots, results/sparc_fit.png.
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import sparc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="data/sparc")
ap.add_argument("--quality", type=int, default=2, help="qualité max (1 haute, 2 moyenne, 3 basse)")
ap.add_argument("--plots", action="store_true")
args = ap.parse_args()

sparc.download_sparc(args.data)
gals = sparc.load_sparc(args.data)
print(f"{len(gals)} galaxies chargées ; qualité <= {args.quality}")

res = {m: sparc.fit_all(gals, m, quality_max=args.quality) for m in ("newton", "mond", "tce")}
k_best, _ = sparc.fit_k_global(gals, quality_max=args.quality)
res["tce_k"] = sparc.fit_all(gals, "tce_k", k_rdisk=k_best, quality_max=args.quality)
nglob = {"newton": 0, "mond": 0, "tce": 0, "tce_k": 1}

print(f"\n{'modèle':8} {'N gal':>6} {'chi2/dof méd.':>14} {'chi2/dof tot.':>14} {'BIC':>10} {'rms (dex)':>10}")
for m, r in res.items():
    s = sparc.summarize(r, nglob[m])
    print(f"{m:8} {s['galaxies']:6d} {s['chi2_red_median']:14.2f} {s['chi2_red_total']:14.2f} "
          f"{s['BIC']:10.0f} {s['rms_dex']:10.3f}")
for m in ("mond", "tce"):
    gb_ = np.concatenate([x["g_bar"] for x in res[m]]); gm_ = np.concatenate([x["g_mod"] for x in res[m]])
    print(f"{m}: fraction de points avec g_modèle < g_baryonique : {(gm_ < 0.99 * gb_).mean():.1%}")
print(f"\nk global (rc = k * Rdisk) : {k_best:.2f}")

rc = np.array([r["rc"] for r in res["tce"]])
rdisk = np.array([g.rdisk for g in gals if g.name in {r['name'] for r in res['tce']}])
names = [r["name"] for r in res["tce"]]
rdisk = np.array([next(g.rdisk for g in gals if g.name == n) for n in names])
ok = np.isfinite(rdisk) & (rdisk > 0) & (rc > 0.11) & (rc < 300)
if ok.sum() > 3:
    corr = np.corrcoef(np.log10(rc[ok]), np.log10(rdisk[ok]))[0, 1]
    print(f"rc libre vs Rdisk : corrélation log-log = {corr:.2f} sur {ok.sum()} galaxies (hors bords de grille)")
print(f"galaxies avec rc au bord de la grille : {(~ok).sum()} / {len(rc)}")

os.makedirs("results", exist_ok=True)
with open("results/sparc_fit.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["galaxie", "modele", "chi2", "dof", "upsilon_disque", "rc_kpc"])
    for m, r in res.items():
        for x in r:
            w.writerow([x["name"], m, f"{x['chi2']:.2f}", x["dof"], f"{x['upsilon']:.3f}", f"{x['rc']:.3f}"])
print("Écrit results/sparc_fit.csv")

if args.plots:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for m, c_ in (("newton", "gray"), ("mond", "C0"), ("tce", "C1")):
        a = np.concatenate([x["g_obs"] for x in res[m]]) * 1e6 / sparc.KPC_M
        b = np.concatenate([x["g_bar"] for x in res[m]]) * 1e6 / sparc.KPC_M
        mod = np.concatenate([x["g_mod"] for x in res[m]]) * 1e6 / sparc.KPC_M
        ok_ = (a > 0) & (b > 0) & (mod > 0)
        if m == "newton":
            ax[0].loglog(b[ok_], a[ok_], ".", ms=2, color="lightgray", label="données")
        ax[0].loglog(b[ok_][::7], mod[ok_][::7], ".", ms=2, color=c_, label=m)
    xx = np.logspace(-13, -8, 50)
    ax[0].loglog(xx, xx, "k:", lw=0.8)
    ax[0].set(xlabel="g_bar (m/s²)", ylabel="g (m/s²)", title="Relation d'accélération radiale"); ax[0].legend(markerscale=4)
    for m, c_ in (("mond", "C0"), ("tce", "C1"), ("tce_k", "C2")):
        ax[1].hist([x["chi2"] / x["dof"] for x in res[m]], bins=np.logspace(-1, 2, 30), alpha=0.5, color=c_, label=m)
    ax[1].set(xscale="log", xlabel="chi2 / dof (par galaxie)", title="Qualité d'ajustement"); ax[1].legend()
    ax[2].loglog(rdisk[ok], rc[ok], "o", ms=4); ax[2].set(xlabel="Rdisk (kpc)", ylabel="rc ajusté (kpc)", title="rc libre vs Rdisk")
    fig.tight_layout(); fig.savefig("results/sparc_fit.png", dpi=130)
    print("Écrit results/sparc_fit.png")
