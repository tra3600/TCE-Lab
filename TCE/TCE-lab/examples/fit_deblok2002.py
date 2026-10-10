"""Test hors échantillon sur de Blok & Bosma (2002) : galaxies communes avec SPARC exclues.

Compare Newton, MOND (a0 fixé), les contrôles MOND à un paramètre libre par galaxie et TCE V1/V2
sur des courbes baryoniques radiales qui ne servent pas à construire les formes V1/V2.
Le rapport masse/luminosité est en bande R : prior log-normale centrée sur --ups (défaut 1,0)
d'écart-type --sigma dex. Aucune marginalisation de D et i (non fournis).

Usage :  python examples/fit_deblok2002.py [--ups 1.0] [--sigma 0.2]
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import deblok2002 as db  # noqa: E402
from tce import littlethings as lt  # noqa: E402
from tce import sparc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="data/deblok2002")
ap.add_argument("--sparc-data", default="data/sparc")
ap.add_argument("--ups", type=float, default=1.0, help="centre de la prior sur M/L en bande R")
ap.add_argument("--sigma", type=float, default=0.2, help="écart-type de la prior (dex)")
args = ap.parse_args()

db.download(args.data)
sparc.download_sparc(args.sparc_data)
allg = db.load(args.data)
sp_names = {lt.normalize_name(g.name) for g in sparc.load_sparc(args.sparc_data)}
overlap = [g.name for g in allg if lt.normalize_name(g.name) in sp_names]
gals = [g for g in allg if lt.normalize_name(g.name) not in sp_names]
print(f"{len(allg)} galaxies ; communes avec SPARC (exclues) : {overlap}")

models = ["newton", "mond", "mond_a0", "mond_rs", "mond_n", "tce", "tce_v1", "tce_v2"]
with sparc.upsilon_prior(args.ups, args.sigma):
    res = {m: sparc.fit_all(gals, m, quality_max=3) for m in models}
names = [x["name"] for x in res["mond"]]
print(f"{len(names)} galaxies indépendantes avec >= {sparc.MIN_POINTS} points : {names}")
print(f"prior sur M/L(R) : {args.ups} ± {args.sigma} dex\n")
print(f"{'modèle':10} {'param./gal.':>12} {'chi2 total':>11} {'chi2/dof':>9} {'BIC':>9} {'rms (dex)':>10}")
for m in models:
    s = sparc.summarize(res[m])
    print(f"{m:10} {0 if m in ('newton', 'mond') else 1:12d} {s['chi2_total']:11.0f} {s['chi2_red_total']:9.2f} "
          f"{s['BIC']:9.0f} {s['rms_dex']:10.3f}")


def sign_test(a, b):
    d = np.array([x["chi2"] for x in res[a]]) - np.array([x["chi2"] for x in res[b]])
    wins, n = int((d < 0).sum()), int((d != 0).sum())
    tail = sum(math.comb(n, k) for k in range(0, min(wins, n - wins) + 1)) / 2**n
    return wins, n, min(1.0, 2 * tail), float(np.median(d))


print("\nTest de signe par galaxie (même nombre de paramètres libres) :")
for a, b in (("tce_v2", "mond_rs"), ("tce_v2", "mond_a0"), ("tce_v2", "mond_n"), ("tce_v1", "mond_rs")):
    w, n, p, med = sign_test(a, b)
    print(f"  {a:8} contre {b:8} : {a} meilleur sur {w}/{n} galaxies, p = {p:.3g}, médiane delta-chi2 = {med:+.2f}")
print("\nPar galaxie (chi2) : " + ", ".join(f"{n}: V2 {res['tce_v2'][i]['chi2']:.0f} / rs {res['mond_rs'][i]['chi2']:.0f}"
                                          for i, n in enumerate(names)))
