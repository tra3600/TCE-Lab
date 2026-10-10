"""Contrôle équitable : TCE contre des familles MOND à UN paramètre libre par galaxie.

Chaque modèle a les mêmes nuisances (Upsilon, et D et i avec --marginalize) et un seul
paramètre libre de plus par galaxie :
  mond_a0   : a0 libre
  mond_rs   : a0 fixé, rayon de transition rs libre (même structure que TCE V2)
  mond_n    : a0 fixé, indice d'interpolation n libre
  tce_v1, tce_v2, tce : rc libre
Comparaison par le chi2 total, le BIC et un test de signe galaxie par galaxie.

Usage :  python examples/fit_sparc_controls.py [--marginalize]
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import sparc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="data/sparc")
ap.add_argument("--marginalize", action="store_true")
args = ap.parse_args()

sparc.download_sparc(args.data)
gals = sparc.load_sparc(args.data)
models = ["mond", "mond_a0", "mond_rs", "mond_n", "tce", "tce_v1", "tce_v2"]
res = {m: sparc.fit_all(gals, m, marginalize=args.marginalize) for m in models}
names = [x["name"] for x in res["mond"]]
print(f"{len(names)} galaxies ; D et i marginalisés : {args.marginalize}\n")
print(f"{'modèle':10} {'param./gal.':>12} {'chi2 total':>11} {'chi2/dof':>9} {'BIC':>9} {'rms (dex)':>10}")
for m in models:
    s = sparc.summarize(res[m])
    print(f"{m:10} {'0' if m == 'mond' else '1':>12} {s['chi2_total']:11.0f} {s['chi2_red_total']:9.2f} "
          f"{s['BIC']:9.0f} {s['rms_dex']:10.3f}")


def sign_test(a, b):
    """Galaxies où a fait mieux que b ; p-value bilatérale exacte (binomiale 1/2)."""
    d = np.array([x["chi2"] for x in res[a]]) - np.array([x["chi2"] for x in res[b]])
    wins, n = int((d < 0).sum()), int((d != 0).sum())
    tail = sum(math.comb(n, k) for k in range(0, min(wins, n - wins) + 1)) / 2**n
    return wins, n, min(1.0, 2 * tail), float(np.median(d))


print("\nTest de signe par galaxie (modèle A contre B, même nombre de paramètres libres) :")
for a, b in (("tce_v2", "mond_rs"), ("tce_v2", "mond_a0"), ("tce_v2", "mond_n"), ("tce_v1", "mond_rs"),
             ("mond_rs", "mond_a0"), ("mond_n", "mond_a0")):
    w, n, p, med = sign_test(a, b)
    print(f"  {a:8} contre {b:8} : {a} meilleur sur {w}/{n} galaxies, p = {p:.3g}, médiane delta-chi2 = {med:+.2f}")

nu = np.array([x["rc"] for x in res["mond_n"]])
print(f"\nIndice d'interpolation ajusté n : médiane {np.median(nu):.2f} "
      f"[16-84 % : {np.percentile(nu, 16):.2f}-{np.percentile(nu, 84):.2f}] (n=1 simple, n=2 standard)")
