"""Ajuste MOND et TCE (V2) sur un échantillon indépendant au format SPARC « *_rotmod.dat ».

Pour THINGS (de Blok+2008), il faut convertir les courbes baryoniques radiales (gaz, étoiles)
au format :  Rad[kpc]  Vobs  errV  Vgas  Vdisk(Upsilon=1)  Vbul  SBdisk  SBbul
avec l'en-tête « # Distance = X Mpc ». Les galaxies déjà dans SPARC sont exclues.

Usage :  python examples/fit_external.py --dir data/things_rotmod [--sparc-data data/sparc]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import littlethings as lt  # noqa: E402
from tce import sparc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="dossier de fichiers *_rotmod.dat")
ap.add_argument("--sparc-data", default="data/sparc")
args = ap.parse_args()

sparc.download_sparc(args.sparc_data)
sp_names = {lt.normalize_name(g.name) for g in sparc.load_sparc(args.sparc_data)}
gals = [g for g in sparc.load_rotmod_directory(args.dir) if lt.normalize_name(g.name) not in sp_names]
print(f"{len(gals)} galaxies indépendantes de SPARC dans {args.dir}")
if not gals:
    sys.exit("Aucun fichier *_rotmod.dat exploitable.")

models = ["newton", "mond", "mond_a0", "mond_rs", "mond_n", "tce_v1", "tce_v2"]
print(f"\n{'modèle':10} {'param./gal.':>12} {'chi2 total':>11} {'chi2/dof':>9} {'BIC':>9} {'rms (dex)':>10}")
for m in models:
    res = sparc.fit_all(gals, m, quality_max=3)
    if not res:
        continue
    s = sparc.summarize(res)
    nfree = 0 if m in ("newton", "mond") else 1
    print(f"{m:10} {nfree:12d} {s['chi2_total']:11.0f} {s['chi2_red_total']:9.2f} {s['BIC']:9.0f} {s['rms_dex']:10.3f}")
print("\nLecture : à nombre de paramètres égal, comparer tce_v2 à mond_rs, mond_a0 et mond_n ;")
print("les formes V1/V2 ayant été écrites après avoir vu SPARC, seul ce test hors échantillon les valide.")
