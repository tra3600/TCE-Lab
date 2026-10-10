"""Test indépendant de l'éq. (7) du manuscrit, v^4 = G M_bar a0, sur LITTLE THINGS (26 naines).

SPARC sert de référence mais les galaxies communes aux deux échantillons sont exclues de
LITTLE THINGS. Seules des quantités intégrées sont disponibles pour LITTLE THINGS (voir
tce/littlethings.py) : ce test porte sur la normalisation (a0) et la pente de la relation
de Tully-Fisher baryonique, pas sur la forme des courbes de rotation.

Usage :  python examples/test_littlethings_btfr.py [--lt-data data/littlethings] [--sparc-data data/sparc]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from tce import littlethings as lt  # noqa: E402
from tce import sparc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--lt-data", default="data/littlethings")
ap.add_argument("--sparc-data", default="data/sparc")
ap.add_argument("--boot", type=int, default=5000)
args = ap.parse_args()

lt.download(args.lt_data)
sparc.download_sparc(args.sparc_data)
ltg = lt.load(args.lt_data)
sp = sparc.load_sparc(args.sparc_data)
sp_names = {lt.normalize_name(g.name) for g in sp}

overlap = [g.name for g in ltg if lt.normalize_name(g.name) in sp_names]
no_star = [g.name for g in ltg if not np.isfinite(g.m_star) and g.name not in overlap]
indep = [g for g in ltg if g.name not in overlap and np.isfinite(g.m_star) and g.v_rmax > 0]
print(f"LITTLE THINGS : {len(ltg)} galaxies ; communes avec SPARC (exclues) : {overlap}")
print(f"sans masse stellaire (exclues) : {no_star}")
print(f"échantillon indépendant : {len(indep)} galaxies\n")

rng = np.random.default_rng(0)


def report(label, v, m):
    r = lt.btfr_residuals(v, m)
    slope, icpt, scat = lt.fit_btfr_slope(v, m)
    idx = rng.integers(0, len(v), size=(args.boot, len(v)))
    slopes = np.array([np.polyfit(np.log10(v[i]), np.log10(m[i]), 1)[0] for i in idx])
    a0_fit = sparc.A0_KMS2_KPC * 10 ** (4 * r.mean())
    a0_si = a0_fit * 1e6 / sparc.KPC_M
    print(f"{label:34} N={len(v):3d}  <log(v_obs/v_pred)> = {r.mean():+.3f} ± {r.std(ddof=1)/np.sqrt(len(v)):.3f} "
          f"(rms {r.std(ddof=1):.3f})  a0 ajusté = {a0_si:.2e} m/s²")
    lo, hi = np.percentile(slopes, [16, 84])
    print(f"{'':34} pente de la BTFR = {slope:.2f} [{lo:.2f}-{hi:.2f}] (prédite : 4), dispersion en log M = {scat:.2f} dex")
    return r


v_s, m_s, _ = lt.sparc_btfr_sample(sp)
report("SPARC (Vflat, qualité 1)", v_s, m_s)

mb = np.array([g.m_bar for g in indep])
report("LITTLE THINGS, V(Rmax)", np.array([g.v_rmax for g in indep]), mb)
report("LITTLE THINGS, Viso(Rmax)", np.array([g.v_iso for g in indep]), mb)

# sous-échantillon dominé par le gaz (moins sensible au rapport masse/luminosité)
gas_dom = [g for g in indep if g.m_gas > 2 * g.m_star]
if len(gas_dom) >= 6:
    report("LITTLE THINGS, gaz dominant", np.array([g.v_rmax for g in gas_dom]), np.array([g.m_bar for g in gas_dom]))

# sensibilité à la masse stellaire : SED contre cinématique
both = [g for g in indep if np.isfinite(g.m_star_sed) and np.isfinite(g.m_star_kin)]
if len(both) >= 6:
    for lab, key in (("M* SED", "m_star_sed"), ("M* cinématique", "m_star_kin")):
        m_ = np.array([g.m_gas + getattr(g, key) for g in both])
        report(f"LITTLE THINGS ({lab}, {len(both)} gal.)", np.array([g.v_rmax for g in both]), m_)

print("\nAvertissement : V(Rmax) peut sous-estimer la vitesse plate si la courbe monte encore ; les")
print("masses stellaires (SED ou cinématique) et gazeuses (hélium inclus ?) ont leurs propres incertitudes.")

# influence : pente de la BTFR en retirant une galaxie à la fois (échantillon indépendant, V(Rmax))
v_all = np.array([g.v_rmax for g in indep])
slopes_loo = np.array([lt.fit_btfr_slope(np.delete(v_all, i), np.delete(mb, i))[0] for i in range(len(indep))])
worst = int(np.argmax(np.abs(slopes_loo - lt.fit_btfr_slope(v_all, mb)[0])))
print(f"\nInfluence : pente en retirant une galaxie = {slopes_loo.min():.2f} à {slopes_loo.max():.2f} ; "
      f"la plus influente est {indep[worst].name} (V = {indep[worst].v_rmax:.0f} km/s, "
      f"M_bar = {indep[worst].m_bar:.1e} M_sun)")
print(f"Étendue des vitesses : {v_all.min():.0f}-{v_all.max():.0f} km/s ({np.log10(v_all.max()/v_all.min()):.1f} dex)")
