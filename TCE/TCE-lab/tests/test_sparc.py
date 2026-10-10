import unittest

import numpy as np

from tce import sparc


def synthetic_galaxy(rc_true=4.0, upsilon=0.5):
    r = np.linspace(0.5, 30.0, 40)
    vdisk = 110.0 * np.sqrt(r / (r + 2.0)) * np.exp(-r / 40.0)   # à Upsilon = 1
    vgas = 40.0 * np.sqrt(r / (r + 5.0))
    gal = sparc.Galaxy("synth", r, np.ones_like(r), np.ones_like(r), vgas, vdisk, np.zeros_like(r))
    gbar = sparc.g_baryon(gal, upsilon)[0]
    gal.vobs = sparc.velocity(sparc.g_tce(gbar, r, rc_true), r)
    gal.err = np.full_like(r, 1.0)
    return gal


class TestSparcFit(unittest.TestCase):
    def test_models_limits(self):
        gb = np.array([1e8, 1.0])  # (km/s)^2/kpc : très au-dessus et très au-dessous de a0
        r = np.array([1.0, 1.0])
        np.testing.assert_allclose(sparc.g_mond_simple(gb, r)[0] / gb[0], 1.0, rtol=1e-2)
        np.testing.assert_allclose(sparc.g_mond_simple(gb, r)[1], np.sqrt(gb[1] * sparc.A0_KMS2_KPC), rtol=2e-2)
        # TCE : Newton pour rc >> r
        np.testing.assert_allclose(sparc.g_tce(gb, r, rc=1e6), gb, rtol=1e-6)

    def test_fit_recovers_rc_and_beats_other_models(self):
        gal = synthetic_galaxy(rc_true=4.0)
        tce = sparc.fit_galaxy(gal, "tce")
        self.assertAlmostEqual(np.log10(tce["rc"]), np.log10(4.0), delta=0.06)
        self.assertAlmostEqual(tce["upsilon"], 0.5, delta=0.06)
        for other in ("newton", "mond"):
            self.assertLess(tce["chi2"], sparc.fit_galaxy(gal, other)["chi2"])

    def test_variants_never_below_baryons_and_have_right_limits(self):
        gb = np.logspace(-2, 8, 60)             # (km/s)^2/kpc, de << a0 à >> a0
        for fn in (sparc.g_tce_v1, sparc.g_tce_v2):
            for r, rc in ((1.0, 1e-3), (1.0, 1.0), (1.0, 1e3), (100.0, 1.0)):
                self.assertTrue(np.all(fn(gb, r, rc) >= gb * (1 - 1e-12)))
            np.testing.assert_allclose(fn(gb, 1.0, 1e6), gb, rtol=1e-6)               # r << rc : Newton
            np.testing.assert_allclose(fn(np.array([1e9]), 1e3, 1.0), 1e9, rtol=5e-3)  # g >> a0 : Newton
            deep = fn(np.array([1e-2]), 1e3, 1.0)                                    # r >> rc, g << a0
            np.testing.assert_allclose(deep, np.sqrt(1e-2 * sparc.A0_KMS2_KPC), rtol=2e-2)

    def test_fit_recovers_rc_for_variant_v2(self):
        gal = synthetic_galaxy(rc_true=4.0)
        gbar = sparc.g_baryon(gal, 0.5)[0]
        gal.vobs = sparc.velocity(sparc.g_tce_v2(gbar, gal.r, 4.0), gal.r)
        fit = sparc.fit_galaxy(gal, "tce_v2")
        self.assertAlmostEqual(np.log10(fit["rc"]), np.log10(4.0), delta=0.06)

    def test_summary_keys(self):
        res = [sparc.fit_galaxy(synthetic_galaxy(), m) for m in ("mond", "tce")]
        s = sparc.summarize(res[:1])
        self.assertEqual(s["galaxies"], 1)
        self.assertGreater(s["points"], 10)


class TestMondControls(unittest.TestCase):
    def test_interpolation_index_recovers_simple_and_standard(self):
        gb = np.logspace(-1, 7, 30)
        np.testing.assert_allclose(sparc.g_mond_n(gb, n=1.0), sparc.g_mond_simple(gb, None), rtol=1e-9)
        y = gb / sparc.A0_KMS2_KPC
        standard = gb * np.sqrt(0.5 + np.sqrt(0.25 + 1.0 / y**2))
        np.testing.assert_allclose(sparc.g_mond_n(gb, n=2.0), standard, rtol=1e-9)

    def test_deep_limit_is_independent_of_n(self):
        gb = np.array([1e-8])
        for n in (0.5, 1.0, 3.0):
            np.testing.assert_allclose(sparc.g_mond_n(gb, n=n), np.sqrt(gb * sparc.A0_KMS2_KPC), rtol=2e-2)

    def test_gated_mond_limits(self):
        gb = np.logspace(0, 5, 20)
        np.testing.assert_allclose(sparc.g_mond_gated(gb, 1.0, 1e6), gb, rtol=1e-6)               # rs >> r : Newton
        np.testing.assert_allclose(sparc.g_mond_gated(gb, 1.0, 1e-6), sparc.g_mond_simple(gb, 1.0), rtol=1e-6)  # rs << r : MOND

    def test_gated_mond_fit_recovers_rs(self):
        gal = synthetic_galaxy()
        gbar = sparc.g_baryon(gal, 0.5)[0]
        gal.vobs = sparc.velocity(sparc.g_mond_gated(gbar, gal.r, 5.0), gal.r)
        fit = sparc.fit_galaxy(gal, "mond_rs")
        self.assertAlmostEqual(np.log10(fit["rc"]), np.log10(5.0), delta=0.08)
        self.assertEqual(fit["dof"], fit["n"] - 2)


class TestNuisance(unittest.TestCase):
    def _galaxy_with_distance_and_inclination_bias(self, fD=1.1, inc_true=65.0):
        gal = synthetic_galaxy(rc_true=4.0)
        gal.distance, gal.e_dist, gal.inc, gal.e_inc = 10.0, 1.0, 60.0, 5.0
        sq = np.sqrt(fD)
        scaled = sparc.Galaxy("s", gal.r * fD, gal.vobs, gal.err, gal.vgas * sq, gal.vdisk * sq, gal.vbul * sq)
        v = sparc.velocity(sparc.g_tce_v2(sparc.g_baryon(scaled, 0.5)[0], scaled.r, 4.0 * fD), scaled.r)
        ratio = np.sin(np.radians(60.0)) / np.sin(np.radians(inc_true))
        gal.vobs = v / ratio
        gal.err = np.full_like(gal.r, 1.0)
        return gal

    def test_marginalization_absorbs_distance_and_inclination_bias(self):
        gal = self._galaxy_with_distance_and_inclination_bias()
        fixed = sparc.fit_galaxy(gal, "tce_v2")
        marg = sparc.fit_galaxy(gal, "tce_v2", marginalize=True)
        self.assertLess(marg["chi2"], 0.3 * fixed["chi2"])
        self.assertAlmostEqual(marg["dist_factor"], 1.1, delta=0.1)

    def test_marginalization_never_worse_than_fixed_up_to_priors(self):
        gal = synthetic_galaxy()
        gal.distance, gal.e_dist, gal.inc, gal.e_inc = 10.0, 1.0, 60.0, 5.0
        a = sparc.fit_galaxy(gal, "mond")["chi2"]
        b = sparc.fit_galaxy(gal, "mond", marginalize=True)["chi2"]
        self.assertLessEqual(b, a + 1e-9)  # le point neutre (z = 0) est dans la grille

    def test_missing_metadata_falls_back_to_no_nuisance(self):
        gal = synthetic_galaxy()
        self.assertEqual(len(sparc._nuisance_grid(gal)), 1)


class TestRcLaws(unittest.TestCase):
    def test_alpha_zero_recovers_rdisk_law(self):
        gal = synthetic_galaxy()
        gal.rdisk, gal.sbdisk = 2.0, 500.0
        self.assertAlmostEqual(sparc.rc_law_surface(gal, 1.7, 0.0), sparc.rc_law_rdisk(gal, 1.7), places=12)

    def test_surface_law_scaling(self):
        gal = synthetic_galaxy()
        gal.rdisk, gal.sbdisk = 2.0, sparc.SIGMA_DAGGER / sparc.UPSILON_STAR   # Sigma_b = 1 Sigma_dagger
        self.assertAlmostEqual(sparc.rc_law_surface(gal, 1.0, 1.0) / (gal.rdisk * 1.0), 1.0, places=9)

    def test_cross_validation_picks_global_minimum_on_toy_table(self):
        # 3 jeux de paramètres, 10 galaxies : le jeu 1 est le meilleur partout
        table = np.array([[5.0] * 10, [1.0] * 10, [3.0] * 10])
        self.assertAlmostEqual(sparc.cross_validate_law(table), 10.0)

    def test_dagger_surface_density_is_physical(self):
        self.assertTrue(500 < sparc.SIGMA_DAGGER < 1500)


if __name__ == "__main__":
    unittest.main()
