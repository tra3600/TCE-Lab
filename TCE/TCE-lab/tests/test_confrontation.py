import unittest

import numpy as np

from tce import chrono_energetics as ce
from tce import confrontation as cf
from tce.constants import c

M = 1e11 * cf.MSUN
RC = 3 * cf.KPC


class TestClocks(unittest.TestCase):
    def test_weak_field_matches_schwarzschild_far_away(self):
        r = np.array([1e9, 1e12]) * 1.0
        np.testing.assert_allclose(
            cf.gr_clock_rate_weak_field(r, cf.MSUN), cf.gr_clock_rate_schwarzschild(r, cf.MSUN), rtol=1e-9
        )

    def test_black_hole_profile_reproduces_gr(self):
        r = np.array([3.0, 10.0, 100.0]) * cf.schwarzschild_radius(cf.MSUN)
        Rv = cf.tce_gr_profile_for_black_hole(r, cf.MSUN)
        np.testing.assert_allclose(cf.tce_clock_rate(Rv), cf.gr_clock_rate_schwarzschild(r, cf.MSUN), rtol=1e-12)

    def test_tce_clock_deviation_is_orders_of_magnitude_above_gr(self):
        out = cf.clock_test_galaxy(np.array([30.0, 100.0]) * cf.KPC, M, RC)
        self.assertTrue(np.all(np.abs(out["ratio_of_deviations"]) > 1e5))


class TestCoupling(unittest.TestCase):
    def test_varying_G_gives_rising_curve_not_flat(self):
        r = np.array([30.0, 100.0]) * cf.KPC
        v = cf.compare_rotation_curves(r, M, RC)["tce_varying_G"]
        self.assertGreater(v[1] / v[0], 1.5)

    def test_sqrt_reading_is_flat(self):
        r = np.array([30.0, 100.0]) * cf.KPC
        v = cf.compare_rotation_curves(r, M, RC)["tce_sqrt"]
        self.assertAlmostEqual(v[1] / v[0], 1.0, delta=0.02)


class TestReferenceModels(unittest.TestCase):
    def test_mond_limits(self):
        for kind in ("simple", "standard"):
            self.assertAlmostEqual(float(cf.mond_g(1e-3, interpolation=kind)) / 1e-3, 1.0, places=3)
            deep = float(cf.mond_g(1e-14, interpolation=kind))
            self.assertAlmostEqual(deep / np.sqrt(1e-14 * ce.A0_MOND), 1.0, delta=0.01)

    def test_nfw_is_finite_and_plausible(self):
        v = cf.nfw_velocity(np.array([10.0, 50.0, 200.0]) * cf.KPC, 1e12 * cf.MSUN)
        self.assertTrue(np.all((v > 5e4) & (v < 4e5)))

    def test_all_models_agree_at_small_radius(self):
        out = cf.compare_rotation_curves(np.array([0.3]) * cf.KPC, M, RC, M200=1e-3 * M)
        for k in ("tce_sqrt", "mond", "gr_nfw"):
            self.assertAlmostEqual(out[k][0] / out["newton"][0], 1.0, delta=0.05)


class TestBridges(unittest.TestCase):
    def test_a0_cosmic_coincidence_within_15_percent(self):
        cand = cf.a0_cosmic_candidates()
        for k in ("Milgrom c*H0/(2 pi)", "Verlinde c*H0/6"):
            self.assertAlmostEqual(cand[k] / ce.A0_MOND, 1.0, delta=0.15)

    def test_jacobson_R0_equals_c_times_manuscript_R0(self):
        self.assertAlmostEqual(cf.jacobson_R0() / (c * ce.vacuum_stiffness_R0()), 1.0, places=9)


if __name__ == "__main__":
    unittest.main()
