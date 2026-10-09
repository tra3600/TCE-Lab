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

    def test_summary_keys(self):
        res = [sparc.fit_galaxy(synthetic_galaxy(), m) for m in ("mond", "tce")]
        s = sparc.summarize(res[:1])
        self.assertEqual(s["galaxies"], 1)
        self.assertGreater(s["points"], 10)


if __name__ == "__main__":
    unittest.main()
