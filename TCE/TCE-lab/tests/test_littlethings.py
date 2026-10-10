import os
import tempfile
import unittest

import numpy as np

from tce import littlethings as lt
from tce import sparc

T1 = "TESTGAL |12 38 39.2 +32 45 41.0| 3.6| 306.2|1.3| 48.4|13.9|66.5| 5.2|-12.4|7.3|0.06|-2.58|0.01|-2.48|0.01\n"
T2 = ("TESTGAL |2.59|2.27| 23.5| 23.7| 13.5|0.24|-0.4|      |11.0| 478.5|       | 16.1| 1.6|2.01| 0.52|   8.19|"
      "  1.62|-1.25| 0.21| |+0.03| 0.27| |  2.91| 0.41| 0.49|8.529|  9.138\n")


class TestLittleThings(unittest.TestCase):
    def test_parse_tables(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "table1.dat"), "w").write(T1)
            open(os.path.join(d, "table2.dat"), "w").write(T2)
            g = lt.load(d)[0]
        self.assertEqual(g.name, "TESTGAL")
        self.assertAlmostEqual(g.distance, 3.6)
        self.assertAlmostEqual(g.v_rmax, 23.5)
        self.assertAlmostEqual(g.m_gas, 2.91e7)
        self.assertAlmostEqual(g.m_star, 0.49e7)          # SED prioritaire
        self.assertAlmostEqual(g.m_bar, 3.40e7)

    def test_stellar_mass_falls_back_to_kinematic(self):
        g = lt.LTGalaxy("x", 1, 1, 10, 10, 1e7, np.nan, 2e6)
        self.assertEqual(g.m_star, 2e6)

    def test_name_normalisation(self):
        for a, b in (("DDO_154", "DDO154"), ("UGC08508", "UGC8508"), ("NGC 2366", "NGC2366")):
            self.assertEqual(lt.normalize_name(a), lt.normalize_name(b))
        self.assertNotEqual(lt.normalize_name("DDO154"), lt.normalize_name("DDO1540"))

    def test_btfr_is_eq7(self):
        m = np.array([1e8, 1e9, 1e10])
        v = lt.btfr_velocity(m)
        np.testing.assert_allclose(v**4, sparc.G_KPC * m * sparc.A0_KMS2_KPC, rtol=1e-12)
        np.testing.assert_allclose(lt.btfr_residuals(v, m), 0.0, atol=1e-12)

    def test_slope_recovers_four_on_exact_relation(self):
        m = np.logspace(7, 11, 30)
        v = lt.btfr_velocity(m)
        s, _, scat = lt.fit_btfr_slope(v, m)
        self.assertAlmostEqual(s, 4.0, places=8)
        self.assertLess(scat, 1e-8)


class TestExternalLoader(unittest.TestCase):
    def test_load_rotmod_directory_round_trip(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "G1_rotmod.dat"), "w") as f:
                f.write("# Distance = 5.00 Mpc\n# Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n"
                        "1.0 20.0 2.0 5.0 10.0 0.0 100.0 0.0\n2.0 30.0 2.0 8.0 15.0 0.0 80.0 0.0\n")
            gal = sparc.load_rotmod_directory(d)[0]
        self.assertEqual(gal.name, "G1")
        self.assertAlmostEqual(gal.distance, 5.0)
        np.testing.assert_allclose(gal.vobs, [20.0, 30.0])
        np.testing.assert_allclose(gal.vdisk, [10.0, 15.0])


if __name__ == "__main__":
    unittest.main()
