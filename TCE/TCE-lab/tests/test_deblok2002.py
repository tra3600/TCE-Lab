import os
import tempfile
import unittest

import numpy as np

from tce import deblok2002 as db

# colonnes à largeur fixe : nom (11) r (5) r_kpc (6) Vgas (7) Vdisk (7) Vrot (8) e_Vrot (7)
LINES = (
    "LSBCF563-1   2.5  0.55   1.02   8.33   15.52  22.84\n"
    "LSBCF563-1   7.5  1.64   3.92  18.05   46.42   9.46\n"
    "NGC1560      1.0  0.10          5.00   10.00   2.00\n"   # Vgas vide = 0
)


class TestDeBlok(unittest.TestCase):
    def test_clean_name(self):
        self.assertEqual(db.clean_name("LSBCF563-1"), "F563-1")
        self.assertEqual(db.clean_name("NGC1560     "), "NGC1560")

    def test_parse(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, db.FILE), "w").write(LINES)
            g = {x.name: x for x in db.load(d)}
        self.assertEqual(set(g), {"F563-1", "NGC1560"})
        np.testing.assert_allclose(g["F563-1"].r, [0.55, 1.64])
        np.testing.assert_allclose(g["F563-1"].vgas, [1.02, 3.92])
        np.testing.assert_allclose(g["F563-1"].vdisk, [8.33, 18.05])
        np.testing.assert_allclose(g["F563-1"].vobs, [15.52, 46.42])
        np.testing.assert_allclose(g["NGC1560"].vgas, [0.0])


if __name__ == "__main__":
    unittest.main()
