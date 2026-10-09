import math
import unittest

import numpy as np

from tce import chrono_energetics as ce
from tce.constants import c, G


class TestEquations(unittest.TestCase):
    def test_eq2_R0(self):
        self.assertAlmostEqual(ce.vacuum_stiffness_R0() / 4.8154538867e42, 1.0, places=8)

    def test_manuscript_R0_value_is_inconsistent(self):
        # le manuscrit imprime 1.21e43 N ; c^4/(8 pi G) = 4.82e42 N
        ratio = ce.check_r0_value()["rapport_manuscrit_sur_calcule"]
        self.assertGreater(abs(ratio - 1.0), 1.0)

    def test_eq1_equals_eq4_with_alpha_one_over_c(self):
        dE = np.array([1e-6, 1.0, 1e9])
        np.testing.assert_allclose(ce.emergent_time_dt(dE), ce.dt_seconds(dE), rtol=1e-12)

    def test_eq3_chronom_is_length(self):
        # 1 J / 1 N = 1 m ; et dt = chi / c
        chi = ce.chronom_from_energy(ce.vacuum_stiffness_R0())
        self.assertAlmostEqual(float(chi), 1.0)
        self.assertAlmostEqual(float(ce.chronom_to_seconds(chi)), 1.0 / c)

    def test_planck_energy_gives_8pi_planck_time(self):
        hbar = 1.054571817e-34
        E_P = math.sqrt(hbar * c**5 / G)
        t_P = math.sqrt(hbar * G / c**5)
        self.assertAlmostEqual(float(ce.dt_seconds(E_P)) / (8 * math.pi * t_P), 1.0, places=9)


class TestRotation(unittest.TestCase):
    M = 1e41  # kg

    def test_newtonian_limit_inside_rc(self):
        r = np.array([1e17])
        out = ce.rotation_curve(r, self.M, rc=1e19)
        np.testing.assert_allclose(out["g_eff"], out["g_N"], rtol=1e-3)

    def test_deep_mond_limit_and_tully_fisher(self):
        r = np.array([1e22])
        out = ce.rotation_curve(r, self.M, rc=3e19)
        np.testing.assert_allclose(out["g_eff"], ce.g_eff_deep_mond(out["g_N"]), rtol=1e-2)
        np.testing.assert_allclose(out["v"], ce.tully_fisher_velocity(self.M), rtol=1e-2)

    def test_vacuum_sags_in_low_acceleration_regime(self):
        r = np.array([1e22])
        gN = ce.g_newton(self.M, r)
        self.assertLess(gN[0], ce.A0_MOND)
        Rv = ce.stiffness_profile(r, 3e19, gN)
        self.assertLess(Rv[0], 0.1 * ce.vacuum_stiffness_R0())

    def test_tully_fisher_eq7(self):
        v = ce.tully_fisher_velocity(self.M)
        self.assertAlmostEqual(v**4 / (G * self.M * ce.A0_MOND), 1.0, places=10)

    def test_velocity_is_flat_at_large_radius(self):
        r = np.array([1e21, 1e22, 1e23])
        v = ce.rotation_curve(r, self.M, rc=3e19)["v"]
        self.assertLess(np.ptp(v) / v.mean(), 0.01)


class TestHysteresis(unittest.TestCase):
    def test_relaxation_matches_ode(self):
        tau, Rs, R = 2.0, 1.0, 5.0
        dt, n = 1e-3, 3000
        R_euler = R
        for _ in range(n):
            R_euler += -(R_euler - Rs) / tau * dt
        R_exact = float(ce.relax_stiffness(R, Rs, dt * n, tau))
        self.assertAlmostEqual(R_euler, R_exact, places=2)
        self.assertAlmostEqual(R_exact, Rs + (R - Rs) * math.exp(-dt * n / tau), places=12)

    def test_instant_response_tracks_gas(self):
        r = ce.simulate_bullet_cluster(tau_v=1e-6)
        self.assertLess(abs(r.offset[-1]), 0.01)

    def test_eq8_alone_lags_the_gas(self):
        # résultat clé : sans transport, le minimum de R_v retarde sur le gaz
        self.assertLess(ce.simulate_bullet_cluster(tau_v=0.3).offset[-1], 0.0)

    def test_inertial_advection_extension_leads_the_gas(self):
        r = ce.simulate_bullet_cluster(tau_v=0.3, advect_with_inertia=True)
        self.assertGreater(r.offset[-1], 0.0)


if __name__ == "__main__":
    unittest.main()
