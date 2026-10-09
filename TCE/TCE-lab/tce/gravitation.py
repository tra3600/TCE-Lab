import numpy as np
from .constants import c, G

class EnergeticSystem:
    def __init__(self, m, E_internal=0.0, E_kinetic=0.0, E_potential=0.0):
        self.m = m
        self.E_internal = E_internal
        self.E_kinetic = E_kinetic
        self.E_potential = E_potential

    def total_energy(self):
        return self.m * c**2 + self.E_internal + self.E_kinetic + self.E_potential

    def available_for_time(self):
        return self.total_energy() - self.E_kinetic - self.E_potential

def gravitational_potential_energy(M, m, r):
    return -G * M * m / r

def gravitational_time_dilation_TCE(system, M, r):
    E_p = gravitational_potential_energy(M, system.m, r)
    system_in_well = EnergeticSystem(
        m=system.m,
        E_internal=system.E_internal,
        E_kinetic=system.E_kinetic,
        E_potential=E_p
    )
    E_time_far = system.available_for_time()
    E_time_near = system_in_well.available_for_time()
    return E_time_near / E_time_far

