import numpy as np
from dataclasses import dataclass

from .constants import c

@dataclass
class SpinFieldEnergyState:
    E: float
    psi: np.ndarray
    psi_prime: np.ndarray

class Dirac1DOperator:
    def __init__(self, dx: float, m: float = 1.0):
        self.dx = dx
        self.m = m
        self.alpha = np.array([[0.0, 1.0],
                               [1.0, 0.0]])
        self.beta = np.array([[1.0, 0.0],
                              [0.0, -1.0]])

    def spatial_derivative(self, psi: np.ndarray) -> np.ndarray:
        grad = np.zeros_like(psi)
        grad[1:-1] = (psi[2:] - psi[:-2]) / (2 * self.dx)
        grad[0] = (psi[1] - psi[0]) / self.dx
        grad[-1] = (psi[-1] - psi[-2]) / self.dx
        return grad

    def __call__(self, psi: np.ndarray) -> np.ndarray:
        dpsi_dx = self.spatial_derivative(psi)
        term_space = np.einsum("ij,xj->xi", self.alpha, dpsi_dx) * c
        term_mass = np.einsum("ij,xj->xi", self.beta, psi) * (self.m * c**2)
        return term_space + term_mass

class SpinFieldEnergyLagrangian:
    def __init__(self, H_spatial_operator: Dirac1DOperator, dx: float):
        self.H_spatial = H_spatial_operator
        self.dx = dx

    def lagrangian_density(self, state: SpinFieldEnergyState) -> np.ndarray:
        psi = state.psi
        psi_prime = state.psi_prime
        Hpsi = self.H_spatial(psi)
        hermitian_time = np.sum(np.conjugate(psi) * psi_prime, axis=-1).real
        hermitian_space = np.sum(np.conjugate(psi) * Hpsi, axis=-1).real
        return hermitian_time - hermitian_space

    def volume_element(self) -> float:
        return self.dx

    def L_energy(self, state: SpinFieldEnergyState) -> float:
        density = self.lagrangian_density(state)
        return np.sum(density) * self.volume_element()
