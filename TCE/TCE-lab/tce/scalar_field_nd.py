import numpy as np
from dataclasses import dataclass

from .constants import E_L

@dataclass
class FieldEnergyStateND:
    E: float
    phi: np.ndarray
    phi_prime: np.ndarray

class ScalarFieldEnergyLagrangianND:
    def __init__(self, shape, spacing, m: float = 1.0):
        self.shape = shape
        self.spacing = spacing
        self.m = m

    def gradient(self, phi: np.ndarray) -> np.ndarray:
        D = len(self.shape)
        grads = []
        for axis in range(D):
            dx = self.spacing[axis]
            grad_axis = np.zeros_like(phi)
            slc_center = [slice(None)] * D
            slc_plus = [slice(None)] * D
            slc_minus = [slice(None)] * D
            slc_center[axis] = slice(1, -1)
            slc_plus[axis] = slice(2, None)
            slc_minus[axis] = slice(None, -2)
            grad_axis[tuple(slc_center)] = (
                phi[tuple(slc_plus)] - phi[tuple(slc_minus)]
            ) / (2 * dx)
            slc0 = [slice(None)] * D
            slc1 = [slice(None)] * D
            slc0[axis] = 0
            slc1[axis] = 1
            grad_axis[tuple(slc0)] = (phi[tuple(slc1)] - phi[tuple(slc0)]) / dx
            slcN1 = [slice(None)] * D
            slcN2 = [slice(None)] * D
            slcN1[axis] = -1
            slcN2[axis] = -2
            grad_axis[tuple(slcN1)] = (phi[tuple(slcN1)] - phi[tuple(slcN2)]) / dx
            grads.append(grad_axis)
        return np.stack(grads, axis=0)

    def lagrangian_density(self, state: FieldEnergyStateND) -> np.ndarray:
        phi = state.phi
        phi_prime = state.phi_prime
        grad_phi = self.gradient(phi)
        grad_sq = np.sum(grad_phi**2, axis=0)
        return 0.5 * phi_prime**2 - 0.5 * grad_sq - 0.5 * self.m**2 * phi**2

    def volume_element(self) -> float:
        return np.prod(self.spacing)

    def L_energy(self, state: FieldEnergyStateND) -> float:
        density = self.lagrangian_density(state)
        return np.sum(density) * self.volume_element()
