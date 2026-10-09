import numpy as np
from typing import Callable, Tuple, List

def rk4_step(f: Callable[[float, np.ndarray], np.ndarray],
             E: float,
             y: np.ndarray,
             dE: float) -> np.ndarray:
    k1 = f(E, y)
    k2 = f(E + 0.5 * dE, y + 0.5 * dE * k1)
    k3 = f(E + 0.5 * dE, y + 0.5 * dE * k2)
    k4 = f(E + dE, y + dE * k3)
    return y + (dE / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

def integrate_field_energy(f: Callable[[float, np.ndarray], np.ndarray],
                           E0: float,
                           y0: np.ndarray,
                           E1: float,
                           n_steps: int) -> Tuple[np.ndarray, List[np.ndarray]]:
    Es = np.linspace(E0, E1, n_steps + 1)
    ys = [y0.copy()]
    y = y0.copy()
    dE = (E1 - E0) / n_steps
    E = E0
    for _ in range(n_steps):
        y = rk4_step(f, E, y, dE)
        E += dE
        ys.append(y.copy())
    return Es, ys
