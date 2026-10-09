import numpy as np
from .constants import E_L

def emergent_time_from_energy(E_values: np.ndarray,
                              E_L_value: float = E_L) -> np.ndarray:
    dE = np.diff(E_values)
    T = [0.0]
    total = 0.0
    for dE_i in dE:
        total += dE_i / E_L_value
        T.append(total)
    return np.array(T)

