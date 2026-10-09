import numpy as np
from .constants import E_L

def dt_from_entropy_change(dS: float, T: float, E_L_value: float = E_L) -> float:
    dE = T * dS
    return dE / E_L_value

def integrate_thermo_time(S_values: np.ndarray,
                          T_values: np.ndarray,
                          E_L_value: float = E_L) -> np.ndarray:
    total_time = 0.0
    times = [0.0]
    for i in range(len(S_values) - 1):
        dS = S_values[i+1] - S_values[i]
        T_mid = 0.5 * (T_values[i] + T_values[i+1])
        total_time += dt_from_entropy_change(dS, T_mid, E_L_value)
        times.append(total_time)
    return np.array(times)
