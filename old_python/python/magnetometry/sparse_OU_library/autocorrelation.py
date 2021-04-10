import numpy as np


def autocorrelation(auto_T, dt, P_max, true_state, t, nT):
    t_autocorrelation = np.zeros(int(auto_T / dt))
    P_max_autocorrelation = np.zeros(int(auto_T / dt))
    true_state_autocorrelation = np.zeros(int(auto_T / dt))
    padding = 200
    auto_start = int(padding / dt)
    auto_stop = int(nT - padding / dt)
    if nT < (auto_T + 2 * padding) / dt:
        print("Not enough time given for autocorrelation")
        return t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation
    for i in range(int(auto_T / dt)):
        P_max_autocorrelation[i] = P_max[auto_start + i:auto_stop] @ P_max[auto_start:auto_stop - i] / len(P_max[auto_start + i:auto_stop])
        true_state_autocorrelation[i] = true_state[auto_start + i:auto_stop] @ true_state[auto_start:auto_stop - i] / len(true_state[auto_start + i:auto_stop])
        t_autocorrelation[i] = i * dt
    return t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation
