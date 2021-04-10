import numpy as np

def ehrenfest_chain(size, base_rate):
    # Ehrenfest chain rates
    W_up = np.zeros(size)
    W_down = np.zeros(size)
    for i in range(0, size):
        W_up[i] = i / (size - 1) * base_rate
        W_down[i] = (1. - i / (size - 1)) * base_rate
    return W_up[1:], W_down[:-1]
