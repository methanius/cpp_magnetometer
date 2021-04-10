import numpy as np
import matplotlib.pyplot as plt
from math import ceil


def ehrenfest_chain(size, occupation):
    # Ehrenfest chain rates
    W_f = np.zeros(size)
    W_b = np.zeros(size)
    for i in range(size):
        W_f[i] = i / (size - 1) * (1 - occupation)
        W_b[i] = (1. - i / (size - 1)) * (1 - occupation)
    return W_f, W_b


T = 10000
dt = 1 / 200
t = np.arange(0, T, dt)
nT = len(t)
size = 16
delta_s = np.linspace(-2, 2, size)
delta_r = np.linspace(0, 0, size)
g = 2 * np.ones(size)
kappa = 0.1 * np.ones(size)
kappa_1 = kappa
beta = np.ones(size)
gamma_dec = np.linspace(0, 2, size)
gamma_phi = np.linspace(0, 2, size)
LOPhi = np.pi / 2
eta = 1
#seed = 902
occupation = 0.0001
emission_startpoint = int(ceil(size/2))

W_f, W_b = ehrenfest_chain(size, occupation)
J_f = np.diag(W_f[1:]**0.5, 1)
J_b = np.diag(W_b[:-1]**0.5, -1)
state = np.zeros((size, size))
state[size - 1 - emission_startpoint, size - 1 - emission_startpoint] = 1

# Run state evolution
state_list = np.zeros(nT)
state_track = emission_startpoint / (size - 1)
del_one = 1 / (size - 1)

for i in range(nT):
	permutation = np.random.rand()
	if permutation < 0.5:
		jump1 = np.random.rand()
		if jump1 < (J_f @ state @ J_f.T).trace() * dt:
			state = J_f @ state @ J_f.T
			state /= state.trace()
			state_track += del_one
		else:
			jump2 = np.random.rand()
			if jump2 < (J_b @ state @ J_b.T).trace() * dt:
				state = J_b @ state @ J_b.T
				state /= state.trace()
				state_track -= del_one
	else:
		jump1 = np.random.rand()
		if jump1 < (J_b @ state @ J_b.T).trace() * dt:
			state = J_b @ state @ J_b.T
			state /= state.trace()
			state_track -= del_one
		else:
			jump2 = np.random.rand()
			if jump2 < (J_f @ state @ J_f.T).trace() * dt:
				state = J_f @ state @ J_f.T
				state /= state.trace()
				state_track += del_one
	state_list[i] = state_track



'''
plt.figure()
plt.plot(t, state_list)
plt.figure()
plt.hist(state_list, bins=size)
plt.show()
'''