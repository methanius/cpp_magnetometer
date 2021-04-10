import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl


N = 100
p = 0.3
T = 5
dt = 1/50
show_plot = 1


U = np.diag(np.ones(N + 1) - (p * dt) * N) + np.diag(np.arange(1, N + 1, 1) * (p * dt), 1) + np.diag(np.arange(1, N + 1, 1)[::-1] * (p * dt), -1)



first_projection = np.zeros(N + 1, dtype=float)
second_projection = np.zeros_like(first_projection)
third_projection = np.zeros_like(first_projection)

first_projection[0] = 1
second_projection[int((N + 1) / 2)] = 1
third_projection[-1] = 1

nT = int((T + 1)/dt)
# Propagation
forward_state = np.ones((nT, N + 1))
backward_state = np.ones_like(forward_state)

for t in range(0, nT):
        if t == 0:
            forward_state[t] = first_projection
        elif t == np.ceil(nT / 2):
            forward_state[t] = second_projection
        elif t == nT - 1:
            forward_state[t] = third_projection
        else:
            forward_state[t] = U @ forward_state[t - 1]
            forward_state[t] /= forward_state[t].sum()

for t in range(0, nT):
        if t == 0:
            backward_state[t] = third_projection
        elif t == np.floor(nT / 2):
            backward_state[t] = second_projection
        elif t == nT - 1:
            backward_state[t] = first_projection
        else:
            backward_state[t] = U.T @ backward_state[t - 1]
            backward_state[t] /= backward_state[t].sum()

backward_state = backward_state[::-1]
forward_backward_state = np.zeros_like(forward_state)
for t in range(nT):
        forward_backward_state[t] = forward_state[t] * backward_state[t] / (forward_state[t] * backward_state[t]).sum()

colourmap = 'coolwarm'
levels = 400

if show_plot == 1:
        x = np.linspace(0, N, N + 1)
        t = np.linspace(0, T, nT, dtype=float)
        plt.figure()
        forward = plt.contourf(t, x, forward_state.T, levels, cmap=colourmap)
        plt.title('Forward estimate of the example record')
        plt.xlabel('Time in arbitrary units')
        plt.ylabel('Hidden states')
        plt.colorbar()
        plt.figure()
        backward = plt.contourf(t, x, backward_state.T, levels, cmap=colourmap)
        plt.title('Backward estimate of the example record')
        plt.xlabel('Time in arbitrary units')
        plt.ylabel('Hidden states')
        plt.colorbar()
        plt.figure()
        forward_backward = plt.contourf(t, x, forward_backward_state.T, levels, cmap=colourmap)
        plt.title('Forward-backward estimate of the example record')
        plt.xlabel('Time in arbitrary units')
        plt.ylabel('Hidden states')
        plt.colorbar()
        plt.show()
