import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom


def q0(p, t):
	return 0.5 * (1 - (1 - 2 * p)**t)


def q1(p, t):
	return 0.5 * (1 + (1 - 2 * p)**t)


N = 10
show_steady_plot = 0
show_autocorrelation_plot = 1
p = 0.03
T = 50 # Integer
time = np.linspace(1, T, T)


# Steady state
x = np.linspace(0, N, N + 1)
P_steady_state = 2**(-N) * binom(N, x)
x = (x / N - 0.5) * 4 

# Autocorrelation
steady_state_weights = x * P_steady_state

PT = np.zeros((T + 1, N + 1))
for i in range(N + 1):
	PT[0] = (P_steady_state * x) 


for t in time:
	for NA in range(N + 1):
		N0A_sum = 0
		for N0A in range(N + 1):
			i_sum = 0
			for i in range(NA + 1):
				i_sum += binom(N0A, i) * binom(N - N0A, NA - i) * (q1(p, t) / q0(p, t))**(2 * i)
			i_sum *= q1(p, t)**N * (q0(p, t) / q1(p, t))**(NA + N0A)
			N0A_sum += i_sum * P_steady_state[N0A] * x[N0A]
		PT[int(t), NA] = N0A_sum * (NA / N - 0.5) * 4

autocorrelation = np.zeros(T + 1)
t_autocorrelation = np.arange(T + 1) / (N * p)
for i in range(T + 1):
	autocorrelation[i] = PT[i].sum()

if show_steady_plot == True:
	plt.figure()
	plt.plot(x, P_steady_state)
	
if show_autocorrelation_plot == True:
	plt.figure()
	plt.plot(t_autocorrelation, autocorrelation)

if show_steady_plot or show_autocorrelation_plot == 1:
	plt.show()