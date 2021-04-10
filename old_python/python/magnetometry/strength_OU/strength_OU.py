import numpy as np
import scipy.linalg as spl
from scipy import sparse
import matplotlib.pyplot as plt
import time
import sys
sys.path.append('..')
import ehrenfest_chain
import sparse_OU_library
import pickle
import datetime



#Parameters
seed = np.random.randint(2**32)
np.random.seed(seed)
T = float(sys.argv[2])
dt = 1 / 200
N = int(sys.argv[3])
t = np.arange(0, T, dt)
nT = len(t)
base_rate = float(sys.argv[4])
LOPhi = np.pi / 2
eta = 1 
delta_s = np.linspace(-2, 2, N)
delta_r = np.linspace(0, 0, N)
g = np.linspace(2, 2, N)
kappa = float(sys.argv[1]) * np.ones(N)
kappa_1 = kappa
beta = np.linspace(3, 3, N)
gamma_dec = np.linspace(0, 2, N)
gamma_phi = np.linspace(0, 2, N)
emission_startpoint = int(N / 2)

W_up, W_down = ehrenfest_chain.ehrenfest_chain(N, base_rate)
c_1, c_2, c_3, c, H, Jup, Jdown = sparse_OU_library.make_operators(N, g, kappa, delta_r, delta_s, kappa_1, beta, gamma_dec, gamma_phi, LOPhi, W_up, W_down) 
dY, true_state = sparse_OU_library.homodyne_emission(N, H, Jup, Jdown, c, c_1, c_2, c_3, dt, nT, eta, emission_startpoint, delta_s)
P, P_f, P_b, P_max, P_f_max = sparse_OU_library.homodyne_PQS(nT, N, H, c, Jdown, Jup, c_1, c_2, c_3, dt, dY, eta, delta_s)


with open('K_{0}_N_{1}_T_{2}_r_{3}_'.format(sys.argv[1], N, T, base_rate) + datetime.datetime.now().strftime("%H_%M_%S_%d_%m_%y") + ".pkl", 'wb') as f:
    pickle.dump([seed, T, dt, N, t, nT, base_rate, LOPhi, eta, delta_s, delta_r, g, kappa, kappa_1, beta, gamma_dec, gamma_phi, emission_startpoint,\
    W_up, W_down, dY, true_state, P, P_f, P_b, P_max, P_f_max], f)#, protocol=pickle.HIGHEST_PROTOCOL)
