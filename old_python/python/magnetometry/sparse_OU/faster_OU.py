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


if __name__ == '__main__':    

    show_plots = 0

    #Parameters

    seed = np.random.randint(2**32)
    np.random.seed(seed)

    T = int(sys.argv[1])
    dt = 1 / 200
    N = int(sys.argv[2])
    t = np.arange(0, T, dt)
    nT = len(t)
    base_rate = float(sys.argv[3])
    LOPhi = np.pi / 2
    eta = 1 
    delta_s = np.linspace(-2, 2, N)
    delta_r = np.linspace(0, 0, N)
    g = np.linspace(2, 2, N)
    kappa = np.linspace(10, 10, N)
    kappa_1 = kappa
    beta = np.linspace(3, 3, N)
    gamma_dec = np.linspace(0.1, 0.1, N)
    gamma_phi = np.linspace(0.1, 0.1, N)
    emission_startpoint = int(N / 2)
    auto_T = float(sys.argv[4])

    W_up, W_down = ehrenfest_chain.ehrenfest_chain(N, base_rate)
    c_1, c_2, c_3, c, H, Jup, Jdown = sparse_OU_library.make_operators(N, g, kappa, delta_r, delta_s, kappa_1, beta, gamma_dec, gamma_phi, LOPhi, W_up, W_down)
    dY, true_state = sparse_OU_library.homodyne_emission(N, H, Jup, Jdown, c, c_1, c_2, c_3, dt, nT, eta, emission_startpoint, delta_s)
"""
    P, P_f, P_b, P_max, P_f_max = sparse_OU_library.homodyne_PQS(nT, N, H, c, Jdown, Jup, c_1, c_2, c_3, dt, dY, eta, delta_s)
    t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation = sparse_OU_library.autocorrelation(auto_T, dt, P_max, true_state, t, nT)
    averaged_width, mean_errorbar = sparse_OU_library.B_PQS_var_vs_width(true_state, P_max, P)

    with open('N_{0}_T_{1}_r_{2}_'.format(N, T, base_rate) + datetime.datetime.now().strftime("%H_%M_%S_%d_%m_%y") + ".pkl", 'wb') as f:
        pickle.dump([seed, T, dt, N, t, nT, base_rate, LOPhi, eta, delta_s, delta_r, g, kappa, kappa_1, beta, gamma_dec, gamma_phi, emission_startpoint,\
        auto_T, W_up, W_down, dY, true_state, P, P_f, P_b, P_max, P_f_max, t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation], f)#, protocol=pickle.HIGHEST_PROTOCOL)


    if show_plots == True:
        plt.figure()
        plt.subplot(4, 1, 1)
        plt.plot(t, true_state)
        plt.plot(t, P_max)
        labels = []
        ticklocation = np.linspace(-2, 2, N)
        for i in range(N):
            labels.append("P{0}".format(i))
        plt.yticks(ticklocation, labels)
        plt.subplot(4, 1, 2)
        for i in range(N):
            plt.plot(t, P[i, :], label='P%s' % (i))
        plt.legend()
        plt.subplot(4, 1, 3)
        for i in range(N):
            plt.plot(t, P_f[i, :], label='P_f%s' % (i))
        plt.legend()
        plt.subplot(4, 1, 4)
        for i in range(N):
            plt.plot(t, P_b[i, :], label='P_b%s' % (i))
        plt.legend()
        plt.show()
"""
