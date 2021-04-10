import numpy as np
from scipy import sparse
import sys
sys.path.append('..')
import sparse_OU_library
import matplotlib.pyplot as plt

'''
def autocorrelation(auto_T, dt, P_max, true_state, t, mod_data):
    auto_correlation_frames = int(auto_T /dt / mod_data)
    t_autocorrelation = np.zeros(auto_correlation_frames)
    P_max_autocorrelation = np.zeros(auto_correlation_frames)
    true_state_autocorrelation = np.zeros(auto_correlation_frames)
    padding = int(200 / dt / mod_data)
    auto_start = int(padding)
    auto_stop = int(- 1 - padding)
    if auto_correlation_frames > (t[-1] + 1) / dt / mod_data - 2 * padding:
        print("Not enough time given for autocorrelation")
        return t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation
    for i in range(auto_correlation_frames):
        P_max_autocorrelation[i] = P_max[auto_start + i:auto_stop] @ P_max[auto_start:auto_stop - i] / t[len(P_max[auto_start + i:auto_stop])]
        true_state_autocorrelation[i] = true_state[auto_start + i:auto_stop] @ true_state[auto_start:auto_stop - i] / t[len(true_state[auto_start + i:auto_stop])]
        t_autocorrelation[i] = i * dt * mod_data
    return t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation
'''


if __name__ == '__main__':    
    #Parameters
    seed = np.random.randint(2**32)
    np.random.seed(seed)

    show_plot = 1

    T = 100 # int(sys.argv[1])
    dt = 1 / 200
    size = 10 # int(sys.argv[2])
    t = np.arange(0, T, dt)
    nT = len(t)

    mod_data = 1 #float(sys.argv[4]) #Saved data slices per unit rescaled time

    frame_steps = np.linspace(0, nT - 1, T * mod_data, dtype=int)
    t_frames = t[frame_steps]

    base_rate = 0.03 # float(sys.argv[3])
    LOPhi = np.pi / 2
    eta = 1
    delta_s = np.linspace(-2, 2, size)
    delta_r = np.linspace(0, 0, size)
    g = np.linspace(2, 2, size)
    kappa = np.linspace(10, 10, size)
    kappa_1 = kappa
    beta = np.linspace(3, 3, size)
    gamma_dec = np.linspace(0, 2, size)
    gamma_phi = np.linspace(0, 2, size)
    emission_startpoint = int(np.random.rand()*size)
    #auto_T = float(sys.argv[4])

    W_up = np.zeros(size)
    W_down = np.zeros(size)
    c_1, c_2, c_3, c, H, Jup, Jdown = sparse_OU_library.make_operators(size, g, kappa, delta_r, delta_s, kappa_1, beta, gamma_dec, gamma_phi, LOPhi, W_up, W_down)
    dY, true_state = sparse_OU_library.homodyne_emission(size, H, Jup, Jdown, c, c_1, c_2, c_3, dt, nT, eta, emission_startpoint, delta_s)
    P, P_f, P_b, P_max, P_f_max = sparse_OU_library.framed_homodyne_PQS(nT, size, H, c, Jdown, Jup, c_1, c_2, c_3, dt, dY, eta, frame_steps, delta_s)

    if show_plot ==1:
        plt.figure()
        plt.subplot(4, 1, 1)
        plt.plot(t, true_state)
        plt.plot(t_frames, P_max)
        labels = []
        ticklocation = np.linspace(-2, 2, size)
        #for i in range(size):
        #    labels.append("P{0}".format(i))
        #plt.yticks(ticklocation, labels)
        plt.subplot(4, 1, 2)
        for i in range(size):
            plt.plot(t_frames, P[i, :], label='P%s' % (i))
        #plt.legend()
        plt.subplot(4, 1, 3)
        for i in range(size):
            plt.plot(t_frames, P_f[i, :], label='P_f%s' % (i))
        #plt.legend()
        plt.subplot(4, 1, 4)
        for i in range(size):
            plt.plot(t_frames, P_b[i, :], label='P_b%s' % (i))
        #plt.legend()
        plt.show()
