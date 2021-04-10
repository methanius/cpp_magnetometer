import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import os

points_per_unit_time = 0.1

image_folder = ''


with open(sys.argv[1], 'rb') as f:
	[seed, T, dt, N, t, nT, base_rate, LOPhi, eta, delta_s, delta_r, g, kappa, kappa_1, beta, gamma_dec, gamma_phi, emission_startpoint,\
        auto_T, W_up, W_down, dY, true_state, P, P_f, P_b, P_max, P_f_max, t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation] \
        = pickle.load(f)

points = int(1 / (points_per_unit_time * dt))


weighted_width_average_of_t = (P * delta_s[:, np.newaxis]).sum(axis=0)
weighted_width_variance_of_t = (P * (delta_s[:, np.newaxis] - weighted_width_average_of_t)**2).sum(axis=0)
weighted_t_averaged_width = (weighted_width_variance_of_t.sum()/len(weighted_width_variance_of_t))**0.5
errorbar = true_state - P_max
mean_errorbar = ((errorbar**2).sum()/len(errorbar))**0.5

weighted_width_average_of_t_forward = (P_f * delta_s[:, np.newaxis]).sum(axis=0)
weighted_width_variance_of_t_forward = (P_f * (delta_s[:, np.newaxis] - weighted_width_average_of_t_forward)**2).sum(axis=0)
weighted_t_averaged_width_forward = (weighted_width_variance_of_t_forward.sum()/len(weighted_width_variance_of_t_forward))**0.5
errorbar_forward = true_state - P_f_max
mean_errorbar_forward = ((errorbar_forward**2).sum()/len(errorbar_forward))**0.5


with open('widths_' + sys.argv[1][:sys.argv[1].find('/')] + '.pkl', 'wb') as f:
    pickle.dump([weighted_t_averaged_width, mean_errorbar, weighted_t_averaged_width_forward, mean_errorbar_forward],f)

