import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import os

points_per_unit_time = 0.05

image_type = ".pdf"
image_folder = sys.argv[1][:-4] + '/'
os.system('mkdir ' + image_folder)


with open(sys.argv[1], 'rb') as f:
        [seed, T, dt, N, t, nT, base_rate, LOPhi, eta, delta_s, delta_r, g, kappa, kappa_1, beta, gamma_dec, gamma_phi, emission_startpoint,\
        auto_T, W_up, W_down, dY, true_state, P, P_f, P_b, P_max, P_f_max, t_autocorrelation, true_state_autocorrelation, P_max_autocorrelation] \
        = pickle.load(f)

points = int(1 / (points_per_unit_time * dt))

clevels = 180
#PQS probability density plot
plt.figure()
PQS_fig = plt.contourf(t[::points], delta_s, P[:, ::points], clevels)
for c in PQS_fig.collections:
    c.set_edgecolor("face")
plt.plot(t[::points], P_max[::points], color='green', label='MLE')
plt.plot(t[::points], true_state[::points], color='red', label='True field')
plt.ylabel('$\Delta_n / \gamma_{dec}$', fontsize=18)
plt.xlabel('$\gamma_{dec} t$', fontsize=18)
plt.xticks([3500, 7000, 10500, 14000], fontsize=18)
plt.yticks(fontsize=18)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=18)
plt.legend(loc="upper right", fontsize=18)
plt.savefig(image_folder + 'PQS' + image_type, transparent=True, bbox_inches='tight')

#Forward-only probability density plot
plt.figure()
forward_fig = plt.contourf(t[::points], delta_s, P_f[:, ::points], clevels)
for c in forward_fig.collections:
    c.set_edgecolor("face")
plt.plot(t[::points], P_f_max[::points], color='green', label='MLE')
plt.plot(t[::points], true_state[::points], color='red', label='True field')
plt.ylabel('$\Delta_n / \gamma_{dec}$', fontsize=18)
plt.xlabel('$\gamma_{dec} t$', fontsize=18)
plt.xticks([3500, 7000, 10500, 14000], fontsize=18)
plt.yticks(fontsize=18)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=18)
plt.legend(loc="upper right", fontsize=18)
plt.savefig(image_folder + 'forward' + image_type, transparent=True, bbox_inches='tight')

#Forward-only probability density plot
plt.figure()
forward_fig = plt.contourf(t[::points], delta_s, P_b[:, ::points], clevels)
for c in forward_fig.collections:
    c.set_edgecolor("face")
plt.plot(t[::points], P_f_max[::points], color='green', label='ML estimator')
plt.plot(t[::points], true_state[::points], color='red', label='True field')
plt.ylabel('$\Delta_n / \gamma_{dec}$', fontsize=18)
plt.xlabel('$\gamma_{dec} t$', fontsize=18)
plt.xticks([3500, 7000, 10500, 14000], fontsize=18)
plt.yticks(fontsize=18)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=18)
plt.legend(loc="upper right", fontsize=18)
plt.savefig(image_folder + 'backward' + image_type, transparent=True, bbox_inches='tight')

#B_true - PQS_max histogram data
PQS_dif_list = true_state - P_max
PQS_heights = np.zeros_like(delta_s)
for i in range(len(delta_s)):
        PQS_heights[i] = (np.around(PQS_dif_list, 2) == np.around(delta_s[i], 2)).sum()*dt


#B_true - PQS_forward historgram data
forward_dif_list = true_state - P_f_max
forward_heights = np.zeros_like(delta_s)
for i in range(len(delta_s)):
        forward_heights[i] = (np.around(forward_dif_list, 2) == np.around(delta_s[i], 2)).sum()*dt

print(PQS_heights.sum())
ss_heights = np.zeros_like(delta_s)
for i in range(len(delta_s)):
        ss_heights[i] = (true_state == delta_s[i]).sum()*dt

#PQS histogram plot
plt.figure()
plt.grid()
plt.bar(delta_s, ss_heights, width=8/(2 * N - 1), color = 'xkcd:orange', alpha=0.4, label='Steady state MLE' )
plt.bar(delta_s, PQS_heights, width=8/(2 * N - 1), alpha=0.6, label='PQS MLE')
plt.xlabel('$\Delta_{n,True} - \Delta_{n,PQS MLE}$', fontsize=18)
plt.ylabel('Time at $\Delta_{n,True} - \Delta_{n,MLE}$ in units of $\gamma t$', fontsize=18)
plt.ylim((0, max(np.max(PQS_heights)*1.02 , np.max(forward_heights)*1.02)))
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.legend(fontsize=12)
plt.savefig(image_folder + 'PQS_histogram' + image_type, transparent=True, bbox_inches='tight')

#Forward histogram plot
plt.figure()
plt.grid()
plt.bar(delta_s, ss_heights, width=8/(2 * N - 1), color = 'xkcd:orange', alpha=0.4, label='Steady state MLE' )
plt.bar(delta_s, forward_heights, width=8/(2 * N - 1), alpha=0.6, label='Filtered state MLE')
plt.xlabel('$\Delta_{n,True} - \Delta_{n,filtered MLE}$', fontsize=18)
plt.ylabel('Time at $\Delta_{n,True} - \Delta_{n,MLE}$ in units of $\gamma t$', fontsize=18)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.ylim((0, max(np.max(PQS_heights)*1.02 , np.max(forward_heights)*1.02)))
plt.legend(fontsize=12)
plt.savefig(image_folder + 'forward_histogram' + image_type, transparent=True, bbox_inches='tight')


with open(image_folder + 'variables.txt','w') as f:
        f.write('seed = {}\n'.format(seed))
        f.write('T = {}\n'.format(T))
        f.write('dt = {}\n'.format(dt))
        f.write('N = {}\n'.format(N))
        f.write('base_rate = {}\n'.format(base_rate))
        f.write('LOPhi = {}\n'.format(LOPhi))
        f.write('eta = {}\n'.format(eta))
        f.write('delta_s = linspace({}, {}, N)\n'.format(np.min(delta_s), np.max(delta_s)))
        f.write('delta_r= linspace({}, {}, N)\n'.format(np.min(delta_r), np.max(delta_r)))
        f.write('g = linspace({}, {}, N)\n'.format(np.min(g), np.max(g)))
        f.write('kappa = linspace({}, {}, N)\n'.format(np.min(kappa), np.max(kappa)))
        f.write('kappa_1 = linspace({}, {}, N)\n'.format(np.min(kappa_1), np.max(kappa_1)))
        f.write('beta = linspace({}, {}, N)\n'.format(np.min(kappa), np.max(kappa)))
        f.write('gamma_dec = linspace({}, {}, N)\n'.format(np.min(gamma_dec), np.max(gamma_dec)))
        f.write('gamma_phi = linspace({}, {}, N)\n'.format(np.min(gamma_phi), np.max(gamma_phi)))
        f.write('emission_startpoint = {}\n'.format(emission_startpoint))
        f.write('auto_T = {}\n'.format(auto_T))

with open(image_folder + 'variances.txt', 'w') as f:
        f.write('var(B) = {}\n'.format(np.var(true_state)))
        f.write('var(B_f) = {}\n'.format(np.var(P_f_max)))
        f.write('var(B_PQS) = {}\n'.format(np.var(P_max)))

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


with open(image_folder + 'width_vs_errorbar.txt', 'w') as f:
        f.write('PQS Averaged width = {}\n'.format(weighted_t_averaged_width))
        f.write('PQS Mean errorbar = {}\n'.format(mean_errorbar))
        f.write('\n')
        f.write('Forward Averaged width = {}\n'.format(weighted_t_averaged_width_forward))
        f.write('Forward Mean errorbar = {}\n'.format(mean_errorbar_forward))


plt.figure()
plt.plot(t, weighted_width_variance_of_t, label='PQS variance')
plt.plot(t, weighted_width_variance_of_t_forward, label='forward variance')
plt.legend(loc="upper right")
plt.title('P(t) variances as functions of time')
plt.savefig(image_folder + 'variance_through_time.pdf')

os.system('mv ' + sys.argv[1] + ' ' + image_folder)
