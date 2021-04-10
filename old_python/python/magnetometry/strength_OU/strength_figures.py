import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import os

points_per_unit_time = 0.1

image_type = ".pdf"

with open(sys.argv[1], 'rb') as f:
	[seed, T, dt, N, t, nT, base_rate, LOPhi, eta, delta_s, delta_r, g, kappa, kappa_1, beta, gamma_dec, gamma_phi, emission_startpoint,\
        W_up, W_down, dY, true_state, P, P_f, P_b, P_max, P_f_max] \
        = pickle.load(f)

image_folder = sys.argv[1][:-4] + '/'

points = int(1 / (points_per_unit_time * dt))

# Make folder for images
os.system('mkdir ' + image_folder)

#PQS probability density plot
plt.figure()
PQS_fig = plt.contourf(t[::points], delta_s, P[:, ::points], 100)
for c in PQS_fig.collections:
    c.set_edgecolor("face")
plt.plot(t[::points], P_max[::points], color='green', label='Maximum likelyhood guess')
plt.plot(t[::points], true_state[::points], color='red', label='"True" simulated B field')
plt.ylabel('$\Delta_s$')
plt.xlabel('Time [$T\gamma$]')
plt.title('PQS-treated probabilities')
plt.colorbar()
plt.legend()
plt.savefig(image_folder + 'PQS' + image_type, transparent=True, bbox_inches='tight')

#Forward-only probability density plot
plt.figure()
forward_fig = plt.contourf(t[::points], delta_s, P_f[:, ::points], 100)
for c in forward_fig.collections:
    c.set_edgecolor("face")
plt.plot(t[::points], P_f_max[::points], color='green', label='Maximum likelyhood guess')
plt.plot(t[::points], true_state[::points], color='red', label='"True" simulated B field')
plt.ylabel('$\Delta_s$')
plt.title('Forwarded propagated probabilities')
plt.xlabel('Time [$T\gamma$]')
plt.colorbar()
plt.legend()
plt.savefig(image_folder + 'forward' + image_type, transparent=True, bbox_inches='tight')



#B_true - PQS_max histogram data
PQS_dif_list = true_state - P_max
PQS_difs = np.zeros(2 * N - 1)
PQS_heights = np.zeros(2 * N - 1)
for i in range(2 * N - 1):
	if i < N - 1:
		PQS_difs[i] = delta_s[0] - delta_s[-1 -i]
	if i > N - 1:
		PQS_difs[i] = delta_s[-1] + delta_s[i - N + 1]
	PQS_heights[i] = (PQS_dif_list == PQS_difs[i]).sum()*dt

#B_true - PQS_forward historgram data
forward_dif_list = true_state - P_f_max
forward_difs = np.zeros(2 * N - 1)
forward_heights = np.zeros(2 * N - 1)
for i in range(2 * N - 1):
	if i < N - 1:
		forward_difs[i] = delta_s[0] - delta_s[-1 -i]
	if i > N - 1:
		forward_difs[i] = delta_s[-1] + delta_s[i - N + 1]
	forward_heights[i] = (forward_dif_list == forward_difs[i]).sum()*dt

#PQS histogram plot
plt.figure()
plt.bar(PQS_difs, PQS_heights, width=8/(2 * N - 1))
plt.xlabel('$B_{true} - B_{PQS ML}$')
plt.ylabel('Time at $B_{true} - B_{PQS ML}$')
plt.ylim((0, max(np.max(PQS_heights)*1.02 , np.max(forward_heights)*1.02)))
plt.title('$B_{true} - B_{PQS ML}$ occupation time')
plt.savefig(image_folder + 'PQS_histogram' + image_type, transparent=True, bbox_inches='tight')

#Forward histogram plot
plt.figure()
plt.bar(forward_difs, forward_heights, width=8/(2 * N - 1))
plt.xlabel('$B_{true} - B_{forward ML}$')
plt.ylabel('Time at $B_{true} - B_{forward ML}$')
plt.ylim((0, max(np.max(PQS_heights)*1.02 , np.max(forward_heights)*1.02)))
plt.title('$B_{true} - B_{forward ML}$ occupation time')
plt.savefig(image_folder + 'forward_histogram' + image_type, transparent=True, bbox_inches='tight')

plt.figure()
plt.plot(t, true_state - P_f_max, label='$B-B_{forward}$')
plt.plot(t, true_state - P_max, label='$B-B_{PQS}$')
plt.legend(loc="upper right")
plt.savefig(image_folder + 'B-Bi' + image_type)
plt.savefig(image_folder + 'forward_histogram' + image_type, transparent=True, bbox_inches='tight')



[seed, T, dt, N, t, nT, base_rate, LOPhi, eta, delta_s, delta_r, g, kappa, kappa_1, beta, gamma_dec, gamma_phi, emission_startpoint,\
        W_up, W_down, dY, true_state, P, P_f, P_b, P_max, P_f_max]


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


with open(image_folder + 'variances.txt', 'w') as f:
	f.write('var(B) = {}\n'.format(np.var(true_state)))
	f.write('var(B_f) = {}\n'.format(np.var(P_f_max)))
	f.write('var(B_PQS) = {}\n'.format(np.var(P_max)))

weighted_width_average_of_t = (P * delta_s[:, np.newaxis]).sum(axis=0)
weighted_width_variance_of_t = (P * (delta_s[:, np.newaxis] - weighted_width_average_of_t)**2).sum(axis=0)
weighted_t_averaged_width = (weighted_width_variance_of_t.sum()/len(weighted_width_variance_of_t))**0.5
errorbar = true_state - P_max
mean_errorbar = ((errorbar**2).sum()/len(errorbar))**0.5

weighted_width_average_of_t_forward = (P_f * delta_s[:, np.newaxis]).var(axis=0)
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
