import numpy as np
import matplotlib.pyplot as plt
import sys
import os
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})
plt.rc("axes", labelsize=18, titlesize=22)   # skriftstørrelse af `xlabel`, `ylabel` og `title`
plt.rc("xtick", labelsize=16, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=16, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=16) # skriftstørrelse af figurers legends

dpi = 500
colormap = 'Blues'
MLE_line_color = 'darkblue'
True_field_line_color = 'crimson' #'firebrick'
contourmax = 0.475
clevels = 25 #100


image_type = '.pdf'
target_directory = os.getcwd() + '/' + sys.argv[1]
if(target_directory[-1] != '/'):
    target_directory += '/'
in_params = target_directory + 'input_parameters/'
dt = np.loadtxt(in_params + 'dt.txt')
delta_q = np.loadtxt(in_params + 'delta_q.txt')
N = np.loadtxt(in_params + 'N.txt')
T = np.loadtxt(in_params + 'T.txt')
frames = np.loadtxt(in_params + 'frames.txt')


in_data = target_directory + 'data/'
t_emission = np.loadtxt(in_data + 't_emission.txt')
t_PQS = np.loadtxt(in_data + 't_PQS.txt')
true_state = np.loadtxt(in_data + 'true_state.txt')
index_frames_in_emission = np.isin(t_emission, t_PQS)

P = np.loadtxt(in_data + 'P.txt')
P_f = np.loadtxt(in_data + 'P_f.txt')
P_b = np.loadtxt(in_data + 'P_b.txt')
P_max = np.loadtxt(in_data + 'P_max.txt').reshape(-1)
P_f_max = np.loadtxt(in_data + 'P_f_max.txt').reshape(-1)

results_folder = target_directory + 'analysis/'
os.system('mkdir ' + results_folder)

for i, state in enumerate(true_state):
    true_state[i] = delta_q[int(state)]

for i in range(len(P_max)):
    P_max[i] = delta_q[int(P_max[i])]
    P_f_max[i] = delta_q[int(P_f_max[i])]


PQS_weighted_averages = [np.average(delta_q, weights = weights) for weights in P.T]
PQS_weighted_variances = np.array([(weights * (delta_q - pdf_mean)**2).sum() for weights, pdf_mean in zip(P.T, PQS_weighted_averages)])
averaged_widths = np.mean(PQS_weighted_variances)**0.5
errorbar = true_state[index_frames_in_emission] - P_max
mean_errorbar = ((errorbar**2).mean())**0.5

forward_weighted_averages = [np.average(delta_q, weights = weights) for weights in P_f.T]
forward_weighted_variances = np.array([(weights * (delta_q - pdf_mean)**2).sum() for weights, pdf_mean in zip(P_f.T, forward_weighted_averages)])
averaged_widths_forward = (forward_weighted_variances.mean())**0.5
errorbar_forward = true_state[index_frames_in_emission] - P_f_max
mean_errorbar_forward = ((errorbar_forward**2).mean())**0.5


with open(results_folder + 'width_vs_errorbar.txt', 'w') as f:
        f.write('PQS Averaged width = {}\n'.format(averaged_widths))
        f.write('PQS Mean errorbar = {}\n'.format(mean_errorbar))
        f.write('\n')
        f.write('Forward Averaged width = {}\n'.format(averaged_widths_forward))
        f.write('Forward Mean errorbar = {}\n'.format(mean_errorbar_forward))


#Forward-only probability density plot
plt.figure()
forward_fig = plt.contourf(t_PQS, delta_q, P_f, clevels, cmap = colormap, vmax = contourmax)
for c in forward_fig.collections:
    c.set_edgecolor("face")
plt.plot(t_PQS, P_f_max, color=MLE_line_color, ls = '-', label='MLE')
plt.plot(t_emission, true_state, color=True_field_line_color, label='True field')
plt.ylabel('$\Delta_n / \gamma$')
plt.xlabel('$\gamma t$')
plt.yticks()
cbar = plt.colorbar()
cbar.ax.tick_params()
plt.legend(loc="upper right")
plt.savefig(results_folder + 'filtering' + image_type, transparent=True, bbox_inches='tight', dpi = dpi)


#PQS probability density plot
plt.figure()
PQS_fig = plt.contourf(t_PQS, delta_q, P,clevels, cmap = colormap, vmax = contourmax)
for c in PQS_fig.collections:
    c.set_edgecolor("face")
plt.plot(t_PQS, P_max, color=MLE_line_color, label='MLE')
plt.plot(t_emission, true_state, color=True_field_line_color, label='True field')
plt.ylabel('$\Delta_n / \gamma$')
plt.xlabel('$\gamma t$')
plt.yticks()
cbar = plt.colorbar()
cbar.ax.tick_params()
plt.legend(loc="upper right")
plt.savefig(results_folder + 'PQS' + image_type, transparent=True, bbox_inches='tight', dpi = dpi)

#Backward-only probability density plot
plt.figure()
forward_fig = plt.contourf(t_PQS, delta_q, P_b, clevels)
for c in forward_fig.collections:
    c.set_edgecolor("face")
plt.plot(t_emission, true_state, color=True_field_line_color, label='True field')
plt.ylabel('$\Delta_n / \gamma$')
plt.xlabel('$\gamma t$')
plt.yticks()
cbar = plt.colorbar()
cbar.ax.tick_params()
plt.legend(loc="upper right")
plt.savefig(results_folder + 'backward' + image_type, transparent=True, bbox_inches='tight', dpi = dpi)

#B_true - PQS_max histogram data
PQS_dif_list = true_state[index_frames_in_emission] - P_max
PQS_heights = np.zeros_like(delta_q)
for i in range(len(delta_q)):
        PQS_heights[i] = (np.around(PQS_dif_list, 1) == np.around(delta_q[i], 1)).sum()*(T / frames)


#B_true - PQS_forward historgram data
forward_dif_list = true_state[index_frames_in_emission] - P_f_max
forward_heights = np.zeros_like(delta_q)
for i in range(len(delta_q)):
        forward_heights[i] = (np.around(forward_dif_list, 2) == np.around(delta_q[i], 2)).sum()*(T / frames)

ss_heights = np.zeros_like(delta_q)
for i in range(len(delta_q)):
        ss_heights[i] = (true_state[index_frames_in_emission] == delta_q[i]).sum()*(T / frames)

#PQS histogram plot
plt.figure()
plt.grid()
plt.bar(delta_q, ss_heights, width=8/(2 * N - 1), color = 'xkcd:orange', alpha=0.4, label='Steady state MLE' )
plt.bar(delta_q, PQS_heights, width=8/(2 * N - 1), alpha=0.6, label='PQS MLE')
plt.xlabel('$\Delta_{n,True} - \Delta_{n,PQS MLE}$')
plt.ylabel('Time at $\Delta_{n,True} - \Delta_{n,MLE}$ in units of $\gamma t$')
plt.ylim((0, max(np.max(PQS_heights)*1.02 , np.max(forward_heights)*1.02)))
plt.yticks()
plt.xticks()
plt.legend()
plt.savefig(results_folder + 'PQS_histogram' + image_type, transparent=True, bbox_inches='tight', )

#Forward histogram plot
plt.figure()
plt.grid()
plt.bar(delta_q, ss_heights, width=8/(2 * N - 1), color = 'xkcd:orange', alpha=0.4, label='Steady state MLE' )
plt.bar(delta_q, forward_heights, width=8/(2 * N - 1), alpha=0.6, label='Filtered state MLE')
plt.xlabel('$\Delta_{n,True} - \Delta_{n,filtered MLE}$', fontsize=18)
plt.ylabel('Time at $\Delta_{n,True} - \Delta_{n,MLE}$ in units of $\gamma t$', fontsize=18)
plt.yticks()
plt.xticks()
plt.ylim((0, max(np.max(PQS_heights)*1.02 , np.max(forward_heights)*1.02)))
plt.legend()
plt.savefig(results_folder + 'forward_histogram' + image_type, transparent=True, bbox_inches='tight')


plt.figure()
plt.plot(t_PQS, PQS_weighted_variances, label='PQS variance')
plt.plot(t_PQS, forward_weighted_variances, label='filtering variance')
plt.legend(loc="upper right")
plt.title('P(t) variances as functions of time')
plt.savefig(results_folder + 'variance_through_time.pdf')
