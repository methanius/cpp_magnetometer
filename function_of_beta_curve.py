import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy as scp
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})
plt.rc("axes", labelsize=18, titlesize=22)   # skriftstørrelse af `xlabel`, `ylabel` og `title`
plt.rc("xtick", labelsize=16, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=16, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=16) # skriftstørrelse af figurers legends

def true_state_mean_deltas(true_state, data, weights):
    return true_state - np.average(data, weights = weights)


def weighted_variance(weights, data):
    weighted_variance = (weights * (data - np.average(data, weights = weights))**2).sum() / weights.sum()
    return weighted_variance

def MSE_errorbars(deltas):
    second_moment = (deltas**2).mean()
    fourth_moment = (deltas**4).mean()
    N = len(deltas)
    return (1 / N * (fourth_moment - (N - 3)/(N - 1) * second_moment**2))**0.5

def mean_and_error_of_the_mean(array):
    # Using the standard sample error of the mean
    return array.mean(), array.std(ddof = 1) / len(array)**0.5


beta_path = os.getcwd() + '/beta_images/'

# Initialize data arrays
N_betas = len(os.listdir(beta_path))
betas = np.zeros(N_betas)

PQS_MSEs = np.zeros(N_betas)
PQS_averaged_variance = np.zeros(N_betas)
PQS_MSE_errorbars = np.zeros(N_betas)
#PQS_MSE_errorbars_ = np.zeros(N_betas)
PQS_averaged_variance_errorbar = np.zeros(N_betas)

filtering_MSEs = np.zeros(N_betas)
filtering_averaged_variance = np.zeros(N_betas)
filtering_MSE_errorbars = np.zeros(N_betas)
#filtering_MSE_errorbars_2 = np.zeros(N_betas)
filtering_averaged_variance_errorbar = np.zeros(N_betas)

# Loading common data only once
common_data = beta_path + os.listdir(beta_path)[0]
common_input = common_data + '/input_parameters/'
common_results = common_data+ '/data/'

t_emission = np.loadtxt(common_results + 't_emission.txt')
t_PQS = np.loadtxt(common_results + 't_PQS.txt')
index_frames_in_emission = np.isin(t_emission, t_PQS)
delta_q = np.loadtxt(common_input + 'delta_q.txt')

beta_directories = os.scandir(beta_path)
for beta_index, beta_dir in enumerate(beta_directories):
    betas[beta_index] = beta_dir.name[-3:]

    print('%%%%%%%%%%  ' + beta_dir.name + '  %%%%%%%%%%' + '\n')

    data_path = beta_path + beta_dir.name + '/'
    betas[beta_index] = np.loadtxt(data_path  + 'input_parameters/' + 'beta.txt')
    true_state = np.array([delta_q[state] for state in \
            np.loadtxt(data_path + 'data/' + 'true_state.txt', dtype = int)])
    true_state = true_state[index_frames_in_emission]
    P_f_max = np.loadtxt(data_path + 'data/' + 'P_f_max.txt').reshape(-1)
    P_max = np.loadtxt(data_path + 'data/' + 'P_max.txt').reshape(-1)

    PQS_weights = np.loadtxt(data_path + 'data/' + 'P.txt')
    filtering_weights = np.loadtxt(data_path + 'data/' + 'P_f.txt')

    for i in range(len(P_max)):
        P_max[i] = delta_q[int(P_max[i])]
        P_f_max[i] = delta_q[int(P_f_max[i])]


    run_PQS_deltas = []
    run_PQS_variances = []
    run_filtering_deltas = []
    run_filtering_variances = []

    for n in range(len(true_state)):
        #run_PQS_deltas.append(true_state_mean_deltas(true_state[n], delta_q, PQS_weights[:, n]))
        run_PQS_deltas.append(true_state[n] - P_max[n])
        run_PQS_variances.append(weighted_variance(PQS_weights[:, n], delta_q))

        #run_filtering_deltas.append(true_state_mean_deltas(true_state[n], delta_q, filtering_weights[:, n]))
        run_filtering_deltas.append(true_state[n] - P_f_max[n])
        run_filtering_variances.append(weighted_variance(filtering_weights[:, n], delta_q))

    run_PQS_deltas = np.array(run_PQS_deltas)
    run_PQS_variances = np.array(run_PQS_variances)
    run_filtering_deltas = np.array(run_filtering_deltas)
    run_filtering_variances = np.array(run_filtering_variances)

    # Filling data arrays for each beta
    PQS_MSEs[beta_index], PQS_MSE_errorbars[beta_index] = mean_and_error_of_the_mean(run_PQS_deltas**2)
    PQS_averaged_variance[beta_index], PQS_averaged_variance_errorbar[beta_index] = mean_and_error_of_the_mean(run_PQS_variances)
    PQS_averaged_variance_errorbar[beta_index] = MSE_errorbars(run_PQS_deltas)

    filtering_MSEs[beta_index], filtering_MSE_errorbars[beta_index] = mean_and_error_of_the_mean(run_filtering_deltas**2)
    filtering_averaged_variance[beta_index], filtering_averaged_variance_errorbar[beta_index] = mean_and_error_of_the_mean(run_filtering_variances)
    filtering_averaged_variance_errorbar[beta_index] = MSE_errorbars(run_filtering_deltas)





sorting = betas.argsort()
betas = betas[sorting]

PQS_MSEs = PQS_MSEs[sorting]**0.5
PQS_averaged_variance = PQS_averaged_variance[sorting]**0.5
PQS_MSE_errorbars = PQS_MSE_errorbars[sorting]**0.5
PQS_averaged_variance_errorbar = PQS_averaged_variance_errorbar[sorting]**0.5

filtering_MSEs = filtering_MSEs[sorting]**0.5
filtering_averaged_variance = filtering_averaged_variance[sorting]**0.5
filtering_MSE_errorbars = filtering_MSE_errorbars[sorting]**0.5
filtering_averaged_variance_errorbar = filtering_averaged_variance_errorbar[sorting]**0.5

fig = plt.figure()
ax = fig.add_subplot()
#plt.errorbar(betas + 0.05, PQS_MSEs, PQS_MSE_errorbars_2, fmt = 'none')
#plt.errorbar(betas + 0.05, filtering_MSEs, filtering_MSE_errorbars_2, fmt = 'none')
ax.plot(betas, PQS_averaged_variance, '--k', label='PQS width', lw = 0.6)
ax.plot(betas, filtering_averaged_variance, '-.b', label='Filtering width', lw = 0.6)
ax.plot(betas, PQS_MSEs, 'sk', label = 'PQS RMSE')
ax.plot(betas, filtering_MSEs, '^b', label = 'Filtering RMSE')
"""
ax.fill_between(betas, PQS_averaged_variance - PQS_averaged_variance_errorbar, PQS_averaged_variance + PQS_averaged_variance_errorbar, facecolor = 'dimgrey', alpha=0.5)
ax.fill_between(betas, filtering_averaged_variance - filtering_averaged_variance_errorbar, filtering_averaged_variance + filtering_averaged_variance_errorbar, facecolor = 'dimgrey', alpha=0.5)
ax.errorbar(betas, PQS_MSEs, PQS_MSE_errorbars, fmt = 'sk', capsize = 1.3, elinewidth = 0.4, ecolor = 'black', label = 'PQS RMSE')
ax.errorbar(betas, filtering_MSEs, filtering_MSE_errorbars, fmt = '^b', capsize = 1.3, elinewidth = 0.4, ecolor = 'black', label = 'Filtering RMSE')
"""
ax.legend(loc = 'lower right', frameon = False)
ax.set_xlabel(r'Coherent drive strength $\beta$ (in units of $\sqrt{\gamma})$', )
ax.set_ylabel(r'Error estimates $\sqrt{(\Delta_n^2 / \gamma^2)}$')
ax.set_ylim(bottom = 0)
plt.savefig('article/drive_strength_figure.pdf', bbox_inches='tight')
plt.show()

