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

def SME_errorbars(deltas):
    second_moment = (deltas**2).mean()
    fourth_moment = (deltas**4).mean()
    N = len(deltas)
    return (1 / N * (fourth_moment - (N - 3)/(N - 1) * second_moment**2))**0.5

subsets = 20


beta_path = os.getcwd() + '/beta_images/'

# Initialize data arrays
N_betas = len(os.listdir(beta_path))
betas = np.zeros(N_betas)

PQS_SMEs = np.zeros(N_betas)
PQS_averaged_variances = np.zeros(N_betas)
PQS_SME_errorbars = np.zeros(N_betas)
#PQS_SME_errorbars_ = np.zeros(N_betas)
PQS_averaged_variance_errorbars = np.zeros(N_betas)

filtering_SMEs = np.zeros(N_betas)
filtering_averaged_variances = np.zeros(N_betas)
filtering_SME_errorbars = np.zeros(N_betas)
#filtering_SME_errorbars_2 = np.zeros(N_betas)
filtering_averaged_variance_errorbars = np.zeros(N_betas)

# Loading common data only once
first_directory_for_common_data = beta_path + os.listdir(beta_path)[0]
common_data_path = first_directory_for_common_data + '/' + os.listdir(first_directory_for_common_data)[0]
common_input = common_data_path + '/input_parameters/'
common_results = common_data_path + '/data/'

t_emission = np.loadtxt(common_results + 't_emission.txt')
t_PQS = np.loadtxt(common_results + 't_PQS.txt')
index_frames_in_emission = np.isin(t_emission, t_PQS)
delta_q = np.loadtxt(common_input + 'delta_q.txt')
run_center = int(len(t_PQS) / 2)

# Below line finds the middle common index of t_emission and t_PQS
true_state_center = np.where(index_frames_in_emission.cumsum() == index_frames_in_emission.sum()/2)[0][0]

# Data analysis loop filling in the above arrays
beta_directories = os.scandir(beta_path)
for beta_index, beta_dir in enumerate(beta_directories):
    betas[beta_index] = beta_dir.name[-3:]

    run_directories = os.scandir(beta_dir)
    print('\n\n' + '%%%%%%%%%%  ' + beta_dir.name + '  %%%%%%%%%%' + '\n')
    N_runs = len(os.listdir(beta_dir))
    if (N_runs % subsets != 0):
        print('The number of subsets chosen must divide the number of runs cleanly.')
        print('Currently: subsets = {0} and N_runs = {1}.'.format(subsets, N_runs))
        sys.exit()
    run_PQS_deltas = np.zeros(N_runs)
    run_PQS_variances = np.zeros(N_runs)
    run_filtering_deltas = np.zeros(N_runs)
    run_filtering_variances = np.zeros(N_runs)

    for n, run_dir in enumerate(run_directories):
        print('- ' + run_dir.name +': -')

        data_path = beta_path + beta_dir.name + '/' + run_dir.name + '/'
        betas[beta_index] = np.loadtxt(data_path  + 'input_parameters/' + 'beta.txt')
        true_state = delta_q[int(np.loadtxt(data_path + 'data/' + 'true_state.txt', skiprows = true_state_center, max_rows = 1))]

        PQS_weights = np.loadtxt(data_path + 'data/' + 'P.txt', usecols = run_center)
        filtering_weights = np.loadtxt(data_path + 'data/' + 'P_f.txt', usecols = run_center)

        run_PQS_deltas[n] = true_state_mean_deltas(true_state, delta_q, PQS_weights)
        run_PQS_variances[n] = weighted_variance(PQS_weights, delta_q)

        run_filtering_deltas[n] = true_state_mean_deltas(true_state, delta_q, filtering_weights)
        run_filtering_variances[n] = weighted_variance(filtering_weights, delta_q)

    # Filling data arrays for each beta
    PQS_SMEs[beta_index] = (run_PQS_deltas**2).mean()
    PQS_averaged_variances[beta_index] = run_PQS_variances.mean()

    filtering_SMEs[beta_index] = (run_filtering_deltas**2).mean()
    filtering_averaged_variances[beta_index] = run_filtering_variances.mean()

    PQS_SME_errorbars[beta_index] = SME_errorbars(run_PQS_deltas)
#    PQS_SME_errorbars_2[beta_index] = (2 * PQS_SMEs[beta_index]**2 / (N_runs - 1))**0.5

    PQS_averaged_variance_errorbars[beta_index] = run_PQS_variances.std(ddof=1) / (N_runs)**0.5

    filtering_SME_errorbars[beta_index] = SME_errorbars(run_filtering_deltas)
#    filtering_SME_errorbars_2[beta_index] = (2 * filtering_SMEs[beta_index]**2 / (N_runs - 1))**0.5

    filtering_averaged_variance_errorbars[beta_index] = run_filtering_variances.std(ddof=1) / (N_runs)**0.5




sorting = betas.argsort()
betas = betas[sorting]
PQS_SMEs = PQS_SMEs[sorting]
PQS_averaged_variances = PQS_averaged_variances[sorting]
PQS_SME_errorbars = PQS_SME_errorbars[sorting]
#PQS_SME_errorbars_2 = PQS_SME_errorbars_2[sorting]
PQS_averaged_variance_errorbars = PQS_averaged_variance_errorbars[sorting]

filtering_SMEs = filtering_SMEs[sorting]
filtering_averaged_variances = filtering_averaged_variances[sorting]
filtering_SME_errorbars = filtering_SME_errorbars[sorting]
#filtering_SME_errorbars_2 = filtering_SME_errorbars_2[sorting]
filtering_averaged_variance_errorbars = filtering_averaged_variance_errorbars[sorting]

fig = plt.figure()
ax = fig.add_subplot()
#plt.errorbar(betas + 0.05, PQS_SMEs, PQS_SME_errorbars_2, fmt = 'none')
#plt.errorbar(betas + 0.05, filtering_SMEs, filtering_SME_errorbars_2, fmt = 'none')
ax.plot(betas, PQS_averaged_variances, '--k', label='PQS Bayes\' width', lw = 0.6)
ax.fill_between(betas, PQS_averaged_variances - PQS_averaged_variance_errorbars, PQS_averaged_variances + PQS_averaged_variance_errorbars, facecolor = 'dimgrey', alpha=0.5)
ax.plot(betas, filtering_averaged_variances, '-.b', label='Filtering Bayes\' width', lw = 0.6)
ax.fill_between(betas, filtering_averaged_variances - filtering_averaged_variance_errorbars, filtering_averaged_variances + filtering_averaged_variance_errorbars, facecolor = 'dimgrey', alpha=0.5)
ax.errorbar(betas, PQS_SMEs, PQS_SME_errorbars, fmt = 'sk', capsize = 1.3, elinewidth = 0.4, ecolor = 'black', label = 'PQS SME')
ax.errorbar(betas, filtering_SMEs, filtering_SME_errorbars, fmt = '^b', capsize = 1.3, elinewidth = 0.4, ecolor = 'black', label = 'Filtering SME')
ax.legend(loc = 'upper center', frameon = False)
ax.set_xlabel(r'Coherent drive strength $(\beta / \sqrt{\gamma})$', )
ax.set_ylabel(r'Error estimates $(\Delta_n^2 / \gamma^2)$')
ax.set_ylim(bottom = 0)
plt.savefig('article/drive_strength_figure.pdf', bbox_inches='tight')
plt.show()


