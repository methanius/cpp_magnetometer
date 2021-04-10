import numpy as np
import matplotlib.pyplot as plt
import os
import sys

path = sys.argv[1]

n_runs = len(os.listdir(path))
averaged_runs = np.arange(n_runs) + 1
PQS_mean_sum_of_prev = 0
PQS_mean_vs_true_SME_averaging = np.zeros(n_runs)
PQS_variance_averaging = np.zeros(n_runs)
filtering_mean_sum_of_prev = 0
filtering_mean_vs_true_SME_averaging = np.zeros(n_runs)
filtering_variance_averaging = np.zeros(n_runs)


# Loading common data only once

common_data_path = path + os.listdir(path)[0]

t_emission = np.loadtxt(common_data_path + '/data/t_emission.txt')
t_PQS = np.loadtxt(common_data_path + '/data/t_PQS.txt')
center_of_run = int(len(t_PQS) / 2)
index_frames_in_emission = np.isin(t_emission, t_PQS)
delta_q = np.loadtxt(common_data_path + '/input_parameters/delta_q.txt')

run_index = 0
for run in os.scandir(path):

    print('Analyzing #{0} out of {1} runs.'.format(run_index + 1, n_runs))

    run_P = np.loadtxt(path + run.name + '/data/' + 'P.txt')
    run_Pf = np.loadtxt(path + run.name + '/data/' + 'P_f.txt')
    true_state = np.loadtxt(path + run.name + '/data/true_state.txt')[index_frames_in_emission]
    true_state_center = delta_q[int(true_state[center_of_run])]

    #PQS SME calculations
    PQS_mean_at_center = (run_P[:, center_of_run] * delta_q).sum()
    PQS_mean_vs_true_square = (PQS_mean_at_center - true_state_center)**2
    PQS_mean_sum_of_prev += PQS_mean_vs_true_square
    PQS_mean_vs_true_SME_averaging[run_index] = PQS_mean_sum_of_prev / (run_index + 1)

    # PQS variance calculation
    PQS_variance = (run_P[:, center_of_run] * (delta_q - PQS_mean_at_center)**2).sum()
    if (run_index) != 0:
        PQS_variance_averaging[run_index] = 1 / (run_index + 1) * (PQS_variance_averaging[run_index - 1] * (run_index) + PQS_variance)
    else:
        PQS_variance_averaging[run_index] = PQS_variance


    #filtering SME calculations
    filtering_mean_at_center = (run_Pf[:, center_of_run] * delta_q).sum()
    filtering_mean_vs_true_square = (filtering_mean_at_center - true_state_center)**2
    filtering_mean_sum_of_prev += filtering_mean_vs_true_square
    filtering_mean_vs_true_SME_averaging[run_index] = (filtering_mean_sum_of_prev / (run_index + 1))

    # filtering variance calculation
    filtering_variance = (run_Pf[:, center_of_run] * (delta_q - filtering_mean_at_center)**2).sum()
    if (run_index) != 0:
        filtering_variance_averaging[run_index] = 1 / (run_index + 1) * (filtering_variance_averaging[run_index - 1] * (run_index) + filtering_variance)
    else:
        filtering_variance_averaging[run_index] = filtering_variance

    #Increment run index
    run_index += 1

plt.figure()
plt.title('Convergence of parameters with n runs')
plt.plot(averaged_runs, filtering_mean_vs_true_SME_averaging, '-.', label='Filtering SME')
plt.plot(averaged_runs, filtering_variance_averaging,'-.', label='Filtering variance')
plt.plot(averaged_runs, PQS_mean_vs_true_SME_averaging, '--', label='PQS SME')
plt.plot(averaged_runs, PQS_variance_averaging,'--', label='PQS variance')
plt.legend()
plt.xlabel('# of runs')
plt.show()
