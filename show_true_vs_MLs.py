import numpy as np
import matplotlib.pyplot as plt
import sys
import os


image_type = '.pdf'
target_directory = os.getcwd() + '/' + sys.argv[1]
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
P_max = np.loadtxt(in_data + 'P_max.txt')#.reshape(-1)
P_f_max = np.loadtxt(in_data + 'P_f_max.txt')#.reshape(-1)

for i in range(len(true_state)):
    true_state[i] = delta_q[int(true_state[i])]

for i in range(len(P_max)):
    P_max[i] = delta_q[int(P_max[i])]
    P_f_max[i] = delta_q[int(P_f_max[i])]

print(delta_q)

plt.figure()
plt.plot(t_PQS, P_max)
plt.plot(t_emission, true_state)
plt.show()
