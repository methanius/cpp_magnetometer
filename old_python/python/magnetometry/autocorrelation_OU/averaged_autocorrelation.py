import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt

B_autocorrelations = []
PQS_autocorrelations = []
t_autocorrelations = []

for i in range(1, len(sys.argv)):
	with open(sys.argv[i], 'rb') as f:
		[t_temp, B_temp, PQS_temp] = pickle.load(f)
		B_autocorrelations.append(B_temp)
		PQS_autocorrelations.append(PQS_temp)
		t_autocorrelations.append(t_temp)

average_B = np.zeros_like(B_autocorrelations[0])
average_PQS = np.zeros_like(PQS_autocorrelations[0])
average_t = np.zeros_like(t_autocorrelations[0])

for i in range(len(B_autocorrelations)):
	average_B += B_autocorrelations[i]
	average_PQS += PQS_autocorrelations[i]
	average_t += t_autocorrelations[i]

average_B /= len(B_autocorrelations)
average_PQS /= len(PQS_autocorrelations)
average_t /= len(t_autocorrelations)



plt.figure()
plt.plot(average_t, average_B, label='$\Delta_{n,True}$')
plt.plot(average_t, average_PQS, label='$\Delta_{n,PQS MLE}$')
plt.grid()
plt.xlabel(r'$\gamma_{dec} \tau$', fontsize=18)
plt.ylabel(r'$R(\gamma_{dec} \tau)$', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.savefig('averaged_aucorrelation.pdf', transparent = 'True', bbox_inches = 'tight')
