import numpy as np
import matplotlib.pyplot as plt


scatterdata = np.loadtxt('scatterplot_data.txt')
bins = int(np.loadtxt('gridnumber.txt'))
dt = np.loadtxt('dt.txt')
true_state = scatterdata[:, 0]
P_max = scatterdata[:, 1]
P_f_max = scatterdata[:, 2]



plt.figure(1)
plt.hist(true_state - (P_max), bins=np.linspace(-2, 2, bins), weights=dt*np.ones_like(true_state))
plt.xlabel('$B_{true} - B_{PQS ML}$')
plt.ylabel('Time spend at dif')
plt.title('$B_{true} - B_{PQS ML}$ occupation through time')
plt.savefig('PQS_histogram.png')


plt.figure(2)
plt.hist(true_state - P_f_max, bins=np.linspace(-2, 2, bins), weights=dt*np.ones_like(true_state))
plt.xlabel('$B_{true} - B_{forward ML}$')
plt.ylabel('Time spend at dif')
plt.title('$B_{true} - B_{forward ML}$ occupation through time')
plt.savefig('forward_histogram.png')
