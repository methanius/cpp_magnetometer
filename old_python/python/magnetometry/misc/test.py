import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def f(v, T):
	return v**3 * (np.exp(v / T) - 1)**-2

v = np.linspace(0, 10, 0.001)
y = f(v, 1)

plt.figure()
plt.plot(v,y)
plt.show()