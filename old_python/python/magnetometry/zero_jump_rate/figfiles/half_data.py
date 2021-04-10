import sys
import numpy as np


data = np.loadtxt(sys.argv[1])

half_name = sys.argv[1][:-4] + '_original.txt'
np.savetxt(half_name, data)
np.savetxt(sys.argv[1], data[:, ::2])
