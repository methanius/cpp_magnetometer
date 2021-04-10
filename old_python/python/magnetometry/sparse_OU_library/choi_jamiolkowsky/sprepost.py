import numpy as np
from scipy import sparse


def sprepost(A, B):
    return sparse.csr_matrix(sparse.kron(A, B.T))
