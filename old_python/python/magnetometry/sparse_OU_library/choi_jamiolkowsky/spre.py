import numpy as np
from scipy import sparse

def spre(A):
    return sparse.csr_matrix(sparse.kron(A, sparse.identity(A.shape[0])))
