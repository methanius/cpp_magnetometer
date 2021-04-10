import numpy as np
from scipy import sparse

def spost(A):
    return sparse.csr_matrix(sparse.kron(sparse.identity(A.shape[0]), A.T))
