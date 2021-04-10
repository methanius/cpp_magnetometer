import numpy as np
from scipy import sparse

def make_operators(size, g, kappa, delta_r, delta_s, kappa_1, beta, gamma_dec, gamma_phi, LOPhi, W_up, W_down):
    sigmam = np.array([[0, 0], [1, 0]])
    sigmaz = np.array([[1, 0], [0, -1]])
    gamma_p = np.zeros(size)
    alpha = np.zeros(size, dtype=complex)
    a_ad = np.zeros((size, 2, 2), dtype=complex)
    eps_s = np.zeros(size)
    H_ap = np.zeros((size, 2, 2))
    for i in range(size):
        gamma_p[i] = 2 * g[i]**2 * kappa[i] / (kappa[i]**2 + (delta_r[i] - delta_s[i])**2)
        alpha[i] = (2 * kappa_1[i])**0.5 * beta[i] / (kappa[i] + 1j * delta_r[i])
        a_ad[i] = alpha[i] * np.identity(2) - 1j * g[i] * sigmam / (kappa[i] + 1j * (delta_r[i] - delta_s[i]))
        eps_s[i] = (delta_r[i] - delta_s[i]) * g[i]**2 * kappa[i] / (kappa[i]**2 + (delta_r[i] - delta_s[i])**2)
        H_ap[i] = (delta_s[i] * sigmaz / 2 + g[i] * (alpha[i] * sigmam.T + alpha[i].conj() * sigmam) - eps_s[i] * sigmam.T @ sigmam).real

    #Block diagonal Hamiltonian and operators
    c_1 = np.zeros((2 * size, 2 * size))
    c_2 = np.zeros((2 * size, 2 * size))
    c_3 = np.zeros((2 * size, 2 * size))
    c = np.zeros((2 * size, 2 * size), dtype=complex)
    H = np.zeros((2 * size, 2 * size))
    Jup = np.zeros((2 * size, 2 * size))
    Jdown = np.zeros((2 * size, 2 * size))
    for i in range(size):
        subrho = np.zeros((size,size))
        subrho[size - i - 1, size - i - 1] = 1
        H += np.kron(subrho, H_ap[i])
        c_1 += np.kron(subrho, gamma_p[i]**0.5 * sigmam)
        c_2 +=  np.kron(subrho, gamma_dec[i]**0.5 * sigmam)
        c_3 += np.kron(subrho, (gamma_phi[i] / 2)**0.5 * sigmaz)
        c += np.kron(subrho, (2 * kappa_1[i])**0.5 * (a_ad[i] - beta[i] * np.identity(2)) * np.exp(-1j * LOPhi))
        if i != size - 1:
            up_element = np.zeros((size, size))
            down_element = np.zeros((size, size))
            up_element[i, i + 1] = 1
            down_element[i + 1, i] = 1
            Jup += np.kron(up_element, W_up[i]**0.5 * np.identity(2))
            Jdown += np.kron(down_element, W_down[i]**0.5 * np.identity(2))
    H = sparse.csr_matrix(H)
    c = sparse.csr_matrix(c)
    c_1 = sparse.csr_matrix(c_1)
    c_2 = sparse.csr_matrix(c_2)
    c_3 = sparse.csr_matrix(c_3)
    Jup = sparse.csr_matrix(Jup)
    Jdown = sparse.csr_matrix(Jdown)
    return c_1, c_2, c_3, c, H, Jup, Jdown
