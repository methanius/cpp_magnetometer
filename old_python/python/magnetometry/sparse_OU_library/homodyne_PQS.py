import numpy as np
import scipy.linalg as spl
from scipy import sparse
from sparse_OU_library.choi_jamiolkowsky.spre import spre
from sparse_OU_library.choi_jamiolkowsky.spost import spost
from sparse_OU_library.choi_jamiolkowsky.sprepost import sprepost


def homodyne_PQS(nT, size, H, c, Jdown, Jup, c_1, c_2, c_3, dt, dY, eta, delta_s):
    rho = sparse.lil_matrix((nT, 4 * size**2), dtype=complex)
    rho[0] = sparse.lil_matrix(sparse.identity(2 * size).reshape(1, 4 * size**2) / (2 * size))
    E = sparse.lil_matrix(rho[::-1])

    # Two Liouvillians
    L = -1j * (spre(H) - spost(H)) \
        + sprepost(Jdown, Jdown.conj().T) - (spre(Jdown.conj().T * Jdown) + spost(Jdown.conj().T * Jdown)) / 2 \
        + sprepost(Jup, Jup.conj().T) - (spre(Jup.conj().T * Jup) + spost(Jup.conj().T * Jup)) / 2 \
        + sprepost(c_1, c_1.conj().T) - (spre(c_1.conj().T * c_1) + spost(c_1.conj().T * c_1)) / 2 \
        + sprepost(c_2, c_2.conj().T) - (spre(c_2.conj().T * c_2) + spost(c_2.conj().T * c_2)) / 2 \
        + sprepost(c_3, c_3.conj().T) - (spre(c_3.conj().T * c_3) + spost(c_3.conj().T * c_3)) / 2
    Ldt = sparse.csr_matrix(spl.expm(sparse.csc_matrix(L) * dt))

    EL = 1j * (spre(H) - spost(H)) \
        + sprepost(Jdown.conj().T, Jdown) - (spre(Jdown.conj().T * Jdown) + spost(Jdown.conj().T * Jdown)) / 2 \
        + sprepost(Jup.conj().T, Jup) - (spre(Jup.conj().T * Jup) + spost(Jup.conj().T * Jup)) / 2 \
        + sprepost(c_1.conj().T, c_1) - (spre(c_1.conj().T * c_1) + spost(c_1.conj().T * c_1)) / 2 \
        + sprepost(c_2.conj().T, c_2) - (spre(c_2.conj().T * c_2) + spost(c_2.conj().T * c_2)) / 2 \
        + sprepost(c_3.conj().T, c_3) - (spre(c_3.conj().T * c_3) + spost(c_3.conj().T * c_3)) / 2
    ELdt = sparse.csr_matrix(spl.expm(sparse.csc_matrix(EL) * dt))

    # Propagating rho and E
    rho_vec = np.zeros((4 * size**2, 1))
    rho_vec = sparse.lil_matrix(rho_vec)
    E_vec = sparse.lil_matrix(rho_vec)
    c_pre = spre(c)
    c_H_pre = spre(c.conj().T)
    c_post = spost(c)
    c_H_post = spost(c.conj().T)

    for i in range(1, nT):
        rho_vec = sparse.csr_matrix(rho[i - 1, :]).T
        rho[i, :] = (Ldt * rho_vec).T + (c_pre * rho_vec + c_H_post * rho_vec).T * dY[i - 1] * eta**0.5
        rho[i, :] /= rho[i, ::(2 * size + 1)].sum()
        E_vec = sparse.csr_matrix(E[nT - i, :]).T
        E[nT - i - 1, :] = (ELdt * E_vec).T + (c_H_pre * E_vec + c_post * E_vec).T * dY[nT - i - 1] * eta**0.5
        E[nT - i - 1, :] /= E[nT - i - 1, ::(2 * size + 1)].sum()

    # Subtracing guess state amplitudes
    P = np.zeros((size, nT))
    P_f = np.zeros((size, nT))
    P_b = np.zeros((size, nT))
    for i in range(nT):

        # Rewrite of square coordinates to linear coordinate, row-wise: (i, j) -> j + len(j) * i
        for n in range(size):
            P[size - 1 - n, i] = (rho[i, 2 * n + 2 * size * 2 * n] * E[i, 2 * n + 2 * size * 2 * n] + 
                       rho[i, 2 * n + 1 + 2 * size * (2 * n + 1)] * E[i, 2 * n + 1 + 2 * size * (2 * n + 1)] +
                       rho[i, 2 * n + 2 * size * (2 * n + 1)] * E[i, 2 * n + 1 + 2 * size * 2 * n] +
                       rho[i, 2 * n + 1 + 2 * size * 2 * n] * E[i, 2 * n + 2 * size * (2 * n + 1)]).real
            P_f[size - 1 - n, i] = (rho[i, 2 * n + 2 * size * 2 * n] + rho[i, 2 * n + 1 + 2 * size * (2 * n + 1)]).real
            P_b[size - 1 - n, i] = (E[i, 2 * n + 2 * size * 2 * n] + E[i, 2 * n + 1 + 2 * size * (2 * n + 1)]).real        
    P[:, :] /= np.sum(P[:, :], axis= 0)

    P_max = np.zeros(nT, dtype=float)
    P_max[:] = np.argmax(P, 0)
    P_f_max = np.zeros(nT, dtype=float)
    P_f_max[:] = np.argmax(P_f, 0)
    for n in range(nT):
        P_max[n] = delta_s[int(P_max[n])]
        P_f_max[n] = delta_s[int(P_f_max[n])]
    return P, P_f, P_b, P_max, P_f_max
