import numpy as np
from scipy import sparse
from sparse_OU_library.choi_jamiolkowsky.spre import spre
from sparse_OU_library.choi_jamiolkowsky.spost import spost
from sparse_OU_library.choi_jamiolkowsky.sprepost import sprepost
from scipy.linalg import expm


def homodyne_emission(size, H, Jup, Jdown, c, c_1, c_2, c_3, dt, nT, eta, emission_startpoint, delta_s):
    data = np.array([1, 1])
    end = 2 * size - 1
    row = np.array([end - 2 * emission_startpoint, end - 2 * emission_startpoint - 1])
    col = np.array([end - 2 * emission_startpoint, end - 2 * emission_startpoint - 1])
    rho = sparse.csr_matrix((data, (row, col)), shape=(2 * size, 2 * size))

    # Liouvillian
    L = -1j * (spre(H) - spost(H)) \
        - (spre(Jup.conj().T * Jup) + spost(Jup.conj().T * Jup)) / 2 \
        - (spre(Jdown.conj().T * Jdown) + spost(Jdown.conj().T * Jdown)) / 2 \
        + sprepost(c_1, c_1.conj().T) - (spre(c_1.conj().T * c_1) + spost(c_1.conj().T * c_1)) / 2 \
        + sprepost(c_2, c_2.conj().T) - (spre(c_2.conj().T * c_2) + spost(c_2.conj().T * c_2)) / 2 \
        + sprepost(c_3, c_3.conj().T) - (spre(c_3.conj().T * c_3) + spost(c_3.conj().T * c_3)) / 2

    Ldt = sparse.csr_matrix(expm(sparse.csc_matrix(L) * dt))
    Ldt.data = np.round(Ldt.data, 16)
    print(Ldt)

    # Signal generation
    dW = np.random.normal(0, 1, nT) * dt**0.5
    dY = np.zeros(nT)
    state_track = emission_startpoint
    true_state = np.zeros(nT)
    for i in range(nT):
        jump1 = np.random.rand()
        test_permutation = np.random.rand()
        if test_permutation <= 0.5:
            if jump1 < (Jup * rho * Jup.conj().T).diagonal().sum() * dt:
                rho = Jup @ rho @ Jup.conj().T
                rho /= rho.diagonal().sum()
                state_track += 1
            else:
                jump2 = np.random.rand()
                if jump2 < (Jdown * rho * Jdown.conj().T).diagonal().sum() * dt:
                    rho = Jdown * rho * Jdown.conj().T
                    rho /= rho.diagonal().sum()
                    state_track -= 1
        else:
            if jump1 < (Jdown * rho * Jdown.conj().T).diagonal().sum() * dt:
                rho = Jdown * rho * Jdown.conj().T
                rho /= rho.diagonal().sum()
                state_track -= 1
            else:
                jump2 = np.random.rand()
                if jump2 < (Jup * rho * Jup.conj().T).diagonal().sum() * dt:
                    rho = Jup * rho * Jup.conj().T
                    rho /= rho.diagonal().sum()
                    state_track += 1
        dY[i] = ((c * rho + rho * c.conj().T).diagonal().sum() * eta**0.5 * dt + dW[i]).real
        rho_vec = sparse.csr_matrix(rho.reshape(4 * size ** 2, 1))
        rho = sparse.csr_matrix((Ldt * rho_vec).reshape(2 * size, 2 * size)) + (c * rho + rho * c.conj().T) * dY[i] * eta**0.5
        rho /= rho.diagonal().sum()
        true_state[i] = state_track
    # Return state of B-field in units of Delta_s
    for n in range(nT):
        true_state[n] = delta_s[int(true_state[n])]
    return dY, true_state
