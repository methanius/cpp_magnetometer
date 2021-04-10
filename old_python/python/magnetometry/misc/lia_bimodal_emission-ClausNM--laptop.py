# This code models a bimodal atom with coherent photon emissions and incoherent jumps between the two states
from typing import Any, Union

import numpy as np
from qutip import *
from hmmlearn import hmm
from numpy.random import seed
from numpy.random import rand
import matplotlib.pyplot as plt
from matplotlib import colors


def rev_vec(v):
    u = []
    for i in range(0, len(v)):
        u.append(0)

    j = len(u) - 1
    for i in range(0, len(v)):
        u[j] = v[i]
        j -= 1

    return u


def update_PQS(time, dt, omega, delta, gamma, eta, W, sd, dN):
    aa = Qobj([[1, 0], [0, 0]])
    bb = Qobj([[0, 0], [0, 1]])
    ab = Qobj([[0, 1], [0, 0]])
    ba = Qobj([[0, 0], [1, 0]])

    sm = tensor(qeye(2), ba)
    sp = tensor(qeye(2), ab)

    rho = tensor(fock_dm(2, 1), aa) + tensor(fock_dm(2, 0), bb)

    H = tensor(omega[0] / 2 * sigmax() + delta[0] / 2 * sigmaz(), aa) + \
        tensor(omega[1] / 2 * sigmax() + delta[1] / 2 * sigmaz(), bb)
    c = np.sqrt(gamma[0]) * tensor(sigmam(), aa) + np.sqrt(gamma[1]) * tensor(sigmam(), bb)

    Jab = np.sqrt(W[0]) * tensor(qeye(2), ba)  # incoherent jumping rates
    Jba = np.sqrt(W[1]) * tensor(qeye(2), ab)

    # Liovillian with coherent photon emission and incoherent jump evolution
    L = -1j * (spre(H) - spost(H)) + 1 / 2 * (2 * sprepost(Jab.dag(), Jab) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (2 * sprepost(Jba.dag(), Jba) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) - \
        1 / 2 * (spre(c.dag() * c) + spost(c.dag() * c))
    ldt = (L * dt).expm()

    seed(sd)
    n = 0
    Pa = []
    Pb = []
    rho_save = []
    for t in time:
        if dN[n] == 0:
            rho = vector_to_operator(ldt * operator_to_vector(rho))
        else:
            rho = c * rho * c.dag()
            rho = rho / rho.tr()

        rho_save.append(rho)
        P1 = np.abs((sp * rho * sp.dag()).tr())
        P2 = np.abs((sm * rho * sm.dag()).tr())
        norm = P1 + P2
        Pa.append(P1 / norm)
        Pb.append(P2 / norm)
        n += 1

    E = tensor(fock_dm(2, 0), aa) + tensor(fock_dm(2, 0), bb)  # initiate E = I
    nb = len(time) - 1
    Pa_back = []
    Pb_back = []
    Lb = 1j * (spre(H) - spost(H)) + 1 / 2 * (2 * sprepost(Jab, Jab.dag()) - spre(Jab * Jab.dag())
                                              - spost(Jab * Jab.dag())) + \
         1 / 2 * (2 * sprepost(Jba, Jba.dag()) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) - \
         1 / 2 * (spre(c.dag() * c) + spost(c.dag() * c))
    lbdt = (Lb * dt).expm()
    #    calculates the effect matrix
    for t in time:
        if dN[nb] == 0:
            E = vector_to_operator(lbdt * operator_to_vector(E))
        else:
            E = c.dag() * E * c
            E = E / E.tr()

        P1 = np.abs((sp * rho_save[nb] * E * sp.dag()).tr())
        P2 = np.abs((sm * rho_save[nb] * E * sm.dag()).tr())
        norm = P1 + P2
        Pa_back.append(P1 / norm)
        Pb_back.append(P2 / norm)
        nb -= 1

    Pa_back = rev_vec(Pa_back)
    Pb_back = rev_vec(Pb_back)

    return Pa, Pb, Pa_back, Pb_back


def update(time, dt, omega, delta, gamma, eta, W, sd, dN):
    aa = Qobj([[1, 0], [0, 0]])
    bb = Qobj([[0, 0], [0, 1]])
    ab = Qobj([[0, 1], [0, 0]])
    ba = Qobj([[0, 0], [1, 0]])

    sm = tensor(qeye(2), ba)
    sp = tensor(qeye(2), ab)

    rho = tensor(fock_dm(2, 1), aa) + tensor(fock_dm(2, 0), bb)

    H = tensor(omega[0] / 2 * sigmax() + delta[0] / 2 * sigmaz(), aa) + \
        tensor(omega[1] / 2 * sigmax() + delta[1] / 2 * sigmaz(), bb)
    c = np.sqrt(gamma[0]) * tensor(sigmam(), aa) + np.sqrt(gamma[1]) * tensor(sigmam(), bb)

    Jab = np.sqrt(W[0]) * tensor(qeye(2), ba)  # incoherent jumping rates
    Jba = np.sqrt(W[1]) * tensor(qeye(2), ab)

    # Liovillian with coherent photon emission and incoherent jump evolution
    L = -1j * (spre(H) - spost(H)) + 1 / 2 * (2 * sprepost(Jab.dag(), Jab) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (2 * sprepost(Jba.dag(), Jba) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) - \
        1 / 2 * (spre(c.dag() * c) + spost(c.dag() * c))
    ldt = (L * dt).expm()

    seed(sd)
    n = 0
    Pa = []
    Pb = []
    for t in time:
        if dN[n] == 0:
            rho = vector_to_operator(ldt * operator_to_vector(rho))
        else:
            rho = c.dag() * rho * c
            rho = rho / rho.tr()

        P1 = np.abs((sp.dag() * rho * sp).tr())
        P2 = np.abs((sm.dag() * rho * sm).tr())
        norm = P1 + P2
        print(P1 / norm)
        Pa.append(P1 / norm)
        Pb.append(P2 / norm)
        n += 1

    return Pa, Pb


def non_obs(time, dt, omega, delta, gamma, eta, W, sd):
    aa = Qobj([[1, 0], [0, 0]])
    bb = Qobj([[0, 0], [0, 1]])
    ab = Qobj([[0, 1], [0, 0]])
    ba = Qobj([[0, 0], [1, 0]])

    rho = tensor(fock_dm(2, 1), aa) + tensor(fock_dm(2, 0), bb)

    H = tensor(omega[0] / 2 * sigmax() + delta[0] / 2 * sigmaz(), aa) + \
        tensor(omega[1] / 2 * sigmax() + delta[1] / 2 * sigmaz(), bb)
    c = np.sqrt(gamma[0]) * tensor(sigmam(), aa) + np.sqrt(gamma[1]) * tensor(sigmam(), bb)

    Jab = np.sqrt(W[0]) * tensor(qeye(2), ba)  # incoherent jumping rates
    Jba = np.sqrt(W[1]) * tensor(qeye(2), ab)

    sm = tensor(qeye(2), ba)
    sp = tensor(qeye(2), ab)

    L = -1j * (spre(H) - spost(H)) + 1 / 2 * (2 * sprepost(Jab.dag(), Jab) - spre(Jab * Jab.dag()) +
                                              spost(Jab * Jab.dag())) \
        + 1 / 2 * (2 * sprepost(Jba.dag(), Jba) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c))
    ldt = (L * dt).expm()

    a_exp = []
    b_exp = []
    for t in time:
        rho = vector_to_operator(ldt * operator_to_vector(rho))
        rho = rho.unit()
        a_exp.append(np.abs((sp.dag() * rho * sp).tr()))
        b_exp.append(np.abs((sm.dag() * rho * sm).tr()))

    return a_exp, b_exp


def emission(time, dt, omega, delta, gamma, eta, W, sd):
    aa = Qobj([[1, 0], [0, 0]])
    bb = Qobj([[0, 0], [0, 1]])
    ab = Qobj([[0, 1], [0, 0]])
    ba = Qobj([[0, 0], [1, 0]])

    rho = tensor(fock_dm(2, 1), aa) + tensor(fock_dm(2, 0), bb)

    H = tensor(omega[0] / 2 * sigmax() + delta[0] / 2 * sigmaz(), aa) + \
        tensor(omega[1] / 2 * sigmax() + delta[1] / 2 * sigmaz(), bb)
    c = np.sqrt(gamma[0]) * tensor(sigmam(), aa) + np.sqrt(gamma[1]) * tensor(sigmam(), bb)

    Jab = np.sqrt(W[0]) * tensor(qeye(2), ba)  # incoherent jumping rates
    Jba = np.sqrt(W[1]) * tensor(qeye(2), ab)

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) - 1 / 2 * (spre(Jab * Jab.dag()) + spost(Jab * Jab.dag())) \
        - 1 / 2 * (spre(Jba * Jba.dag()) + spost(Jba * Jba.dag())) - \
        1 / 2 * (spre(c.dag() * c) + spost(c.dag() * c))

    ldt = (L * dt).expm()

    sm = tensor(qeye(2), ba)
    sp = tensor(qeye(2), ab)
    dN = []
    dNtime = []
    dNtime.append(0)
    seed(sd)
    a_exp = []
    b_exp = []
    change = []
    for t in time:
        trans = rand()
        epsilon = rand()

        if trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
            rho = (Jab * rho * Jab.dag())
            change.append(t)
            dN.append(0)
            if np.abs(rho.tr()) > 0:
                rho = rho / rho.tr()

        elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
            rho = (Jba * rho * Jba.dag())
            change.append(t)
            dN.append(0)
            if np.abs(rho.tr()) > 0:
                rho = rho / rho.tr()

        elif epsilon < np.abs((c.dag() * c * rho).tr() * dt):
            rho = (c * rho * c.dag())
            dN.append(1)
            dNtime.append(t)
            if np.abs(rho.tr()) > 0:
                rho = rho / rho.tr()

        else:
            rho = vector_to_operator(ldt * operator_to_vector(rho))
            rho = rho.unit()
            dN.append(0)

        a = np.abs((sp * rho * sp.dag()).tr())
        b = np.abs((sm * rho * sm.dag()).tr())
        norm = a + b
        a_exp.append(a / norm)
        b_exp.append(b / norm)

    return dN, a_exp, b_exp, change, dNtime


def grid_bayes_twolevel(time, dt, gamma, delta_list, eta, omega, dN):
    c = np.sqrt(gamma) * sigmam()  # decay operator

    rho_test = []
    ldt_test = []
    follow = []
    for om in delta_list:
        rho_test.append(fock_dm(2, 0))
        H = (1 / 2) * (omega * sigmax() + om * sigmaz())  # test hamiltonian
        L = -1j * (spre(H) - spost(H)) - 1 / 2 * (spre(c.dag() * c) + spost(c.dag() * c))  # liouvilian

        ldt_test.append((L * dt).expm())
        follow.append(True)

    Prob_matrix = np.ones((len(delta_list), len(time)))
    n = 0
    for t in time:
        trans = rand()
        # ----------- bayesian update ----------------------
        i = 0
        norm = 0
        ext = 0
        for om in delta_list:
            if follow[i] is True:
                if dN[n] == 0:
                    rho_test[i] = vector_to_operator(ldt_test[i] * operator_to_vector(rho_test[i]))
                else:
                    rho_test[i] = c * rho_test[i] * c.dag()

                P = np.abs(rho_test[i].tr())

                if np.isnan(P) is True or P < 1e-5:
                    follow[i] = False
                    P = 0
                    Prob_matrix[i, n] = P
                else:
                    follow[i] = True
                    Prob_matrix[i, n] = P
                    norm += P
            else:
                Prob_matrix[i, n] = 0

            i += 1

        j = 0
        check_sum = 0
        for om in delta_list:
            Prob_matrix[j, n] = Prob_matrix[j, n] / (norm + ext)
            check_sum = check_sum + Prob_matrix[j, n]
            j += 1

        n += 1

    best_curve = []
    for t in range(0, len(time)):
        print(t)
        P = Prob_matrix[0, t]
        omega_guess = 1
        i = 0
        for om in delta_list:
            if Prob_matrix[i, t] < P:
                P = Prob_matrix[i, t]
                omega_guess = om
            i += 1
        best_curve.append(omega_guess)

    return Prob_matrix, best_curve


def rev_vec(v):
    u = []
    for i in range(0, len(v)):
        u.append(0)

    j = len(u) - 1
    for i in range(0, len(v)):
        u[j] = v[i]
        j -= 1

    return u


def off_diag(i):
    S = []
    for j in range(0, i):
        M = np.zeros((i, i))
        M[j, i - 1 - j] = 1
        S.append(Qobj(M))
    return S


def J_matrix(i):
    J_for = []
    J_bac = []
    for j in range(0, i):
        M = np.zeros((i, i))
        N = np.zeros((i, i))
        if 0 < j < i - 1:
            M[j, j + 1] = 1
            N[j, j - 1] = 1
        elif j == i - 1:
            N[j, j - 1] = 1
        elif j == 0:
            M[j, j + 1] = 1

        J_for.append(Qobj(M))
        J_bac.append(Qobj(N))

    return J_for, J_bac


def twolevel_emission(time, dt, omega, delta, gamma, eta, W_f, W_b, sd):
    print(int(len(delta) / 2))
    H_ap = []
    for i in range(0, len(delta)):
        H_ap.append(delta[i] * sigmaz() / 2 + omega[i] * sigmax())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c = 0
    rho = 0
    H = 0
    rho2 = 0
    for i in range(0, len(delta)):
        H = H + tensor(H_ap[i], fock_dm(len(delta), i))
        c = c + tensor(np.sqrt(gamma[i]) * sigmam(), fock_dm(len(delta), i))
        if i == int(len(delta) / 2):
            rho = rho + tensor(fock_dm(2, 0), fock_dm(len(delta), i))
            rho2 = rho2 + tensor(fock_dm(2, 0), fock_dm(len(delta), i))
        else:
            rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta), i))
            rho2 = rho2 + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta))
    for i in range(0, len(delta), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (- spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (- spre(Jba * Jba.dag()) - spost(Jba * Jba.dag()))

    ldt = (L * dt).expm()

    L2 = -1j * (spre(H) - spost(H)) + \
         1 / 2 * (2 * sprepost(Jab, Jab.dag()) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
         1 / 2 * (2 * sprepost(Jba, Jba.dag()) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag()))

    ldt2 = (L * dt).expm()

    seed(sd)
    a_exp = []
    b_exp = []
    dW = []  # generates random noise
    dY = []  # homodyne current
    n = 0
    s = off_diag(len(delta))
    P_list = np.zeros((len(delta), len(time)))
    follow_traj = []
    state_list = []
    state_list2 = []
    for t in time:
        print(t)
        trans = rand()
        test = np.random.randint(0, high=100, size=None)
        if test > 50:
            w = 0
        else:
            w = 1
        if w == 0:
            if trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

            elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
        elif w == 1:
            if trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
            elif trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

        rho = vector_to_operator(ldt * operator_to_vector(rho))
        rho = rho.unit()

        rho2 = vector_to_operator(ldt2 * operator_to_vector(rho2))
        rho2 = rho2.unit()

        norm = 0
        P = []
        for i in range(0, len(delta)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        test = 0
        follow = 0
        for i in range(0, len(delta)):
            P_list[i, n] = P[i]
            if test < P[i]:
                test = P[i]
                follow = i

        follow_traj.append(delta[follow])

        for j in range(0, len(delta)):
            x = 2 * np.real(rho[j, j + len(delta)])
            y = - 2 * np.imag(rho[j, j + len(delta)])
            z = np.real(rho[j, j] - rho[j + len(delta), j + len(delta)])

            if x != 0 or y != 0 or z != 0:
                state = [x, y, z]

        state_list.append(state)

        norm = 0
        P2 = []
        for i in range(0, len(delta)):
            sm = tensor(qeye(2), s[i])
            P2.append(np.abs((sm * rho2 * sm.dag()).tr()))
            norm += np.abs((sm * rho2 * sm.dag()).tr())

        P2 = P2 / norm
        test = 0
        follow = 0
        for i in range(0, len(delta)):
            if test < P2[i]:
                test = P2[i]
                follow = i

        for j in range(0, len(delta)):
            x = 2 * np.real(rho2[j, j + len(delta)])
            y = - 2 * np.imag(rho2[j, j + len(delta)])
            z = np.real(rho2[j, j] - rho2[j + len(delta), j + len(delta)])

            if x != 0 or y != 0 or z != 0:
                state = [x, y, z]

        state_list2.append(state)

        n += 1

    return dY, P_list, follow_traj, state_list, state_list2


def twolevel_emission2(time, dt, omega, delta, gamma, eta, W_f, W_b, sd):
    print(int(len(delta) / 2))
    H_ap = []
    for i in range(0, len(delta)):
        H_ap.append(delta[i] * sigmaz() / 2 + omega[i] * sigmax())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c = 0
    rho = 0
    H = 0
    rho2 = 0
    for i in range(0, len(delta)):
        H = H + tensor(H_ap[i], fock_dm(len(delta), i))
        c = c + tensor(np.sqrt(gamma[i]) * sigmam(), fock_dm(len(delta), i))
        if i == int(len(delta) / 2):
            rho = rho + tensor(fock_dm(2, 0), fock_dm(len(delta), i))
            rho2 = rho2 + tensor(fock_dm(2, 0), fock_dm(len(delta), i))
        else:
            rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta), i))
            rho2 = rho2 + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta))
    for i in range(0, len(delta), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (- spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (- spre(Jba * Jba.dag()) - spost(Jba * Jba.dag()))

    ldt = (L * dt).expm()

    seed(sd)
    a_exp = []
    b_exp = []
    dW = []  # generates random noise
    dY = []  # homodyne current
    n = 0
    s = off_diag(len(delta))
    P_list = np.zeros((len(delta), len(time)))
    follow_traj = []
    state_list = []
    state_list2 = []
    for t in time:
        print(t)
        trans = rand()
        test = np.random.randint(0, high=100, size=None)
        if test > 50:
            w = 0
        else:
            w = 1
        if w == 0:
            if trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

            elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
        elif w == 1:
            if trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
            elif trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

        rho = vector_to_operator(ldt * operator_to_vector(rho))
        rho = rho.unit()

        norm = 0
        P = []
        for i in range(0, len(delta)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        test = 0
        follow = 0
        for i in range(0, len(delta)):
            P_list[i, n] = P[i]
            if test < P[i]:
                test = P[i]
                follow = i

        follow_traj.append(delta[follow])

        for j in range(0, len(delta)):
            x = 2 * np.real(rho[j, j + len(delta)])
            y = - 2 * np.imag(rho[j, j + len(delta)])
            z = np.real(rho[j, j] - rho[j + len(delta), j + len(delta)])

            if x != 0 or y != 0 or z != 0:
                state = [x, y, z]

        state_list.append(state)

        n += 1

    return dY, P_list, follow_traj, state_list


def twolevel_PQS(time, dt, omega, delta, gamma, eta, W_f, W_b, sd, dY):
    H_ap = []
    for i in range(0, len(delta)):
        H_ap.append(delta[i] * sigmaz() / 2 + omega[i] * sigmax())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c = 0
    rho = 0
    H = 0
    for i in range(0, len(delta)):
        H = H + tensor(H_ap[i], fock_dm(len(delta), i))
        c = c + tensor(np.sqrt(gamma[i]) * sigmam(), fock_dm(len(delta), i))
        rho = rho + tensor(fock_dm(2, 0), fock_dm(len(delta), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta))
    for i in range(0, len(delta), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (2 * sprepost(Jab.dag(), Jab) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (2 * sprepost(Jba.dag(), Jba) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c))

    ldt = (L * dt).expm()

    seed(sd)
    n = 0
    Pa = []
    Pb = []
    rho_save = []
    s = off_diag(len(delta))
    P_list = np.zeros((len(delta), len(time)))
    follow_traj = []
    for t in time:
        rho = vector_to_operator(ldt * operator_to_vector(rho)) + (c * rho + rho * c.dag()) * dY[n] * np.sqrt(eta)
        rho = rho.unit()
        rho_save.append(rho)

        norm = 0
        P = []
        for i in range(0, len(delta)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        for i in range(0, len(delta)):
            P_list[i, n] = P[i]
        n += 1

        test = 0
        follow = 0
        for i in range(0, len(delta)):
            if test < P[i]:
                test = P[i]
                follow = i
        follow_traj.append(delta[follow])

    # Propagating the past quantum state
    E = 0
    for i in range(0, len(delta)):
        E = E + tensor(qeye(2), fock_dm(len(delta), i))

    Pb_list = np.zeros((len(delta), len(time)))
    Lb = 1j * (spre(H) - spost(H)) + \
         1 / 2 * (2 * sprepost(Jab, Jab.dag()) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
         1 / 2 * (2 * sprepost(Jba, Jba.dag()) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
         1 / 2 * (2 * sprepost(c, c.dag()) - spre(c.dag() * c) - spost(c.dag() * c))

    lbdt = (Lb * dt).expm()
    #    calculates the effect matrix
    nb = len(time) - 1
    follow_trajb = []
    for t in time:
        E = vector_to_operator(lbdt * operator_to_vector(E)) + (c.dag() * E + E * c) * dY[nb] * np.sqrt(eta)
        # E = E.unit()
        P = []
        norm = 0
        for i in range(0, len(delta)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho_save[nb] * E * sm.dag()).tr()))
            norm += np.abs((sm * rho_save[nb] * E * sm.dag()).tr())
            # print(norm)

        if np.isnan(norm):
            print('ahh, Im breaking')
            print(norm)
            # Propagating the past quantum state
            E = 0
            norm = 0
            for i in range(0, len(delta)):
                E = E + tensor(qeye(2), fock_dm(len(delta), i)) * Pb_list[i, nb + 1]
                P[i] = Pb_list[i, nb + 1]
                print(Pb_list[i, nb + 1])
                norm = norm + P[i]
            print(norm)
        elif norm > 0:
            P = P / norm

        for i in range(0, len(delta)):
            Pb_list[i, nb] = P[i]

        test = 0
        follow = 0
        for i in range(0, len(delta)):
            if test < P[i]:
                test = P[i]
                follow = i
        follow_trajb.append(delta[follow])

        nb -= 1

    follow_trajb = rev_vec(follow_trajb)

    return P_list, Pb_list, follow_traj, follow_trajb


def twolevel_emission_dY(time, dt, omega, delta, gamma, eta, W_f, W_b, sd):
    print(int(len(delta) / 2))
    H_ap = []
    for i in range(0, len(delta)):
        H_ap.append(delta[i] * sigmaz() / 2 + omega[i] * sigmax())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c = 0
    rho = 0
    H = 0
    for i in range(0, len(delta)):
        H = H + tensor(H_ap[i], fock_dm(len(delta), i))
        c = c + tensor(np.sqrt(gamma[i]) * sigmam(), fock_dm(len(delta), i))
        if i == int(len(delta) / 2):
            rho = rho + tensor(fock_dm(2, 0), fock_dm(len(delta), i))
        else:
            rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta))
    for i in range(0, len(delta), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c)) + \
        1 / 2 * (- spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (- spre(Jba * Jba.dag()) - spost(Jba * Jba.dag()))

    ldt = (L * dt).expm()

    seed(sd)
    a_exp = []
    b_exp = []
    dW = []  # generates random noise
    dY = []  # homodyne current
    n = 0
    s = off_diag(len(delta))
    P_list = np.zeros((len(delta), len(time)))
    follow_traj = []
    state_list = []
    state_list2 = []
    for t in time:
        print(t)
        trans = rand()
        test = np.random.randint(0, high=100, size=None)
        if test > 50:
            w = 0
        else:
            w = 1
        if w == 0:
            if trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

            elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
        elif w == 1:
            if trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
            elif trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                # print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

        dW = np.random.normal(loc=0.0, scale=1.0, size=None) * np.sqrt(dt)
        dY.append((c * rho + rho * c.dag()).tr() * np.sqrt(eta) * dt + dW)

        rho = vector_to_operator(ldt * operator_to_vector(rho)) + (c * rho + rho * c.dag()) * dY[n] * np.sqrt(eta)
        rho = rho.unit()

        norm = 0
        P = []
        for i in range(0, len(delta)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        test = 0
        follow = 0
        for i in range(0, len(delta)):
            P_list[i, n] = P[i]
            if test < P[i]:
                test = P[i]
                follow = i

        follow_traj.append(delta[follow])

        for j in range(0, len(delta)):
            x = 2 * np.real(rho[j, j + len(delta)])
            y = - 2 * np.imag(rho[j, j + len(delta)])
            z = np.real(rho[j, j] - rho[j + len(delta), j + len(delta)])

            if x != 0 or y != 0 or z != 0:
                state = [x, y, z]

        state_list.append(state)

        n += 1

    return dY, follow_traj, state_list


def main():
    t_max = 50
    time = np.linspace(0, t_max, 5000)
    dt = time[1] - time[0]
    '''
    W = [.1, .1]
    omega = [1, 1]
    delta = [1.5, .5]
    gamma = [4, .1]
    eta = 1.
    sd = 11

    dN, P_reala, P_realb, change, dNtime = emission(time, dt, omega, delta, gamma, eta, W, sd)
    Pa, Pb, Pa_pass, Pb_pass = update_PQS(time, dt, omega, delta, gamma, eta, W, sd, dN)
    #rho_a, rho_b = non_obs(time, dt, omega, delta, gamma, eta, W, sd)


    fig2 = plt.figure(2)
    plt.title('photon counting signal')
    plt.plot(time, dN)
    plt.xlabel('$t \Omega $')
    plt.gca().axes.get_yaxis().set_visible(False)
    fig2.savefig('bimodal_dN.eps')

    fig1 = plt.figure(1)
    plt.ylim((0, 1))
    plt.title('$Tr(\\rho_{1})$')
    plt.fill_between(time, 0, Pb, color='b')
    plt.fill_between(time, 0, P_realb, color='grey', alpha=.5)
    plt.ylabel('$P_1$')
    plt.gca().axes.get_yaxis().set_visible(True)
    plt.xlabel('$t \\kappa$')
    fig2.savefig('bimodal_P1.eps')

    fig3 = plt.figure(3)
    plt.ylim((0, 1))
    plt.title('$Tr(\\rho_{1})$')
    plt.fill_between(time, 0, Pb_pass, color='b')
    plt.fill_between(time, 0, P_realb, color='grey', alpha=.5)
    plt.ylabel('$P_1$')
    plt.xlabel('$t \Omega $')
    fig3.savefig('bimodal_P1_PQS.eps')

    bin = 100
    dN_bin = []
    time_bin = []
    for i in range(0, len(dN), bin):
        count = 0
        for j in range(0, bin):
            if dN[i+j] == 1:
                count += 1
        dN_bin.append(count)
        time_bin.append(i * dt)


    fig4 =  plt.figure(4)
    plt.plot(time_bin, dN_bin)
    plt.xlabel('$t \Omega $')
    '''
    # ---- PC bayes estimate in a grid-----------------------------

    '''
    W = [0, 0]
    omega = [.5, .5]
    delta = [.5, .5]
    gamma = [1, 1]
    eta = 1.
    sd = 11

    dN, P_reala, P_realb, change, dNtime = emission(time, dt, omega, delta, gamma, eta, W, sd)

    gamma = gamma[0]
    delta = delta[0]
    delta_list = np.linspace(0, 2, 21)
    omega = omega[0]
    P, best_curve = grid_bayes_twolevel(time, dt, gamma, delta_list, eta, omega, dN)
    tlist = []
    omlist = []
    clist = []
    i = 0

    fig2 = plt.figure(2)
    plt.title('photon counting signal')
    plt.plot(time, dN)
    plt.xlabel('$t \Omega $')
    plt.gca().axes.get_yaxis().set_visible(False)

    fig6 = plt.figure(6)
    plt.contourf(time, delta_list, P, levels=len(delta_list)-1)
    plt.colorbar()
    '''

    # ---- many trajectories for Ehrenfest detuning ---------------

    size = 21
    delta = 0.5 * np.linspace(-1, 1, size)
    omega = 0.1 * np.ones(size)
    gamma = 0.1 * np.ones(size)

    eta = 1.
    sd = 11

    W_f1 = []
    W_b1 = []
    seed(sd)
    size2 = int(size)
    for i in range(0, size2):
        # ehrenfest chain
        s = 0.995
        #W_f1.append((i / (size2 - 1)) * (1 - s))
        #W_b1.append((1. - i / (size2 - 1)) * (1 - s))
        W_f1.append(0)
        W_b1.append(0)

    dY, follow_traj, state_list = twolevel_emission_dY(time, dt, omega, delta, gamma, eta, W_f1, W_b1, sd)
    P_list, Pb_list, follow_traj, follow_trajb = twolevel_PQS(time, dt, omega, delta, gamma, eta, W_f1, W_b1, sd, dY)




    '''
    W_f = []
    W_b = []
    for i in range(0, size):
        W_f.append(W_f1[i + size2 - size])
        W_b.append(W_b1[i])
    print(W_f)
    print(W_b)
    '''
    '''
    dY, P_list, follow_traj, state_list = twolevel_emission_dY(time, dt, omega,
                                                               delta, gamma, eta, W_f1, W_b1, sd)

    P_list, Pb_list, follow_traj, follow_trajb = twolevel_PQS(time, dt, omega, delta, gamma, eta, W_f1, W_b1, sd, dY)

    fig1 = plt.figure(1)
    plt.contourf(time, delta, P_list, levels=100)
    plt.plot(time, follow_traj, color='r')
    plt.colorbar()

    fig2 = plt.figure(2)
    plt.contourf(time, delta, Pb_list, levels=100)
    plt.plot(time, follow_traj, color='r')
    plt.colorbar()

    exp_sx = []
    exp_sy = []
    exp_sz = []

    for t in range(0, len(time)):
        state = state_list[t]
        norm = state[0] ** 2 + state[1] ** 2 + state[2] ** 2

        exp_sx.append(state[0] / np.sqrt(norm))
        exp_sy.append(state[1] / np.sqrt(norm))
        exp_sz.append(state[2] / np.sqrt(norm))

    fig3 = plt.figure(3)
    plt.plot(time, exp_sz)
    plt.plot(time, follow_traj)
    '''
    '''
    dY, P_list, follow_traj, state_list, state_list2 = twolevel_emission(time, dt, omega,
                                                                         delta, gamma, eta, W_f, W_b, sd)

    state1 = []
    for i in range(0, 10):
        sd = i
        dY, P_list, follow_traj, state_list = twolevel_emission2(time, dt, omega, delta, gamma, eta, W_f, W_b, sd)

        state1.append(state_list)

    exp_sx = []
    exp_sy = []
    exp_sz = []

    exp_sx2 = []
    exp_sy2 = []
    exp_sz2 = []

    for t in range(0, len(time)):
        state = state_list[t]
        norm = state[0]**2 + state[1]**2 + state[2]**2

        exp_sx.append(state[0]/np.sqrt(norm))
        exp_sy.append(state[1]/np.sqrt(norm))
        exp_sz.append(state[2]/np.sqrt(norm))

        state = state_list2[t]
        norm = state[0] ** 2 + state[1] ** 2 + state[2] ** 2

        exp_sx2.append(state[0] / np.sqrt(norm))
        exp_sy2.append(state[1] / np.sqrt(norm))
        exp_sz2.append(state[2] / np.sqrt(norm))



    b = Bloch()
    b.view = [-40, 30]
    b.add_points([exp_sx[:len(time)], exp_sy[:len(time)], exp_sz[:len(time)]], meth='l')
    b.save(dirc='temp')

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    plt.rcParams.update({'font.size': 20})

    fig4 = plt.figure(4)
    plt.plot(time, follow_traj, color='r')
    plt.xlabel('$t/\gamma$')
    plt.ylabel('$\Delta/0.1\gamma$')

    fig5 = plt.figure(5)
    ax3 = plt.subplot(313)
    plt.ylabel('z')
    plt.plot(time, exp_sz, alpha=0.4)
    plt.plot(time, exp_sz2, color='r')
    for i in range(0, 10):
        exp_sx3 = []
        exp_sy3 = []
        exp_sz3 = []
        state_time = state1[i]
        for t in range(0, len(time)):
            state = state_time[t]
            norm = state[0] ** 2 + state[1] ** 2 + state[2] ** 2
            exp_sx3.append(state[0] / np.sqrt(norm))
            exp_sy3.append(state[1] / np.sqrt(norm))
            exp_sz3.append(state[2] / np.sqrt(norm))
        plt.plot(time, exp_sz3, alpha=0.4)

    plt.xlabel('$t/\gamma$')
    ax1 = plt.subplot(311, sharex=ax3)
    plt.ylabel('x')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot(time, exp_sx, alpha=0.4)
    plt.plot(time, exp_sx2, color='r')
    for i in range(0, 10):
        exp_sx3 = []
        exp_sy3 = []
        exp_sz3 = []
        state_time = state1[i]
        for t in range(0, len(time)):
            state = state_time[t]
            norm = state[0] ** 2 + state[1] ** 2 + state[2] ** 2
            exp_sx3.append(state[0] / np.sqrt(norm))
            exp_sy3.append(state[1] / np.sqrt(norm))
            exp_sz3.append(state[2] / np.sqrt(norm))
        plt.plot(time, exp_sx3, alpha=0.4)
    ax2 = plt.subplot(312, sharex=ax3)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.ylabel('y')
    plt.plot(time, exp_sy, alpha=0.4)
    plt.plot(time, exp_sy2, color='r')
    for i in range(0, 10):
        exp_sx3 = []
        exp_sy3 = []
        exp_sz3 = []
        state_time = state1[i]
        for t in range(0, len(time)):
            state = state_time[t]
            norm = state[0] ** 2 + state[1] ** 2 + state[2] ** 2
            exp_sx3.append(state[0] / np.sqrt(norm))
            exp_sy3.append(state[1] / np.sqrt(norm))
            exp_sz3.append(state[2] / np.sqrt(norm))
        plt.plot(time, exp_sy3, alpha=0.4)
        
        
    '''
    '''
    # -------- OU exact solution -------------------------------#
    sigma = 1
    mu = 0
    theta = 0.05
    sigma = 1
    x = []
    omega = 1
    delta = 0
    phi = fock(2, 0)
    H = omega / 2 * sigmax() + delta / 2 * sigmaz()
    state_save = []
    for t in range(0, len(time)):
        print(t)
        if t == 0:
            x.append(0)
        elif t > 0:
            x.append(- theta * (x[t-1] - mu) * dt + sigma * np.random.randn() * np.sqrt(dt))
        L = - 1j * H + x[t] * sigmaz()
        ldt = (L*dt).expm()
        phi = ldt * phi
        phi = phi.unit()
        state_save.append(np.asscalar(phi[0]))

    fig1 = plt.figure(1)
    plt.plot(time, state_save)
    plt.plot(time, x)
    '''
    plt.show()


if __name__ == '__main__':
    main()
