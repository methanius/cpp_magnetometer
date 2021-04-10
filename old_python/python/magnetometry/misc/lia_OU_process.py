'''
import numpy as np
from numpy.random import seed
from numpy.random import rand
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits import axes_grid1
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def transition_matrix(M):
    forw = []
    bac = []
    for i in range(0, len(M) - 1):
        forw.append(M[i + 1]/M[i])

    for i in range(len(M) - 1, 0, -1):
        bac.append(M[i - 1]/M[i])

    return forw, bac

def main():
    th = 1
    mu = 1.2
    sig = 1
    t_max = 10000
    time = np.linspace(0, t_max, 30000)
    dt = time[1] - time[0]

    x = []
    x4 = []
    x3 = []
    x.append(1)
    x4.append(0)
    x3.append(3)
    sd = 3
    s = []
    seed(sd)
    events = 0
    
    for i in range(0, len(time) - 1):
        sd = np.random.normal(loc=0.0, scale=dt, size=None)
        x.append(x[i] + th * (mu - x[i]) * dt + sig * sd)
        events += 1
        if i > len(time) / 2:
            s.append(mu - x[i] + th * (mu - x[i]) * dt + sig * sd)

    fig1 = plt.figure(1)
    n, bins, patches = plt.hist(s, 40, density=True)

    hist, bin_edges = np.histogram(s, 40, density=False)
    hist = hist/events
    
    ## Let us check the Ehrenfest chain
    N = 100
    W_b = []
    W_f = []
    for i in range(0, N):
        # ehrenfest chain
        s = 0.9
        W_b.append((i / (N - 1)) * (1 - s))
        W_f.append((1. - i / (N - 1)) * (1 - s))

    x_eh = []
    whereami = 58
    for t in range(0, len(time)):
        trans = rand()
        test = np.random.randint(0, high=100, size=None)
        if test > 50:
            w = 0
        else:
            w = 1
        if w == 0:
            if trans < W_f[whereami] * dt:
                whereami = whereami + 1

            elif trans < W_b[whereami] * dt:
                whereami = whereami - 1
        elif w == 1:
            if trans < W_b[whereami] * dt:
                whereami = whereami - 1
            elif trans < W_f[whereami] * dt:
                whereami = whereami + 1

        x_eh.append(whereami)
        print(whereami)

    bin_list = []
    for i in range(35, 66):
        bin_list.append(i)
        bin_list.append(i+1)
    fig3 = plt.figure(3)
    n, bins, patches = plt.hist(x_eh, bins=bin_list, density=True)
    plt.title('position in the chain')

    fig4 = plt.figure(4)
    plt.plot(time, x_eh)

    plt.show()
    
    for i in range(0, len(time) - 1):
        sd = np.random.normal(loc=0.0, scale=dt, size=None)
        x.append(x[i] + th * (mu - x[i]) * dt + sig * sd)
        sd = np.random.normal(loc=0.0, scale=dt, size=None)
        x4.append(x4[i] + th * (mu - x4[i]) * dt + sig * sd)
        sd = np.random.normal(loc=0.0, scale=dt, size=None)
        x3.append(x3[i] + th * (mu - x3[i]) * dt + sig * sd)
        if i > len(time)/2:
            s.append(mu - x[i] + th * (mu - x[i]) * dt + sig * sd)
            s.append(mu - x4[i] + th * (mu - x4[i]) * dt + sig * sd)
            s.append(mu - x3[i] + th * (mu - x3[i]) * dt + sig * sd)

    fig1, ax = plt.subplots()
    plt.plot(time, x4)
    plt.plot(time, x)
    plt.plot(time, x3)
    plt.xlabel('t')
    plt.ylabel('$\mu$')
    plt.axhline(y=mu, linestyle=':')
    axins = zoomed_inset_axes(ax, 5, loc=4)
    x1, x2, y1, y2 = 290, 320, 1.1, 1.3  # specify the limits
    axins.set_xlim(x1, x2)  # apply the x-limits
    axins.set_ylim(y1, y2)  # apply the y-limits
    axins.plot(time, x)
    axins.plot(time, x4)
    axins.plot(time, x3)
    axins.axhline(y=mu, linestyle=':')
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.5")
    #ax = fig1.add_subplot(121)
    inaxis = inset_axes(ax,
                            width="50%",  # width = 30% of parent_bbox
                            height=1.5,  # height : 1 inch
                            loc=2)
    plt.yticks(visible=False)
    n, bins, patches = plt.hist(s, 400, density=True)
    fig1.savefig('ou_path.eps')

    xra = []
    gaus = []
    xq = -2
    dx = - 4 * xq / 1000
    for i in range(0, 1000):
        xq += dx
        print(xq)
        xra.append(xq)
        y = np.sqrt(th / (np.pi * sig * sig)) * np.exp(-th * (xq - mu) * (xq - mu) / (sig * sig))
        gaus.append(y)

    fig2 = plt.figure(2)
    plt.plot(xra, gaus)
    plt.show()

if __name__ == '__main__':
    main()
    
'''
import numpy as np
from qutip import *
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


def homodyne_emission_simple(time, dt, delta_s, delta_r, g, gamma_dec, gamma_phi, kappa_1, kappa,
                             eta, sd, beta, Phi):
    alpha = np.sqrt(2 * kappa_1) * beta / (kappa + 1j * delta_r)
    eps_s = (delta_r - delta_s) * g ** 2 * kappa / (kappa ** 2 + (delta_r - delta_s) ** 2)
    gamma_p = 2 * g ** 2 * kappa / (kappa ** 2 + (delta_r - delta_s) ** 2)
    H = delta_s * sigmaz() / 2 + g * (alpha * sigmam().dag() + alpha.conjugate() * sigmam()) - \
        eps_s * sigmam().dag() * sigmam()
    a_ad = alpha - 1j * g * sigmaz() / (kappa + 1j * (delta_r - delta_s))

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    rho = (fock_dm(2, 1) + fock_dm(2, 0)).unit()
    c_1 = np.sqrt(gamma_p) * sigmam()
    c_2 = np.sqrt(gamma_dec) * sigmam()
    c_3 = np.sqrt(gamma_phi) * sigmam()
    c = np.sqrt(2 * kappa_1) * (a_ad - beta) * np.exp(-1j * Phi)

    Rabi = 1
    delta_s = 0

    H = delta_s * sigmax() / 2 + Rabi * sigmax()

    c = np.sqrt(2 * kappa_1) * sigmam() * np.exp(-1j * Phi)

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c)) + \
        1 / 2 * (2 * sprepost(c_1.dag(), c_1) - spre(c_1 * c_1.dag()) - spost(c_1 * c_1.dag())) + \
        1 / 2 * (2 * sprepost(c_2.dag(), c_2) - spre(c_2 * c_2.dag()) - spost(c_2 * c_2.dag())) + \
        1 / 2 * (2 * sprepost(c_3.dag(), c_3) - spre(c_3 * c_3.dag()) - spost(c_3 * c_3.dag()))

    ldt = (L * dt).expm()

    # sm = tensor(qeye(2), ba)
    # sp = tensor(qeye(2), ab)
    seed(sd)
    a_exp = []
    b_exp = []
    dW = []  # generates random noise
    dY = []  # homodyne current
    n = 0
    rho_list = []
    for t in time:
        print(t)
        dW.append(np.random.normal(loc=0.0, scale=1.0, size=None) * np.sqrt(dt))
        dY.append((c * rho + rho * c.dag()).tr() * np.sqrt(eta) * dt + dW[n])

        rho = vector_to_operator(ldt * operator_to_vector(rho)) + \
              (c * rho + rho * c.dag() - (c * rho + rho * c.dag()).tr() * rho) * dY[n] * np.sqrt(eta)
        rho = rho.unit()
        print(rho)
        rho_list.append(rho)
        n += 1

    return dY, rho_list


def homodyne_emission(time, dt, delta_s, delta_r, g, gamma_dec, gamma_phi, kappa_1, kappa,
                      eta, W_f, W_b, sd, beta, Phi):
    gamma_p = []  # represents the Purcell enhanced damping of the spin
    alpha = []
    a_ad = []
    eps_s = []
    H_ap = []
    for i in range(0, len(delta_s)):
        gamma_p.append(2 * g[i] ** 2 * kappa[i] / (kappa[i] ** 2 + (delta_r[i] - delta_s[i]) ** 2))
        alpha.append(np.sqrt(2 * kappa_1[i]) * beta[i] / (kappa[i] + 1j * delta_r[i]))
        a_ad.append(alpha[i] - 1j * g[i] * sigmaz() / (kappa[i] + 1j * (delta_r[i] - delta_s[i])))
        eps_s.append((delta_r[i] - delta_s[i]) * g[i] ** 2 * kappa[i] / (kappa[i] ** 2
                                                                         + (delta_r[i] - delta_s[i]) ** 2))
        H_ap.append(delta_s[i] * sigmaz() / 2 + g[i] * (
                alpha[i] * sigmam().dag() + alpha[i].conjugate() * sigmam()) - eps_s[i] * sigmam().dag() * sigmam())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c_1 = 0
    c_2 = 0
    c_3 = 0
    c = 0
    rho = 0
    H = 0
    for i in range(0, len(delta_s)):
        H = H + tensor(H_ap[i], fock_dm(len(delta_s), i))
        c_1 = c_1 + tensor(np.sqrt(gamma_p[i]) * sigmam(), fock_dm(len(delta_s), i))
        c_2 = c_2 + tensor(np.sqrt(gamma_dec[i]) * sigmam(), fock_dm(len(delta_s), i))
        c_3 = c_3 + tensor(np.sqrt(gamma_phi[i] / 2) * sigmam(), fock_dm(len(delta_s), i))
        c = c + tensor((np.sqrt(2 * kappa_1[i]) * a_ad[i] - beta[i]) * np.exp(-1j * Phi), fock_dm(len(delta_s), i))
        if i == int(len(delta_s) / 2):  # or i == int(3 * len(delta_s) / 2):
            rho = rho + tensor(fock_dm(2, 1), fock_dm(len(delta_s), i))
        else:
            rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta_s), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta_s))
    for i in range(0, len(delta_s), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c)) + \
        1 / 2 * (- spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (- spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
        1 / 2 * (2 * sprepost(c_1.dag(), c_1) - spre(c_1 * c_1.dag()) - spost(c_1 * c_1.dag())) + \
        1 / 2 * (2 * sprepost(c_2.dag(), c_2) - spre(c_2 * c_2.dag()) - spost(c_2 * c_2.dag())) + \
        1 / 2 * (2 * sprepost(c_3.dag(), c_3) - spre(c_3 * c_3.dag()) - spost(c_3 * c_3.dag()))

    ldt = (L * dt).expm()

    # sm = tensor(qeye(2), ba)
    # sp = tensor(qeye(2), ab)
    seed(sd)
    a_exp = []
    b_exp = []
    dW = []  # generates random noise
    dY = []  # homodyne current
    n = 0
    s = off_diag(len(delta_s))
    P_list = np.zeros((len(delta_s), len(time)))
    follow_traj = []
    rho_list = []
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
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

            elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
        elif w == 1:
            if trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
            elif trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

        dW.append(np.random.normal(loc=0.0, scale=1.0, size=None) * np.sqrt(dt))
        dY.append((c * rho + rho * c.dag()).tr() * np.sqrt(eta) * dt + dW[n])

        rho = vector_to_operator(ldt * operator_to_vector(rho)) + (c * rho + rho * c.dag()) * dY[n] * np.sqrt(eta)
        rho = rho / rho.tr()
        #rho_list.append(rho)  # GET THIS BACK FOR ESTIMATION ON THE BLOCH SPHERE
        norm = 0
        P = []
        for i in range(0, len(delta_s)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        test = 0
        follow = 0
        for i in range(0, len(delta_s)):
            P_list[i, n] = P[i]
            if test < P[i]:
                test = P[i]
                follow = i
        follow_traj.append(follow)
        n += 1

    return dY, P_list, follow_traj, rho_list
    # return dY, P_list, follow_traj


def homodyne_emission_steady(time, dt, delta_s, delta_r, g, gamma_dec, gamma_phi, kappa_1, kappa,
                       eta, W_f, W_b, sd, beta, Phi):
    gamma_p = []  # represents the Purcell enhanced damping of the spin
    alpha = []
    a_ad = []
    eps_s = []
    H_ap = []
    for i in range(0, len(delta_s)):
        gamma_p.append(2 * g[i] ** 2 * kappa[i] / (kappa[i] ** 2 + (delta_r[i] - delta_s[i]) ** 2))
        alpha.append(np.sqrt(2 * kappa_1[i]) * beta[i] / (kappa[i] + 1j * delta_r[i]))
        a_ad.append(alpha[i] - 1j * g[i] * sigmaz() / (kappa[i] + 1j * (delta_r[i] - delta_s[i])))
        eps_s.append((delta_r[i] - delta_s[i]) * g[i] ** 2 * kappa[i] / (kappa[i] ** 2
                                                                         + (delta_r[i] - delta_s[i]) ** 2))
        H_ap.append(delta_s[i] * sigmaz() / 2 + g[i] * (
                alpha[i] * sigmam().dag() + alpha[i].conjugate() * sigmam()) - eps_s[i] * sigmam().dag() * sigmam())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c_1 = 0
    c_2 = 0
    c_3 = 0
    c = 0
    rho = 0
    H = 0
    for i in range(0, len(delta_s)):
        H = H + tensor(H_ap[i], fock_dm(len(delta_s), i))
        c_1 = c_1 + tensor(np.sqrt(gamma_p[i]) * sigmam(), fock_dm(len(delta_s), i))
        c_2 = c_2 + tensor(np.sqrt(gamma_dec[i]) * sigmam(), fock_dm(len(delta_s), i))
        c_3 = c_3 + tensor(np.sqrt(gamma_phi[i] / 2) * sigmam(), fock_dm(len(delta_s), i))
        c = c + tensor((np.sqrt(2 * kappa_1[i]) * a_ad[i] - beta[i]) * np.exp(-1j * Phi), fock_dm(len(delta_s), i))
        if i == int(len(delta_s) / 2):  # or i == int(3 * len(delta_s) / 2):
            rho = rho + tensor(fock_dm(2, 1), fock_dm(len(delta_s), i))
        else:
            rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta_s), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta_s))
    for i in range(0, len(delta_s), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c)) + \
        1 / 2 * (- spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (- spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
        1 / 2 * (2 * sprepost(c_1.dag(), c_1) - spre(c_1 * c_1.dag()) - spost(c_1 * c_1.dag())) + \
        1 / 2 * (2 * sprepost(c_2.dag(), c_2) - spre(c_2 * c_2.dag()) - spost(c_2 * c_2.dag())) + \
        1 / 2 * (2 * sprepost(c_3.dag(), c_3) - spre(c_3 * c_3.dag()) - spost(c_3 * c_3.dag()))

    ldt = (L * dt).expm()

    # sm = tensor(qeye(2), ba)
    # sp = tensor(qeye(2), ab)
    seed(sd)
    a_exp = []
    b_exp = []
    dW = []  # generates random noise
    dY = []  # homodyne current
    n = 0
    s = off_diag(len(delta_s))
    P_list = np.zeros((len(delta_s), len(time)))
    follow_traj = []
    rho_list = []
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
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

            elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
        elif w == 1:
            if trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                rho = (Jba * rho * Jba.dag())
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()
            elif trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                rho = (Jab * rho * Jab.dag())
                print(t)
                if np.abs(rho.tr()) > 0:
                    rho = rho / rho.tr()

        dW.append(np.random.normal(loc=0.0, scale=1.0, size=None) * np.sqrt(dt))
        dY.append((c * rho + rho * c.dag()).tr() * np.sqrt(eta) * dt + dW[n])

        rho = vector_to_operator(ldt * operator_to_vector(rho))  # + (c * rho + rho * c.dag()) * dY[n] * np.sqrt(eta)
        rho = rho / rho.tr()
        rho_list.append(rho)  # GET THIS BACK FOR ESTIMATION ON THE BLOCH SPHERE
        norm = 0
        P = []
        for i in range(0, len(delta_s)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        test = 0
        follow = 0
        for i in range(0, len(delta_s)):
            P_list[i, n] = P[i]
            if test < P[i]:
                test = P[i]
                follow = i
        follow_traj.append(follow)
        n += 1

    return dY, P_list, follow_traj, rho_list
    # return dY, P_list, follow_traj


def homodyne_update_PQS(time, dt, delta_s, delta_r, g, gamma_dec, gamma_phi, kappa_1, kappa,
                        eta, W_f, W_b, sd, beta, Phi, dY):
    gamma_p = []  # represents the Purcell enhanced damping of the spin
    alpha = []
    a_ad = []
    eps_s = []
    H_ap = []
    for i in range(0, len(delta_s)):
        gamma_p.append(2 * g[i] ** 2 * kappa[i] / (kappa[i] ** 2 + (delta_r[i] - delta_s[i]) ** 2))
        alpha.append(np.sqrt(2 * kappa_1[i]) * beta[i] / (kappa[i] + 1j * delta_r[i]))
        a_ad.append(alpha[i] - 1j * g[i] * sigmaz() / (kappa[i] + 1j * (delta_r[i] - delta_s[i])))
        eps_s.append((delta_r[i] - delta_s[i]) * g[i] ** 2 * kappa[i] / (kappa[i] ** 2
                                                                         + (delta_r[i] - delta_s[i]) ** 2))
        H_ap.append(delta_s[i] * sigmaz() / 2 + g[i] * (
                alpha[i] * sigmam().dag() + alpha[i].conjugate() * sigmam()) - eps_s[i] * sigmam().dag() * sigmam())

    # defining the block diagonal hamiltonian and operators ----------------------------------------
    c_1 = 0
    c_2 = 0
    c_3 = 0
    c = 0
    rho = 0
    H = 0
    for i in range(0, len(delta_s)):
        H = H + tensor(H_ap[i], fock_dm(len(delta_s), i))
        c_1 = c_1 + tensor(np.sqrt(gamma_p[i]) * sigmam(), fock_dm(len(delta_s), i))
        c_2 = c_2 + tensor(np.sqrt(gamma_dec[i]) * sigmam(), fock_dm(len(delta_s), i))
        c_3 = c_3 + tensor(np.sqrt(gamma_phi[i] / 2) * sigmam(), fock_dm(len(delta_s), i))
        c = c + tensor((np.sqrt(2 * kappa_1[i]) * a_ad[i] - beta[i]) * np.exp(-1j * Phi), fock_dm(len(delta_s), i))
        rho = rho + tensor(fock_dm(2, 0), fock_dm(len(delta_s), i))

    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(delta_s))
    for i in range(0, len(delta_s)):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (2 * sprepost(Jab.dag(), Jab) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (2 * sprepost(Jba.dag(), Jba) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
        1 / 2 * (2 * sprepost(c.dag(), c) - spre(c.dag() * c) - spost(c.dag() * c)) + \
        1 / 2 * (2 * sprepost(c_1.dag(), c_1) - spre(c_1 * c_1.dag()) - spost(c_1 * c_1.dag())) + \
        1 / 2 * (2 * sprepost(c_2.dag(), c_2) - spre(c_2 * c_2.dag()) - spost(c_2 * c_2.dag())) + \
        1 / 2 * (2 * sprepost(c_3.dag(), c_3) - spre(c_3 * c_3.dag()) - spost(c_3 * c_3.dag()))

    ldt = (L * dt).expm()

    seed(sd)
    n = 0
    Pa = []
    Pb = []
    rho_save = []
    s = off_diag(len(delta_s))
    P_list = np.zeros((len(delta_s), len(time)))
    follow_traj = []
    for t in time:
        rho = vector_to_operator(ldt * operator_to_vector(rho)) + (c * rho + rho * c.dag()) * dY[n] * np.sqrt(eta)
        rho = rho.unit()
        rho_save.append(rho)

        norm = 0
        P = []
        for i in range(0, len(delta_s)):
            sm = tensor(qeye(2), s[i])
            P.append(np.abs((sm * rho * sm.dag()).tr()))
            norm += np.abs((sm * rho * sm.dag()).tr())

        P = P / norm
        for i in range(0, len(delta_s)):
            P_list[i, n] = P[i]
        n += 1

        test = 0
        follow = 0
        for i in range(0, len(delta_s)):
            if test < P[i]:
                test = P[i]
                follow = i
        follow_traj.append(delta_s[follow])

    # Propagating the past quantum state
    E = 0
    for i in range(0, len(delta_s)):
        E = E + tensor(qeye(2), fock_dm(len(delta_s), i))

    Pb_list = np.zeros((len(delta_s), len(time)))
    Lb = 1j * (spre(H) - spost(H)) + \
         1 / 2 * (2 * sprepost(Jab, Jab.dag()) - spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
         1 / 2 * (2 * sprepost(Jba, Jba.dag()) - spre(Jba * Jba.dag()) - spost(Jba * Jba.dag())) + \
         1 / 2 * (2 * sprepost(c, c.dag()) - spre(c.dag() * c) - spost(c.dag() * c)) + \
         1 / 2 * (2 * sprepost(c_1, c_1.dag()) - spre(c_1 * c_1.dag()) - spost(c_1 * c_1.dag())) + \
         1 / 2 * (2 * sprepost(c_2, c_2.dag()) - spre(c_2 * c_2.dag()) - spost(c_2 * c_2.dag())) + \
         1 / 2 * (2 * sprepost(c_3, c_3.dag()) - spre(c_3 * c_3.dag()) - spost(c_3 * c_3.dag()))

    lbdt = (Lb * dt).expm()
    #    calculates the effect matrix
    nb = len(time) - 1
    follow_trajb = []
    for t in time:
        E = vector_to_operator(lbdt * operator_to_vector(E)) + (c.dag() * E + E * c) * dY[nb] * np.sqrt(eta)
        # E = E.unit()
        P = []
        norm = 0
        for i in range(0, len(delta_s)):
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
            for i in range(0, len(delta_s)):
                E = E + tensor(qeye(2), fock_dm(len(delta_s), i)) * Pb_list[i, nb + 1]
                P[i] = Pb_list[i, nb + 1]
                print(Pb_list[i, nb + 1])
                norm = norm + P[i]
            print(norm)
        elif norm > 0:
            P = P / norm

        for i in range(0, len(delta_s)):
            Pb_list[i, nb] = P[i]

        test = 0
        follow = 0
        for i in range(0, len(delta_s)):
            if test < P[i]:
                test = P[i]
                follow = i
        follow_trajb.append(delta_s[follow])

        nb -= 1

    follow_trajb = rev_vec(follow_trajb)

    return P_list, Pb_list, follow_traj, follow_trajb


def const_tensorspace(A, s):
    c = 0
    for i in range(0, s):
        c = c + tensor(A, fock_dm(s, i))

    return c



def main():

    t_max = 1000
    time = np.linspace(0, t_max, 100000)
    dt = time[1] - time[0]

    '''
    follow_traj_05_2.txt
    follow_traj_08_2.txt
    follow_traj_096_2.txt
    Pb_list_05_2.txt
    Pb_list_08_2.txt
    Pb_list_096_2.txt
    P_list_05_2.txt
    P_list_08_2.txt
    P_list_096_2.txt
    '''
    '''
    time = np.linspace(0, 1000, 10)
    x_pos = np.linspace(0, 4, 10)


    plt.rcParams.update({'font.size': 20})

    random_data1 = np.random.random_sample(size=(len(time), len(x_pos)))*0.28
    random_data2 = np.random.random_sample(size=(len(time), len(x_pos)))*0.35

    fig1 = plt.figure(1)
    ax1 = plt.subplot(211)
    plt.contourf(time, x_pos, random_data1, levels=100, cmap='viridis')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.colorbar()
    plt.ylabel('$\Delta_s$')
    ax2 = plt.subplot(212)
    plt.contourf(time, x_pos, random_data2, levels=100, cmap='viridis')
    plt.xlabel('$\gamma t$')
    plt.ylabel('$\Delta_s$')
    plt.setp(ax2.get_xticklabels(), visible=True)
    plt.colorbar()

    plt.show()
    '''
    '''

    t_max = 1000
    time = np.linspace(0, t_max, 100000)
    size = 21
    delta_s = np.linspace(-2, 2, size)
    
    follow_traj_05 = np.loadtxt(fname='follow_traj_05_2.txt')
    Pb_list_05 = np.loadtxt(fname='Pb_list_05_2.txt')

    follow_traj_08 = np.loadtxt(fname='follow_traj_08_2.txt')
    Pb_list_08 = np.loadtxt(fname='Pb_list_08_2.txt')

    follow_traj_096 = np.loadtxt(fname='follow_traj_096_2.txt')
    Pb_list_096 = np.loadtxt(fname='Pb_list_096_2.txt')
    #follow_traj_096 = np.loadtxt(fname='traj_eta1_6.txt')
    #Pb_list_096 = np.loadtxt(fname='Pb_list_eta1_6.txt')
    

    Pb_list_096[0, 0] = 1
    Pb_list_08[0, 0] = 1
    Pb_list_05[0, 0] = 1

    plt.rcParams.update({'font.size': 20})

    fig3 = plt.figure(3)
    plt.contourf(time, delta_s, Pb_list_05, levels=100, cmap='nipy_spectral')
    plt.colorbar()


    fig1 = plt.figure(1)
    ax1 = plt.subplot(311)
    #plt.contour(time, delta_s, Pb_list_05, levels=100, cmap='nipy_spectral')
    plt.contourf(time, delta_s, Pb_list_05, levels=100, cmap='nipy_spectral')
    plt.plot(time, follow_traj_05, label='p=0.50', color='r')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    #plt.colorbar()
    ax2 = plt.subplot(312)
    #plt.contour(time, delta_s, Pb_list_08, levels=100, cmap='nipy_spectral')
    plt.contourf(time, delta_s, Pb_list_08, levels=100, cmap='nipy_spectral')
    plt.plot(time, follow_traj_08, label='p=0.80', color='r')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    #plt.colorbar()
    ax3 = plt.subplot(313)
    #plt.contour(time, delta_s, Pb_list_096, levels=100, cmap='nipy_spectral')
    plt.contourf(time, delta_s, Pb_list_096, levels=100, cmap='nipy_spectral')
    plt.plot(time, follow_traj_096, label='p=0.93', color='r')
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    plt.xlabel('$\gamma t$')
    #plt.colorbar()


    follow_traj_05 = np.loadtxt(fname='follow_traj_05_2.txt')
    P_list_05 = np.loadtxt(fname='P_list_05_2.txt')

    follow_traj_i08 = np.loadtxt(fname='follow_traj_08_2.txt')
    P_list_08 = np.loadtxt(fname='P_list_08_2.txt')

    #follow_traj_096 = np.loadtxt(fname='traj_eta1.txt')
    P_list_096 = np.loadtxt(fname='P_list_096_2.txt')
    #P_list_096 = np.loadtxt(fname='P_list_eta1_6.txt')

    P_list_05[0, 0] = 1
    P_list_08[0, 0] = 1
    P_list_096[0, 0] = 1

    # CMRmap
    fig2 = plt.figure(2)
    ax1 = plt.subplot(311)
    #plt.contour(time, delta_s, P_list_05, levels=100, cmap='nipy_spectral')
    plt.contourf(time, delta_s, P_list_05, levels=100, cmap='nipy_spectral')
    plt.plot(time, follow_traj_05, label='p=0.50', color='r')
    plt.legend(loc=1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.ylabel('$\Delta_s$')
    #plt.colorbar()
    ax2 = plt.subplot(312)
    plt.contour(time, delta_s, P_list_08, levels=100, cmap='nipy_spectral')
    plt.contourf(time, delta_s, P_list_08, levels=100, cmap='nipy_spectral')
    plt.plot(time, follow_traj_08, label='p=0.80', color='r')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    #plt.colorbar()
    ax3 = plt.subplot(313)
    #plt.contour(time, delta_s, P_list_096, levels=100, cmap='nipy_spectral')
    plt.contourf(time, delta_s, P_list_096, levels=100, cmap='nipy_spectral')
    plt.plot(time, follow_traj_096, label='p=0.93', color='r')
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    plt.xlabel('$\gamma t$')
    #plt.colorbar()
    



    plt.show()

    '''

    '''
    size = 5
    delta_s = np.linspace(0, 4, size)
    delta_r = np.zeros(size)
    g = 2 * np.ones(size)
    kappa = 10 * np.ones(size)
    kappa_1 = kappa
    beta = 1 * np.ones(size)
    gamma_dec = np.linspace(0, 2, size)
    gamma_phi = np.linspace(0, 2, size)
    Phi = np.pi / 2

    eta = 1.
    sd = 17

    W_f1 = []
    W_b1 = []
    seed(sd)
    size2 = int(size * 3 / 2)
    for i in range(0, size2):
        # ehrenfest chain
        s = 0.97
        W_f1.append((i / (size2 - 1)) * (1 - s))
        W_b1.append((1. - i / (size2 - 1)) * (1 - s))

    W_f = []
    W_b = []
    for i in range(0, size):
        #W_f.append(W_f1[int(i + size / 2)])
        #W_b.append(W_b1[int(i + size / 2)])
        W_f.append(0)
        W_b.append(0)
    print(W_f)
    print(W_b)

    '''

    # THIS IS WERE ARE WORKING!!!!!!

    size = 31
    delta_s = np.linspace(-2, 2, size)
    delta_r = np.zeros(size)
    g = 2 * np.ones(size)
    kappa = 10 * np.ones(size)
    kappa_1 = kappa
    beta = 3 * np.ones(size)
    gamma_dec = np.linspace(0, 2, size)
    gamma_phi = np.linspace(0, 2, size)
    Phi = np.pi / 2


    eta = 1.
    sd = 902

    W_f = []
    W_b = []
    seed(sd)
    size2 = int(size * 100/100)
    for i in range(0, size2):
        # ehrenfest chain
        s = 0.90
        W_f.append((i / (size2 - 1)) * (1 - s))
        W_b.append((1. - i / (size2 - 1)) * (1 - s))

    dY, P_real, traj, rho_list = homodyne_emission(time, dt, delta_s, delta_r, g,
                                                   gamma_dec, gamma_phi, kappa_1, kappa, eta, W_f, W_b, sd, beta,
                                                   Phi)
    P_list, Pb_list, f_traj, fb_traj = homodyne_update_PQS(time, dt, delta_s, delta_r,
                                                           g, gamma_dec, gamma_phi, kappa_1,
                                                           kappa, eta, W_f, W_b, sd, beta, Phi, dY)

    np.savetxt('P_list_90_article.txt', P_list)
    np.savetxt('Pb_list_90_article.txt', Pb_list)
    np.savetxt('traj_90_article.txt', traj)
    '''


    P_list = np.loadtxt(fname='P_list_95_article.txt')
    Pb_list = np.loadtxt(fname='Pb_list_95_article.txt')
    traj = np.loadtxt(fname='traj_95_article.txt')

    for i in range(0, len(time)):
        j = int(traj[i])
        traj[i] = delta_s[j]


    fig2 = plt.figure(2)
    plt.contourf(time, delta_s, P_list, levels=100, cmap='nipy_spectral')
    plt.plot(time, traj, label='p=0.50', color='r')
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    plt.colorbar()

    fig3 = plt.figure(3)
    plt.contourf(time, delta_s, Pb_list, levels=100, cmap='nipy_spectral')
    plt.plot(time, traj, label='p=0.50', color='r')
    plt.legend(loc=1)
    plt.ylabel('$\Delta_s$')
    plt.colorbar()

    plt.show()
    '''
    '''
    sx = const_tensorspace(sigmax(), len(delta_s))
    sy = const_tensorspace(sigmay(), len(delta_s))
    sz = const_tensorspace(sigmaz(), len(delta_s))
    exp_sx = []
    exp_sy = []
    exp_sz = []
    b = Bloch()
    for t in range(0, len(time)):
        i = traj[t]
        norm = ((sx * rho_list[t]).tr()) * ((sx * rho_list[t]).tr()) + \
               ((sy * rho_list[t]).tr()) * ((sy * rho_list[t]).tr()) + ((sz * rho_list[t]).tr()) * (
                   (sz * rho_list[t]).tr())
        norm = abs(np.sqrt(norm))
        print(norm)
        exp_sx.append((sx * rho_list[t]).tr() / norm)
        exp_sy.append((sy * rho_list[t]).tr() / norm)
        exp_sz.append((sz * rho_list[t]).tr() / norm)

    dY1, P_real, traj, rho_list = homodyne_emission_steady(time, dt, delta_s, delta_r, g, gamma_dec, gamma_phi,
                                                           kappa_1, kappa, eta, W_f, W_b, sd, beta, Phi)

    sx = const_tensorspace(sigmax(), len(delta_s))
    sy = const_tensorspace(sigmay(), len(delta_s))
    sz = const_tensorspace(sigmaz(), len(delta_s))
    exp_sx1 = []
    exp_sy1 = []
    exp_sz1 = []
    b = Bloch()
    for t in range(0, len(time)):
        i = traj[t]
        norm = ((sx * rho_list[t]).tr()) * ((sx * rho_list[t]).tr()) + \
               ((sy * rho_list[t]).tr()) * ((sy * rho_list[t]).tr()) + ((sz * rho_list[t]).tr()) * (
                   (sz * rho_list[t]).tr())
        norm = abs(np.sqrt(norm))
        print(norm)
        exp_sx1.append((sx * rho_list[t]).tr() / norm)
        exp_sy1.append((sy * rho_list[t]).tr() / norm)
        exp_sz1.append((sz * rho_list[t]).tr() / norm)

    fig1 = plt.figure(1)
    plt.plot(time, exp_sx, label='$\langle \sigma_x \\rangle$')
    plt.plot(time, exp_sy, label='$\langle \sigma_y \\rangle$')
    plt.plot(time, exp_sz, label='$\langle \sigma_z \\rangle$')
    plt.legend(loc=1)

    fig2 = plt.figure(2)
    plt.plot(time, exp_sx1, label='$\langle \sigma_x \\rangle$')
    plt.plot(time, exp_sy1, label='$\langle \sigma_y \\rangle$')
    plt.plot(time, exp_sz1, label='$\langle \sigma_z \\rangle$')
    plt.legend(loc=1)
    '''

    '''
    dec_list = np.linspace(-.5, .5, 10)
    delta = 0
    omega = .2

    size = len(dec_list)
    size2 = int(size * 3/2)
    W_f = []
    W_b = []
    for i in range(0, size2):
        # ehrenfest chain
        s = 0.05
        W_f.append((i / (size2 - 1)) * (1 - s))
        W_b.append((1. - i / (size2 - 1)) * (1 - s))

    delta_s = dec_list
    H = 0
    rho = 0
    rho1 = 0
    for i in range(0, len(dec_list)):
        H = H + tensor(omega * sigmax() + dec_list[i] * sigmaz(), fock_dm(len(dec_list), i))
        if i == int(len(delta_s) / 2):  # or i == int(3 * len(delta_s) / 2):
            rho = rho + tensor(fock_dm(2, 1), fock_dm(len(delta_s), i))
            rho1 = rho1 + tensor(fock_dm(2, 1), fock_dm(len(delta_s), i))
        else:
            rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta_s), i))
            rho1 = rho1 + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta_s), i))
    Jab = 0
    Jba = 0
    J_for, J_bac = J_matrix(len(dec_list))
    for i in range(0, len(dec_list), 1):
        Jab = Jab + tensor(np.sqrt(W_f[i]) * qeye(2), J_for[i])
        Jba = Jba + tensor(np.sqrt(W_b[i]) * qeye(2), J_bac[i])

    # Liovillian with coherent photon emission and coherent jump evolution
    L = -1j * (spre(H) - spost(H)) + \
        1 / 2 * (- spre(Jab * Jab.dag()) - spost(Jab * Jab.dag())) + \
        1 / 2 * (- spre(Jba * Jba.dag()) - spost(Jba * Jba.dag()))

    ldt = (L * dt).expm()
    sx = const_tensorspace(sigmax(), len(delta_s))
    sy = const_tensorspace(sigmay(), len(delta_s))
    sz = const_tensorspace(sigmaz(), len(delta_s))
    #sx = const_tensorspace(sigmap(), len(delta_s))
    e_ops = [sx, sy, sz]

    Lindblad_solve = mesolve(H, rho, time, Jab, e_ops)


    ex_t = []
    ey_t = []
    ez_t = []

    for j in range(0, 5):

        rho = 0
        for i in range(0, len(dec_list)):
            if i == int(len(delta_s) / 2):  # or i == int(3 * len(delta_s) / 2):
                rho = rho + tensor(fock_dm(2, 1), fock_dm(len(delta_s), i))
            else:
                rho = rho + tensor(Qobj([[0, 0], [0, 0]]), fock_dm(len(delta_s), i))

        sd = j
        seed(sd)

        n = 0
        s = off_diag(len(dec_list))
        rho_list = []
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
                    print(t)
                    if np.abs(rho.tr()) > 0:
                        rho = rho / rho.tr()

                elif trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                    rho = (Jba * rho * Jba.dag())
                    print(t)
                    if np.abs(rho.tr()) > 0:
                        rho = rho / rho.tr()
            elif w == 1:
                if trans < np.abs((Jba.dag() * Jba * rho).tr() * dt):
                    rho = (Jba * rho * Jba.dag())
                    print(t)
                    if np.abs(rho.tr()) > 0:
                        rho = rho / rho.tr()
                elif trans < np.abs((Jab.dag() * Jab * rho).tr() * dt):
                    rho = (Jab * rho * Jab.dag())
                    print(t)
                    if np.abs(rho.tr()) > 0:
                        rho = rho / rho.tr()

            rho = vector_to_operator(ldt * operator_to_vector(rho))
            rho = rho / rho.tr()
            rho_list.append(rho)  # GET THIS BACK FOR ESTIMATION ON THE BLOCH SPHERE
            norm = 0

            n += 1

        traj = rho_list


        delta_s = dec_list
        sx = const_tensorspace(sigmax(), len(delta_s))
        sy = const_tensorspace(sigmay(), len(delta_s))
        sz = const_tensorspace(sigmaz(), len(delta_s))
        exp_sx = []
        exp_sy = []
        exp_sz = []
        for t in range(0, len(time)):
            i = traj[t]
            norm = ((sx * rho_list[t]).tr()) * ((sx * rho_list[t]).tr()) + \
                   ((sy * rho_list[t]).tr()) * ((sy * rho_list[t]).tr()) + ((sz * rho_list[t]).tr()) * (
                       (sz * rho_list[t]).tr())
            norm = abs(np.sqrt(norm))
            print(norm)
            exp_sx.append((sx * rho_list[t]).tr() / norm)
            exp_sy.append((sy * rho_list[t]).tr() / norm)
            exp_sz.append((sz * rho_list[t]).tr() / norm)

        ex_t.append(exp_sx)
        ey_t.append(exp_sy)
        ez_t.append(exp_sz)


    sum_x = ex_t[0]
    sum_y = ey_t[0]
    sum_z = ez_t[0]
    for i in range(1, 5):
        x = ex_t[i]
        y = ey_t[i]
        z = ez_t[i]
        for t in range(0, len(time)):
            sum_x[t] = (sum_x[t] + x[t]) / 2
            sum_y[t] = (sum_y[t] + y[t]) / 2
            sum_z[t] = (sum_z[t] + z[t]) / 2


    fig1 = plt.figure(1)
    ax1 = plt.subplot(311)
    for i in range(0, 4):
        plt.plot(time, ex_t[i], label='$\langle \sigma_x \\rangle$', alpha=0.5)
    plt.plot(time, sum_x, label='$\langle \sigma_x \\rangle$')
    plt.plot(time, Lindblad_solve.expect[0], linestyle='dashed', label='Unmonitored')

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.ylabel('$\Delta_s$')
    ax1 = plt.subplot(312)
    for i in range(0, 4):
        plt.plot(time, ey_t[i], label='$\langle \sigma_x \\rangle$', alpha=0.5)
    plt.plot(time, sum_y, label='$\langle \sigma_x \\rangle$')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot(time, Lindblad_solve.expect[1], linestyle='dashed', label='Unmonitored')

    plt.ylabel('$\Delta_s$')
    ax2 = plt.subplot(313)
    for i in range(0, 4):
        plt.plot(time, ez_t[i], label='$\langle \sigma_x \\rangle$', alpha=0.5)
    plt.plot(time, sum_z, label='$\langle \sigma_x \\rangle$')
    plt.plot(time, Lindblad_solve.expect[2], linestyle='dashed', label='Unmonitored')
    plt.xlabel('$\gamma t$')
    plt.ylabel('$\Delta_s$')
    plt.setp(ax2.get_xticklabels(), visible=True)

    fig2 = plt.figure(2)
    for i in range(0, 4):
        plt.plot(time, ex_t[i], label='$\langle \sigma_x \\rangle$', alpha=0.5)
    plt.plot(time, sum_x, label='$\langle \sigma_x \\rangle$')

    plt.show()
    '''

if __name__ == '__main__':
    main()
