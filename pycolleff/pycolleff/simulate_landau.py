#!/usr/bin/env python-sirius

import numpy as np
import matplotlib.pyplot as plt

# c = 299792458
# wrf = 2*np.pi*500e6
# h = 864
#
# C0 = 2*np.pi/wrf*h*c
# E0 = 3e9
# It = 350e-3

## MAX IV
c = 299792458
wrf = 2*np.pi*99.93e6
h = 176

C0 = 2*np.pi/wrf*h*c
E0 = 3e9
U0 = 856e3/E0
It = 500e-3


def get_potentials(wr=wrf, Rs=15e6, Q=1e6,
                   Ibs=np.ones(h)*350e-3/h,
                   z=np.linspace(-0.1, 0.1, 1000)):

    T0 = C0/c

    alpha = wr/2/Q
    wrb = np.sqrt(wr*wr-alpha*alpha)
    beta = (alpha + 1j*wrb)/c

    Gl = 1/(1-np.exp(-beta*C0))  # Greater
    Ll = Gl*np.exp(-beta*C0)     # Lesser
    El = Ll + 1/2                # Equal

    V = np.zeros([len(Ibs), len(z)])
    for n in range(len(Ibs)):
        Vhn = 0 + 1j*0
        for i, ib in enumerate(Ibs):
            P = ib*np.exp(beta*(n-i)*C0/h)
            if i > n:
                P *= Gl
            elif i == n:
                P *= El
            elif i < n:
                P *= Ll
            Vhn += P

        Vt = np.exp(beta*z)*Vhn

        V[n] = T0/E0 * 2*alpha*Rs*(Vt.real + alpha/wrb*Vt.imag)
    return V


if __name__ == "__main__":
    wr = 3*wrf + 2*np.pi*28.43e3
    Rs = 2.017e6
    Q = 21600

    Ibs = np.zeros(h)

    nb = 176
    start = 0
    Ibs[start:start+nb] = It/nb

    plt.plot(Ibs)
    plt.show()

    z = np.linspace(-0.06, 0.06, 1001)
    V = get_potentials(wr=wr, Rs=Rs, Q=Q, Ibs=Ibs, z=z)

    offset = V[0, 500]
    plt.plot(V[:, 0])
    plt.show()
    # plt.plot(z, V.T)
    plt.plot(z, -V.T[:, 1:10] + offset)
    plt.plot(z, -V.T[:, 39:49] + offset)
    plt.plot(z, -V.T[:, 78:88] + offset)

    V0 = 1.63e6/E0
    phi0 = np.pi - np.arcsin(U0/V0)
    krf = 2*np.pi*h/C0

    Vrf = V0*np.sin(phi0 + krf*z) - U0
    plt.plot(z, Vrf, linewidth=4, color='k')
    plt.grid(True)
    plt.show()
