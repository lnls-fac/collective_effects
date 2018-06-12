#!/usr/bin/env python-sirius

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.integrate as scy_int
# from pycolleff.colleff import Ring

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
n = 3

C0 = 2*np.pi/wrf*h*c
E0 = 3e9
U0 = 856e3/E0
It = 500e-3

def get_potentials(wr=wrf, Rs=15e6, Q=1e6,
                   Ibs=np.ones(h)*350e-3/h,
                   z=np.linspace(-0.1, 0.1, 1000)):

    T0 = C0/c
    w0 = 2*np.pi/T0
    alpha = wr/2/Q
    wrb = np.sqrt(wr*wr-alpha*alpha)
    beta = (alpha - 1j*wrb)/c

    Ll = 1/(1-np.exp(-beta*C0))  # Lesser
    Gl = Ll*np.exp(-beta*C0)     # Greater

    A = Ll*np.tri(h,h,-1) + Gl*np.tri(h,h).T
    l = np.arange(h)
    dif = l[:,None] - l[None,:]
    V1 = np.exp(-beta*C0*dif/h)
    Vv = np.dot(A*V1,Ibs)

    # sigz = 20e-3**2 * 2 /1.225)**(4/3)
    # dist = np.exp(-(z-z0)**4/sigz)
    # dist /= np.trapz(dist, x=z)
    sigz = 4*10e-3;
    dist = 1/np.sqrt(2*np.pi)/sigz*np.exp(-(z-z0)**2/2/sigz**2)
    sd  = scy_int.cumtrapz(np.exp(beta*z)*dist, x=z, initial=0*1j)
    sd *= np.exp(-beta*z)
    aux = np.ones(h)
    S = aux[:,None]*sd[None,:]
    Vt = np.exp(-beta*z)[None,:] * Vv[:, None] + It/h*S
    return -T0/E0 * 2*alpha*Rs*(Vt.real - alpha/wrb*Vt.imag)

def get_potentials_imp(p0=n*h, wr=wrf, Rs=15e6, Q=1e6,
                   Ibs=np.ones(h)*350e-3/h,
                   z=np.linspace(-0.1, 0.1, 1000)):

    n_terms = 30
    T0 = C0/c
    w0 = 2*np.pi/T0
    l = np.arange(h)
    dif = C0/h*(l[:,None] - l[None,:])
    Vt = np.zeros([h, len(z)])
    for p in range(p0-n_terms, p0+n_terms+1):
        wp = p*w0
        V1 = np.exp(1j*wp*dif/c)
        Vv = np.dot(V1,Ibs)
        Zp  = Rs/(1-1j*Q*(wr/wp - wp/wr))
        V = w0/2/np.pi * Zp*np.exp(1j*wp*z/c)[None,:] * Vv[:, None]
        Vt += 2 * V.real
    return -T0/E0 * Vt

if __name__ == "__main__":
    zlim = 1.50
    z = np.linspace(-zlim, zlim, 10001)

    V0 = 1.63e6/E0
    r = U0/V0
    phi0 = np.pi - np.arcsin(r)
    print('Unperturbed synch phase {0:7.3f}°'.format(phi0/np.pi*180));
    krf = wrf/c
    Vrf = V0*np.sin(phi0 + krf*z) - U0
    k = np.sqrt(1/n**2 - 1/(n**2 - 1)*r**2)
    print('Vc/Vrf = {0:6.4f}'.format(k))
    rc = n**2/(n**2-1)*r
    phi = np.pi - np.arcsin(rc)
    print('synch phase {0:7.3f}°'.format(phi/np.pi*180));
    phih = -(1/n)*np.arctan(n*r/np.sqrt((n**2 - 1)**2 - n**4*r**2))
    psi = np.pi/2 - n*phih
    print('Harm. cav. phase {0:7.3f}°'.format(psi*180/np.pi))
    Rs = k*V0*E0/2/It/np.abs(np.cos(psi))
    # Rs = 2.017e6
    Q = 21600
    print('Shunt Impedance {0:6.4f}M'.format(Rs*1e-6))
    F = k*V0*E0/2/It/np.abs(np.cos(psi))/Rs
    #F = 0.998
    #F = 0.9435
    print('Form factor {0:6.4f}'.format(F))
    df = n*wrf/(2*np.pi)/2/Q*np.tan(psi)
    print('Harm. cav. detuning {0:7.3f}kHz'.format(df*1e-3))

    wr = n*wrf-2*np.pi*df

    Ibs = np.zeros(h, dtype=complex)
    nb = 176
    start = 0
    z0 = (phi-phi0)/wrf*c
    alpha = wr/2/Q
    wrb = np.sqrt(wr*wr-alpha*alpha)
    Ibs[start:start+nb] = F*It/nb*np.exp(-1j*wrb*z0/c)

    #plt.plot(Ibs)
    #plt.show()
    z_n = z - z0
    V_bassi = -2*It*Rs*F*np.cos(psi)*np.cos(n*wrf*z_n/c-psi)/E0

    V = get_potentials(wr=wr, Rs=Rs, Q=Q, Ibs=Ibs, z=z)
    Vf = V.T[:, 0]

    Ibs[start:start+nb] = F*It/nb * np.exp(-1j*n*wrf*z0/c)
    V_i = get_potentials_imp(wr=wr, Rs=Rs, Q=Q, Ibs=Ibs, z=z)
    V_i = V_i.T[:, 0]

    fig = plt.figure(figsize=(10,14))
    gs = gridspec.GridSpec(3, 1)
    gs.update(left=0.10,right=0.95,bottom=0.10, top=0.97,wspace=0.35, hspace=0.25)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0],sharex=ax1, sharey=ax1)
    ax3 = plt.subplot(gs[2,0],sharex=ax1, sharey=ax1)

    Vt_wake = Vrf + Vf
    Vt_imp = Vrf + V_i
    Vt_bassi = Vrf + V_bassi
    ax1.plot(z*krf/np.pi*180, Vf, label='Wake')
    ax1.plot(z*krf/np.pi*180, Vrf, label='Cav')
    ax1.plot(z*krf/np.pi*180, Vt_wake,label='Cav + Wake')
    #ax1.plot(z*krf/np.pi*180, Vt_imp, label='Cav + Imp')
    #ax1.plot(z*krf/np.pi*180, Vt_bassi,label='Cav + Bassi')
    #ax1.plot(z*krf/np.pi*180, (V_i - Vf)*1e4, label='Cav + Wake')
    ax2.plot(z*krf/np.pi*180, V_i, label='Imp')
    ax2.plot(z*krf/np.pi*180, Vrf, label='Cav')
    ax2.plot(z*krf/np.pi*180, Vt_imp, label='Cav + Imp')
    ax3.plot(z*krf/np.pi*180, V_bassi, label='Bassi')
    ax3.plot(z*krf/np.pi*180, Vrf, label='Cav')
    ax3.plot(z*krf/np.pi*180, Vt_bassi,label='Cav + Bassi')
    # ax1.set_xlim([-5.5, -3.5])
    # ax1.set_ylim([-5e-9, 5e-9])
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax3.legend(loc='best')
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    plt.show()
