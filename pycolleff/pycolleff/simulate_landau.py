#!/usr/bin/env python-sirius
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import scipy.integrate as scy_int
import mathphys as _mp
import resource
import sys

c = _mp.constants.light_speed

def memory_limit():
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024 / 2, hard))

def get_memory():
    with open('/proc/meminfo', 'r') as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                free_memory += int(sline[1])
    return free_memory

class Ring:
    def __init__(self):
        self.E0          = 0      #Nominal Energy [eV]
        self.harm_num    = 0      #Harmonic Number
        self.en_lost_rad = 0      #Energy loss per turn [eV]
        self.nom_cur     = 0      #Total current [A]
        self.bunlen      = 0      #Bunch length [m]   #VER COMO property
        self.f_rf        = 0      #Rf frequency [Hz]
        self.mom_cmpct   = 0      #Momentum compaction factor
        self.en_spread   = 0      #Relative Energy Spread
        self.peak_rf     = 0      #Peak RF Voltage [V]
        self.frf         = 0      #Rf frequency [Hz]

    @property
    def synch_phase(self):
        return np.pi - np.arcsin(self.en_lost_rad/self.peak_rf)

    @property
    def C0(self):
        return c*self.harm_num/self.frf

    @property
    def T0(self):
        return self.harm_num/self.frf

    @property
    def w0(self):
        return 2*np.pi/(self.T0)

    @property
    def wrf(self):
        return 2*np.pi*self.frf

    def get_Vrf(self, z):
        V0 = self.peak_rf
        phi0 = self.synch_phase
        krf = self.wrf/c
        return (V0*np.sin(phi0 + krf*z) - self.en_lost_rad) / self.E0


class HarmCav:
    def __init__(self, wrf, n, Q, r=0):
        self._r = r
        self._wrf = wrf
        self.num = n
        self.quality_fact = Q
        self.shunt_imp = 0
        self.psi = 0
        self.form_fact = 0

    @property
    def delta(self):
        return np.tan(self.psi)/2/self.quality_fact

    @delta.setter
    def delta(self, delta):
        self.psi = np.arctan(delta*2*self.quality_fact)

    @property
    def harm_phase(self):
        return (np.pi/2 - self.psi)/self.num

    @harm_phase.setter
    def harm_phase(self, val):
        self.psi = np.pi/2 - self.num*val

    @property
    def wr(self):
        delta = self.delta
        return self.num*self._wrf*(-delta + np.sqrt(delta*delta + 1))

    @property
    def alpha(self):
        return self.wr/2/self.quality_fact

    @property
    def wrb(self):
        wr = self.wr
        alpha = self.alpha
        return np.sqrt(wr*wr-alpha*alpha)

    @property
    def beta(self):
        return (self.alpha - 1j*self.wrb)/c

    @property
    def k(self):
        return np.sqrt(1/self.num**2 - 1/(self.num**2 - 1)*self._r**2)

    @property
    def pert_synch_phase(self):
        rc = self.num**2/(self.num**2-1)*self._r
        return np.pi - np.arcsin(rc)

    @property
    def phase_shift(self):
        return self.pert_synch_phase - (np.pi - np.arcsin(self._r))

    @property
    def length_shift(self):
        return self.phase_shift/self._wrf*c

    def calc_flat_potential(self, ring, F=1):
        r   = ring.en_lost_rad/ring.peak_rf
        wrf = ring.wrf
        It  = ring.nom_cur
        n   = self.num
        Q   = self.quality_fact

        self._r = r
        phih = -(1/n)*np.arctan(n*r/np.sqrt((n**2 - 1)**2 - n**4*r**2))
        self.psi = np.pi/2 - n*phih
        self.shunt_imp = self.k*ring.peak_rf/2/It/np.abs(np.cos(self.psi))/F
        self.form_fact = F

    def print_param(self, ring):
        print('Unperturbed synch phase {0:7.3f}°'.format(ring.synch_phase/np.pi*180));
        print('Vc/Vrf = {0:6.4f}'.format(self.k))
        print('synch phase {0:7.3f}°'.format(self.pert_synch_phase/np.pi*180));
        print('phase shift {0:7.3f}° ({1:7.3f} mm)'.format(self.phase_shift/np.pi*180,
                                                           self.length_shift*1e3))
        print('Harm. cav. phase {0:7.3f}°'.format(self.psi*180/np.pi))
        print('Shunt Impedance {0:6.4f}M'.format(self.shunt_imp*1e-6))
        print('Form factor {0:6.4f}'.format(self.form_fact))
        print('Harm. cav. detuning {0:7.3f}kHz'.format(self.delta*self.wr/2/np.pi*1e-3))


class Lambda:
    def __init__(self, z, I=None, ring=None, lamb=None):
        if lamb is not None and isinstance(lamb, Lambda):
            self.z = np.array(lamb.z)
            self.I = np.array(lamb.I)
            self.dist = np.array(lamb.dist)
        else:
            self.z = np.array(z)
            self.I = np.array(I)
        if ring is not None:
            self.dist = self.create_dist(ring.get_Vrf(z), ring)
        else:
            self.dist = np.zeros([len(I), len(z)], dtype=float)

    def get_gauss_dist(self, z0, sigz):
        lamb = np.exp(-(self.z-z0)**2/2/sigz**2)/np.sqrt(2*np.pi)/sigz
        return np.ones(self.I.shape)[:, None] * lamb[None, :]

    def get_dist(self, V, ring):
        if len(V.shape) < 2:
            V = V[None,:]
        sigE = ring.en_spread
        a    = ring.mom_cmpct
        U = -scy_int.cumtrapz(V, self.z, initial=0)
        U -= np.min(U, axis=1)[:,None]
        rho = np.exp(-U/2/sigE**2/a/ring.C0)
        rho /= np.trapz(rho,self.z)[:, None]
        return rho

    def get_synch_phases(self):
        zl = self.z[None,:]
        return np.trapz(self.dist * zl, self.z)

    def get_bun_lens(self):
        zl = self.z[None,:]
        z_m = self.get_synch_phases()
        z_2 = np.trapz(self.dist * zl * zl, self.z)
        return np.sqrt(z_2 -z_m*z_m)


def get_potentials_wake(ring, harm_cav, lamb):
    E0   = ring.E0
    It   = ring.nom_cur
    wrf  = ring.wrf
    C0   = ring.C0
    T0   = ring.T0
    h    = ring.harm_num
    sigz = ring.bunlen

    Q     = harm_cav.quality_fact
    Rs    = harm_cav.shunt_imp
    alpha = harm_cav.alpha
    wr    = harm_cav.wr
    n     = harm_cav.num
    beta  = harm_cav.beta
    z0    = harm_cav.length_shift
    z     = lamb.z
    dist  = lamb.dist*lamb.I[:,None]

    Ll = 1/(1-np.exp(-beta*C0))  # Lesser
    Gl = Ll*np.exp(-beta*C0)     # Greater

    A = Ll*np.tri(h,h,-1) + Gl*np.tri(h,h).T
    l = np.arange(h)
    dif = l[:,None] - l[None,:]
    B = np.exp(-beta*C0*dif/h)
    lamb_w = np.trapz(dist*np.exp(beta*z)[None,:],z)
    V = np.dot(A*B,lamb_w)

    sd  = scy_int.cumtrapz(np.exp(beta*z)[None,:]*dist, x=z, initial=0*1j)
    Vt = np.exp(-beta*z)[None,:] * (V[:, None] + sd)

    return -T0/E0 * 2*alpha*Rs*(Vt.real - alpha/harm_cav.wrb*Vt.imag)


def get_potentials_imp(ring, harm_cav, lamb, n_terms=20):

    E0   = ring.E0
    It   = ring.nom_cur
    wrf  = ring.wrf
    C0   = ring.C0
    T0   = ring.T0
    w0   = ring.w0
    h    = ring.harm_num
    sigz = ring.bunlen

    Q  = harm_cav.quality_fact
    Rs = harm_cav.shunt_imp
    wr = harm_cav.wr
    n  = harm_cav.num
    beta = harm_cav.beta
    z = lamb.z

    p0 = h*n

    l = np.arange(h)
    dif =l[:,None] - l[None,:]
    Vt = np.zeros([h, len(z)])
    for p in range(p0-n_terms, p0+n_terms+1):
        wp = p*w0
        B = np.exp(1j*wp*dif*C0/h/c)
        lamb_w = lamb.I*np.trapz(lamb.dist*np.exp(-1j*wp*z/c)[None,:],z)
        V = np.dot(B,lamb_w)
        Zp  = Rs/(1-1j*Q*(wr/wp - wp/wr))
        V_hat = (w0/2/np.pi *Zp * np.exp(1j*wp*z/c))[None,:] * V[:, None]
        Vt += 2 * V_hat.real
    return -T0/E0 * Vt


def get_potential_analytic_uni_fill(ring, harm_cav, lamb):
    It = ring.nom_cur
    Rs = harm_cav.shunt_imp
    F  = harm_cav.form_fact
    psi = harm_cav.psi
    n = harm_cav.num
    z0 = harm_cav.length_shift
    wrf = ring.wrf
    z_n = lamb.z-z0
    return -2*It*Rs*F*np.cos(psi)*np.cos(n*wrf*z_n/c-psi)/ring.E0


def main():

    ring = Ring()

    ring.E0 = 3e9
    ring.nom_cur = 500e-3
    ring.en_lost_rad = 856e3
    ring.frf = 99.93e6
    ring.peak_rf = 1.63e6
    ring.harm_num = 176
    ring.bunlen = 10.1e-3
    ring.en_spread = 7.82e-4
    ring.mom_cmpct = 3.07e-4
    wrf = ring.wrf
    krf = wrf/c
    h = ring.harm_num

    hc_flat = HarmCav(wrf, n=3, Q=21600)
    # F = np.exp(-(ring.bunlen/c)**2/2*(hc_flat.num*wrf)**2)
    F = 0.94348
    hc_flat.calc_flat_potential(ring, F=F)
    # hc_flat.print_param(ring)

    # hc_flat.shunt_imp = 4.2e6
    # hc_flat.psi = 96.558*np.pi/180

    zlim = 0.50
    npoints = 10001
    z = np.linspace(-zlim, zlim, npoints)

    Ib = np.zeros(h, dtype=float)
    nfill = h
    start = h-nfill
    Ib[start:start+nfill] = ring.nom_cur/nfill

    lamb = Lambda(z, Ib)

    z0 = hc_flat.length_shift
    sigz = ring.bunlen

    Vrf = ring.get_Vrf(z)
    # ideal_sigz = c/hc_flat.num/wrf*np.sqrt(2*np.log(1/F))
    dist_old = lamb.get_gauss_dist(z0, sigz)
    lamb.dist = np.array(dist_old)
    V_i = get_potentials_imp(ring, hc_flat, lamb)
    Vt_imp = Vrf[None, :] + V_i
    dist_old = lamb.get_dist(Vt_imp, ring)
    epsilon = 1e-6
    param_conv = 7
    n_iters = 100

    plt.figure(figsize=(10, 14))
    gs = gridspec.GridSpec(3, 1)
    gs.update(
        left=0.10, right=0.95, bottom=0.10,
        top=0.97, wspace=0.35, hspace=0.25)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0], sharex=ax1)
    ax3 = plt.subplot(gs[2, 0])
    ph = z*krf/np.pi*180

    dist_new = np.zeros(dist_old.shape)
    cmap = cm.brg(np.linspace(0, 1, n_iters+1))
    for i in range(n_iters):
        lamb.dist = np.array(dist_old)
        # V_w = get_potentials_wake(ring, hc_flat, lamb)
        V_i = get_potentials_imp(ring, hc_flat, lamb)
        # V_analytic = get_potential_analytic_uni_fill(ring, hc_flat, lamb)
        Vt_imp = Vrf[None, :] + V_i

        dist_new = lamb.get_dist(Vt_imp, ring)
        dif = dist_new - dist_old
        res = np.trapz(dif * dif, z)
        conv = np.sqrt(np.mean(res))
        print('{0:03d}: {1:f}'.format(i, conv))
        if conv < epsilon:
            break
        dist_old *= param_conv
        dist_old += dist_new
        dist_old /= param_conv + 1
        lab = ''
        if not i % 10:
            lab = '{0:02d}'.format(i)
            ax1.plot(ph, Vt_imp[0, :], color=cmap[i], label=lab)
            ax2.plot(ph, dist_new[0, :], color=cmap[i], label=lab)

    lamb.dist = np.array(dist_new)
    sigma_z_imp = lamb.get_bun_lens()
    bl_imp = np.mean(sigma_z_imp)
    z_ave = lamb.get_synch_phases()
    z_ave_ave = np.mean(z_ave)

    V_i = get_potentials_imp(ring, hc_flat, lamb)
    Vt_imp = Vrf + V_i[0, :]

    print('sync phase: {0:7.3f}'.format(z_ave_ave*1e3))
    print('bun length: {0:7.3f}'.format(bl_imp*1e3))

    ax1.plot(ph, Vt_imp, label='Final')
    ax2.plot(ph, lamb.dist[0, :], label='Final')
    ax3.plot(sigma_z_imp*1e3)
    ax3.plot(z_ave*1e3)

    # ax1.set_xlim([-0.5, 0.5])
    # ax1.set_ylim([-5e-7, 5e-7])
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    # ax3.legend(loc='best')
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    plt.show()

if __name__ == "__main__":
    memory_limit() # Limitates maximun memory usage to half
    try:
        main()
    except MemoryError:
        sys.stderr.write('\n\nERROR: Memory Exception\n')
        sys.exit(1)
