#!/usr/bin/env python-sirius
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# import matplotlib.cm as cm
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
        self.E0 = 0               # Nominal Energy [eV]
        self.harm_num = 0         # Harmonic Number
        self.en_lost_rad = 0      # Energy loss per turn [eV]
        self.nom_cur = 0          # Total current [A]
        self.f_rf = 0             # Rf frequency [Hz]
        self.mom_cmpct = 0        # Momentum compaction factor
        self.en_spread = 0        # Relative Energy Spread
        self.peak_rf = 0          # Peak RF Voltage [V]
        self.frf = 0              # Rf frequency [Hz]
        self.synch_tune = 0       # Synchrotron Tune
        self.peak_cur = 0         # Peak bunch current

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
    def ws(self):
        return self.synch_tune*self.w0

    @property
    def wrf(self):
        return 2*np.pi*self.frf

    @property
    def bunlen(self):
        return self.en_spread * self.mom_cmpct*c / self.ws

    def get_Vrf(self, z):
        V0 = self.peak_rf
        phi0 = self.synch_phase
        krf = self.wrf/c
        return (V0*np.sin(phi0 + krf*z) - self.en_lost_rad) / self.E0


class HarmCav:
    def __init__(self, wrf, n, Q, Rs=None, r=0):
        self._r = r
        self._wrf = wrf
        self.num = n
        self.quality_fact = Q
        self.psi = 0
        self.form_fact = 0
        if Rs is not None:
            self.shunt_imp = np.array(Rs)
        else:
            self.shunt_imp = 0

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
    def harm_phase(self, phih):
        self.psi = np.pi/2 - self.num*phih

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
        It = ring.nom_cur
        self.psi = np.pi - np.abs(np.arccos(self.k * ring.peak_rf / 2 / It /
                                            self.shunt_imp/F))
        self.form_fact = F

    def print_param(self, ring):
        print('Unperturbed synch phase {0:7.3f}°'.format(ring.synch_phase /
                                                          np.pi*180))
        print('Vc/Vrf = {0:6.4f}'.format(self.k))
        print('synch phase {0:7.3f}°'.format(self.pert_synch_phase/np.pi*180))
        print('phase shift {0:7.3f}° ({1:7.3f} mm)'.format(self.phase_shift /
                                                            np.pi*180,
                                                            self.length_shift *
                                                            1e3))
        print('Harm. cav. phase {0:7.3f}°'.format(self.harm_phase * 180/np.pi))
        print('Shunt Impedance {0:6.4f}M'.format(self.shunt_imp*1e-6))
        print('Form factor {0:6.4f}'.format(self.form_fact))
        print('Harm. cav. detuning {0:7.3f}kHz'.format(self.delta * self.wr / 2
                                                       / np.pi*1e-3))


class Lambda:
    def __init__(self, z, cur=None, ring=None, lamb=None):
        if lamb is not None and isinstance(lamb, Lambda):
            self.z = np.array(lamb.z)
            self.cur = np.array(lamb.cur)
            self.dist = np.array(lamb.dist)
        else:
            self.z = np.array(z)
            self.cur = np.array(cur)
        if ring is not None:
            self.dist = self.get_dist(ring.get_Vrf(z), ring)
        else:
            self.dist = np.zeros([len(cur), len(z)], dtype=float)

    def get_gauss_dist(self, z0, sigz):
        lamb = np.exp(-(self.z-z0)**2/2/sigz**2)/np.sqrt(2*np.pi)/sigz
        return np.ones(self.cur.shape)[:, None] * lamb[None, :]

    def get_dist(self, V, ring):
        if len(V.shape) < 2:
            V = V[None, :]
        sigE = ring.en_spread
        a = ring.mom_cmpct
        pos = np.array(self.z)
        npo = len(pos)//2
        pos_par = np.delete(pos, npo)

        z_plus = np.split(pos_par, 2)[1]
        z_minus = np.split(pos_par, 2)[0]
        z_plus = np.append(0, z_plus)
        z_minus = np.append(z_minus, 0)
        # U = -scy_int.cumtrapz(V, self.z, initial=0)
        # U -= np.min(U, axis=1)[:, None]
        # rho = np.exp(-U/sigE/sigE/a/ring.C0/2)
        # rho /= np.trapz(rho, self.z)[:, None]
        U_minus = -scy_int.cumtrapz(V[:, 0:npo+1], z_minus, initial=0)
        U_plus = -scy_int.cumtrapz(V[:, npo:2*npo+1], z_plus, initial=0)
        U_plus = np.array(U_plus)
        U_minus = np.array(U_minus)
        U_plus -= np.absolute(U_minus[0, npo])
        U = np.concatenate((U_minus, U_plus), axis=1)
        U_new = np.delete(U, npo, axis=1)
        rho = np.exp(-U_new/sigE/sigE/a/ring.C0)
        rho /= np.trapz(rho, self.z)[:, None]
        return rho

    def get_synch_phases(self):
        zl = self.z[None, :]
        return np.trapz(self.dist * zl, self.z)

    def get_bun_lens(self):
        zl = self.z[None, :]
        z_m = self.get_synch_phases()
        z_2 = np.trapz(self.dist * zl * zl, self.z)
        return np.sqrt(z_2 - z_m * z_m)


def get_potentials_wake(ring, harm_cav, lamb):
    E0 = ring.E0
    It = ring.nom_cur
    wrf = ring.wrf
    C0 = ring.C0
    T0 = ring.T0
    h = ring.harm_num
    sigz = ring.bunlen

    Q = harm_cav.quality_fact
    Rs = harm_cav.shunt_imp
    alpha = harm_cav.alpha
    wr = harm_cav.wr
    n = harm_cav.num
    beta = harm_cav.beta
    z0 = harm_cav.length_shift
    z = lamb.z
    dist = lamb.dist*lamb.cur[:, None]

    Ll = 1/(1-np.exp(-beta*C0))  # Lesser
    Gl = Ll*np.exp(-beta*C0)     # Greater or Equal

    A = Ll*np.tri(h, h, -1) + Gl*np.tri(h, h).T
    ind_b = np.arange(h)
    dif = ind_b[:, None] - ind_b[None, :]
    B = np.exp(-beta*C0*dif/h)
    lamb_w = np.trapz(dist*np.exp(beta*z)[None, :], z)
    V = np.dot(A*B, lamb_w)

    Sn = scy_int.cumtrapz(np.exp(beta*z)[None, :] * dist, x=z, initial=0*1j)
    Vt = np.exp(-beta*z)[None, :] * (V[:, None] + Sn)

    return -T0/E0 * 2*alpha*Rs*(Vt.real - alpha/harm_cav.wrb*Vt.imag)


def get_potentials_imp(ring, harm_cav, lamb, n_terms=20):

    E0 = ring.E0
    C0 = ring.C0
    T0 = ring.T0
    w0 = ring.w0
    h = ring.harm_num

    Q = harm_cav.quality_fact
    Rs = harm_cav.shunt_imp
    wr = harm_cav.wr
    n = harm_cav.num
    z = lamb.z

    p0 = h*n

    ind_b = np.arange(h)
    dif = ind_b[:, None] - ind_b[None, :]
    Vt = np.zeros([h, len(z)])
    for p in range(p0-n_terms, p0+n_terms+1):
        wp = p*w0
        B = np.exp(-1j*wp*dif*C0/h/c)
        lamb_w = lamb.cur*np.trapz(lamb.dist*np.exp(-1j*wp*z/c)[None, :], z)
        V = np.dot(B, lamb_w)
        Zp = Rs/(1 - 1j * Q * (wr / wp - wp / wr))
        V_hat = (w0 / 2 / np.pi * Zp * np.exp(1j*wp*z/c))[None, :] * V[:, None]
        Vt += 2 * V_hat.real
    return -T0/E0 * Vt


def get_potential_analytic_uni_fill(ring, harm_cav, z):
    It = ring.nom_cur
    Rs = harm_cav.shunt_imp
    F = harm_cav.form_fact
    psi = harm_cav.psi
    n = harm_cav.num
    z0 = harm_cav.length_shift
    wrf = ring.wrf
    z_n = z - z0
    return -2*It*Rs*F*np.cos(psi)*np.cos(n*wrf*z_n/c-psi)/ring.E0


def main():

    ring = Ring()

    # MAX IV
    ring.E0 = 3e9
    ring.nom_cur = 500e-3
    ring.en_lost_rad = 856e3
    ring.frf = 99.931e6
    ring.peak_rf = 1.63e6
    ring.harm_num = 176
    # ring.bunlen = 10.1e-3
    ring.en_spread = 7.82e-4
    ring.mom_cmpct = 3.07e-4
    ring.synch_tune = 1.994e-3

    # # MAX II
    # ring.E0 = 1.5e9
    # ring.en_lost_rad = 133.4e3
    # ring.frf = 99.9602e6
    # ring.peak_rf = 358.19e3
    # ring.nom_cur = 140e-3
    # ring.harm_num = 30
    # # ring.bunlen = 10.1e-3
    # ring.en_spread = 7.01e-4
    # ring.mom_cmpct = 3.82e-3
    # ring.synch_tune = 2*np.pi*6699/ring.w0

    # # SLS
    # ring.E0 = 2.4e9
    # # ring.nom_cur = 500e-3
    # ring.en_lost_rad = 530e3
    # ring.frf = 499.632e6
    # ring.peak_rf = 2.1e6
    # ring.harm_num = 480
    # # ring.bunlen = 10.1e-3
    # ring.en_spread = 9e-4
    # ring.mom_cmpct = 6e-4

    r = ring.en_lost_rad/ring.peak_rf

    wrf = ring.wrf
    krf = wrf/c
    h = ring.harm_num
    It = ring.nom_cur

    hc = HarmCav(wrf, n=3, Q=21600, Rs=4.2e6, r=r)  # MAX IV
    # hc = HarmCav(wrf, n=5, Q=21720, Rs=1.57e6, r=r)  # MAX II
    F = np.exp(-(ring.bunlen/c)**2/2*(hc.num*wrf)**2)
    hc.calc_flat_potential(ring, F=F)

    # hc.shunt_imp = 2.017e6
    # hc.psi = 103.717*np.pi/180

    # hc.shunt_imp = 4.2e6
    # hc.psi = 96.558*np.pi/180

    zlim = 0.30
    npoints = 501
    z = np.linspace(-zlim, zlim, npoints)

    Ib = np.zeros(h, dtype=float)
    nfill = h
    start = h-nfill
    Ib[start:start+nfill] = It/h

    z0 = hc.length_shift
    sigz = ring.bunlen

    lamb = Lambda(z, Ib, ring)

    Vrf = ring.get_Vrf(z)
    dist_old = lamb.get_gauss_dist(z0, sigz)
    lamb.dist = np.array(dist_old)

    V_i = get_potentials_imp(ring, hc, lamb)
    Vt_imp = Vrf[None, :] + V_i
    dist_old = lamb.get_dist(Vt_imp, ring)
    F = np.trapz(dist_old*np.exp(1j*hc.num*wrf*z/c), z)
    F /= np.trapz(lamb.dist, z)

    # V_w = get_potentials_wake(ring, hc, lamb)
    # Vt_wake = Vrf[None, :] + V_w
    # dist_old = lamb.get_dist(Vt_wake, ring)
    # F = np.trapz(dist_old*np.exp(1j*hc.num*wrf*z/c), z)
    # F /= np.trapz(lamb.dist, z)

    epsilon = 1e-5
    param_conv = 15
    n_iters = 500

    plt.figure(figsize=(10, 14))
    gs = gridspec.GridSpec(3, 1)
    gs.update(
        left=0.10, right=0.95, bottom=0.10,
        top=0.97, wspace=0.35, hspace=0.25)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0], sharex=ax1)
    ax3 = plt.subplot(gs[2, 0])
    ph = z*krf
    # ph0 = z0*krf

    dist_new = np.zeros(dist_old.shape)
    # cmap = cm.brg(np.linspace(0, 1, n_iters+1))
    for i in range(n_iters):
        lamb.dist = np.array(dist_old)
        lamb.form_fact = np.array(F)
        F_new = np.trapz(lamb.dist*np.exp(1j*hc.num*wrf*z/c), z)
        F_new /= np.trapz(lamb.dist, z)
        F_mod = np.absolute(F_new[0])

        V_i = get_potentials_imp(ring, hc, lamb)
        Vt_imp = Vrf[None, :] + V_i
        dist_new = lamb.get_dist(Vt_imp, ring)

        # V_w = get_potentials_wake(ring, hc, lamb)
        # Vt_wake = Vrf[None, :] + V_w
        # dist_new = lamb.get_dist(Vt_wake, ring)

        dif = dist_new - dist_old
        res = np.trapz(dif * dif, z)
        conv = np.sqrt(np.mean(res))
        print('{0:03d}: {1:f}, F: {2:f}, Psi: {3:f}'.format(i, conv, F_mod,
                                                            hc.psi*180/np.pi))
        if conv < epsilon:
            break

        dist_old *= param_conv
        dist_old += dist_new
        dist_old /= param_conv + 1
        # lab = ''
        # if not i % 10:
        #     lab = '{0:02d}'.format(i)
        #     ax1.plot(ph, Vt_imp[0, :], color=cmap[i], label=lab)
        #     ax2.plot(ph, dist_new[0, :], color=cmap[i], label=lab)

    lamb.dist = np.array(dist_new)
    sigma_z_imp = lamb.get_bun_lens()
    bl_imp = np.mean(sigma_z_imp)
    z_ave_i = lamb.get_synch_phases()
    z_ave_ave_i = np.mean(z_ave_i)

    print('WAKE')
    hc.print_param(ring)
    V_i = get_potentials_imp(ring, hc, lamb)
    Vt_imp = Vrf + V_i[0, :]

    # V_w = get_potentials_wake(ring, hc, lamb)
    # Vt_wake = Vrf + V_w[0, :]

    print('sync phase imp: {0:7.3f} mm'.format(z_ave_ave_i*1e3))
    print('bun length imp: {0:7.3f} mm, ({1:7.3f} ps)'.format(bl_imp*1e3,
                                                              bl_imp*1e12/c))

    print('=============================')

    hc_flat = HarmCav(wrf, n=3, Q=21600, Rs=4.2e6, r=r)

    hc_flat.calc_flat_potential(ring, F=F_mod)

    print('ANALYTICAL')
    hc_flat.print_param(ring)
    V_a = get_potential_analytic_uni_fill(ring, hc_flat, z)
    Vt_a = Vrf + V_a
    lamb.dist = lamb.get_dist(Vt_a, ring)
    sigma_z_a = lamb.get_bun_lens()
    bl_a = np.mean(sigma_z_a)
    z_ave_a = lamb.get_synch_phases()
    z_ave_ave_a = np.mean(z_ave_a)

    print('sync phase analy: {0:7.3f}'.format(z_ave_ave_a*1e3))
    print('bun length analy: {0:7.3f} mm, ({1:7.3f} ps)'.format(bl_a*1e3,
                                                                bl_a*1e12/c))

    ax1.plot(ph, Vt_imp, label='Final Imp')
    # ax1.plot(ph, Vt_wake, label='Final Wake')
    ax1.plot(ph, Vt_a, label='Final Analytic')

    ax2.plot(ph, dist_new[0, :], '--', label='Final Wake')
    # ax2.plot(ph, dist_new[50, :], label='Final Imp - 50')
    # ax2.plot(ph, dist_new[100, :], label='Final Imp - 100')
    # ax2.plot(ph, dist_new[150, :], label='Final Imp - 150')
    ax2.plot(ph, lamb.dist[0, :], label='Final Analytic')

    ax3.plot(np.sqrt(sigma_z_imp)*1e3, label='Imp')
    ax3.plot(z_ave_i*1e3, label='Imp')
    ax3.plot(np.sqrt(sigma_z_a)*1e3, label='Analytic')
    ax3.plot(z_ave_a*1e3, label='Analytic')

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
    memory_limit()   # Limitates maximun memory usage to half
    try:
        main()
    except MemoryError:
        sys.stderr.write('\n\nERROR: Memory Exception\n')
sys.exit(1)
