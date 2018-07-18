import numpy as np
import scipy.integrate as scy_int
import mathphys as _mp
import resource
import simulate_landau as landau_cav

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
        self.frf = 0              # Rf frequency [Hz]
        self.mom_cmpct = 0        # Momentum compaction factor
        self.en_spread = 0        # Relative Energy Spread
        self.peak_rf = 0          # Peak RF Voltage [V]
        self.peak_cur = 0         # Peak bunch current [A]

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
    def synch_tune(self):
        return np.sqrt(self.mom_cmpct * self.peak_rf * self.wrf *
                       np.abs(np.cos(self.synch_phase))
                       / self.T0 / self.E0) / self.w0

    @property
    def ws(self):
        return self.synch_tune*self.w0

    @property
    def wrf(self):
        return 2*np.pi*self.frf

    @property
    def bunlen(self):
        return self.en_spread * self.mom_cmpct * c / self.ws

    def get_Vrf(self, z):
        V0 = self.peak_rf
        phi0 = self.synch_phase
        krf = self.wrf/c
        return (V0*np.sin(phi0 + krf*z) - self.en_lost_rad) / self.E0


class HarmCav:
    def __init__(self, wrf, n, Q, Rs=0, r=0):
        self._r = r
        self._wrf = wrf
        self.num = n
        self.quality_fact = Q
        self.psi = 0
        self.form_fact = 0
        self.shunt_imp = Rs

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
    def dw(self):
        return self.num * self._wrf - self.wr

    @dw.setter
    def dw(self, dw):
        self.wr = self.num * self._wrf - dw

    @property
    def wr(self):
        delta = self.delta
        return self.num*self._wrf*(-delta + np.sqrt(delta*delta + 1))

    @wr.setter
    def wr(self, wr):
        x = wr / self.num / self._wrf
        self.delta = 1 / 2 * (x - 1 / x)

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
        self.harm_phase = (1/self.num)*np.arctan(-self.num*self._r /
                                                 (np.sqrt((self.num**2-1)**2 -
                                                  self.num**4*self._r**2)))
        self.psi = np.pi/2 - self.num * self.harm_phase
        # self.psi = np.pi - np.abs(np.arccos(self.k * ring.peak_rf / 2 / It /
                                            # self.shunt_imp/F))
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
        print('Harm. cav. detuning {0:7.3f}kHz'.format(self.dw / 2
                                                       / np.pi*1e-3))
        print('Harm. cav. psi {0:7.3f}º'.format(self.psi * 180 / np.pi))


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

        U = -scy_int.cumtrapz(V, self.z, initial=0)
        U -= np.min(U, axis=1)[:, None]

        rho = np.exp(-U/sigE/sigE/a/ring.C0)
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
    C0 = ring.C0
    T0 = ring.T0
    h = ring.harm_num
    Rs = harm_cav.shunt_imp
    alpha = harm_cav.alpha
    beta = harm_cav.beta
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


def get_potentials_imp(ring, harm_cav, lamb, n_terms=5):

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
    cur = np.array(lamb.cur)

    p0 = h*n

    ind_b = np.arange(h)
    dif = ind_b[:, None] - ind_b[None, :]
    Vt = np.zeros([h, len(z)])

    if np.all(cur[1:] == cur[:-1]):
        n_terms = 0

    for p in range(p0-n_terms, p0+n_terms+1):
        wp = p*w0
        B = np.exp(-1j * wp * dif * C0 / h / c)
        lamb_w = cur * np.trapz(lamb.dist*np.exp(-1j * wp * z / c)[None, :], z)
        V = np.dot(B, lamb_w)
        Zp = Rs/(1 - 1j * Q * (wr / wp - wp / wr))
        V_hat = (w0 / 2 / np.pi * Zp * np.exp(1j*wp*z/c))[None, :] * V[:, None]
        Vt += 2 * V_hat.real
    return - T0 / E0 * Vt


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


def calc_equilibrium_potential(ring, lamb, hc, z, epsilon=1e-6, param_conv=15,
                               n_iters=1000):

    Vrf = ring.get_Vrf(z)
    dist_old = lamb.get_gauss_dist(hc.length_shift, ring.bunlen)
    lamb.dist = np.array(dist_old)

    V_i = landau_cav.get_potentials_imp(ring, hc, lamb)
    Vt_imp = Vrf[None, :] + V_i
    dist_old = lamb.get_dist(Vt_imp, ring)

    dist_new = np.zeros(dist_old.shape)

    for i in range(n_iters):
        lamb.dist = np.array(dist_old)
        F_new = np.trapz(lamb.dist*np.exp(1j*hc.num*ring.wrf*z/c), z)
        F_new /= np.trapz(lamb.dist, z)
        F_mod = np.absolute(np.mean(F_new))

        V_i = landau_cav.get_potentials_imp(ring, hc, lamb)
        Vt_imp = Vrf[None, :] + V_i
        dist_new = lamb.get_dist(Vt_imp, ring)

        # V_w = landau_cav.get_potentials_wake(ring, hc, lamb)
        # Vt_wake = Vrf[None, :] + V_w
        # dist_new = lamb.get_dist(Vt_wake, ring)

        dif = dist_new - dist_old
        res = np.trapz(np.absolute(dif), z)
        conv = np.mean(res)
        print('{0:03d}: {1:f}, F: {2:f}'.format(i, conv, F_mod))

        if conv < epsilon:
            break

        dist_old *= param_conv
        dist_old += dist_new
        dist_old /= param_conv + 1

    return V_i, Vt_imp, dist_new
