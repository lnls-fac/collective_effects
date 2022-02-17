#!/usr/bin/env python-sirius

import matplotlib.pyplot as plt
import matplotlib.gridspec as mpl_gs

import numpy as _np
import mathphys

import pycolleff.impedances as imp
import pycolleff.sirius as si
import scipy.integrate as scy_int

_c = mathphys.constants.light_speed
_pi = _np.pi
_electron_rest_energy = mathphys.constants.electron_rest_energy
_electron_charge = mathphys.constants.elementary_charge


class Params:
    """."""

    def __init__(self):
        """."""
        self._U0 = None
        self._Vrf = None
        self._frf = None
        self._h = None
        self.I0 = None
        self.E0 = None
        self.alpha = None
        self.sync_phase = None
        self.tunes = None
        self.bunlen = None
        self.espread = None
        self.Q = None
        self.Rs = None
        self.nharm = None

    def sirius_params(self):
        """."""
        self.U0 = 871e3
        self.Vrf = 3e6
        self.frf = 499.663824e6
        self.h = 864
        self.I0 = 350e-3
        self.E0 = 3e9
        self.alpha = 1.6446539e-4
        self.sync_phase = 163.1219 * _pi/180
        self.tunes = 4.6520358017681e-3
        self.bunlen = 2.46e-3
        self.espread = 8.43589e-4
        self.nharm = 3
        self.dampte = 13e-3

    def maxiv_params(self):
        """."""
        self.U0 = 856e3
        self.Vrf = 1.63e6
        self.frf = 99.93e6
        self.h = 176
        self.I0 = 500e-3
        self.E0 = 3e9
        self.alpha = 3.07e-4
        self.sync_phase = 148.32 * _pi/180
        self.tunes = 1.994e-3
        self.bunlen = 10.1e-3
        self.espread = 7.82e-4
        self.Q = 21600
        self.nharm = 3
        self.dampte = 25.6e-3

    @property
    def U0(self):
        """."""
        return self._U0

    @U0.setter
    def U0(self, value):
        self._U0 = value

    @property
    def Vrf(self):
        """."""
        return self._Vrf

    @Vrf.setter
    def Vrf(self, value):
        self._Vrf = value

    @property
    def frf(self):
        """."""
        return self._frf

    @frf.setter
    def frf(self, value):
        self._frf = value

    @property
    def h(self):
        """."""
        return self._h

    @h.setter
    def h(self, value):
        self._h = value

    @property
    def wrf(self):
        """."""
        return 2*_pi*self.frf

    @property
    def Trf(self):
        """."""
        return 1/self.frf

    @property
    def w0(self):
        """."""
        return self.wrf/self.h

    @property
    def f0(self):
        """."""
        return self.w0/2/_pi

    @property
    def T0(self):
        """."""
        return 1/self.f0

    @property
    def L0(self):
        """."""
        return self.beta * _c * self.T0

    @property
    def gamma(self):
        """."""
        return self.E0/(_electron_rest_energy/_electron_charge)

    @property
    def beta(self):
        """."""
        return (1-1/self.gamma**2)**(1/2)

    @property
    def buntime(self):
        """."""
        return self.beta * _c * self.bunlen

    @property
    def ws(self):
        """."""
        return self.tunes * self.w0

    @property
    def Ib(self):
        """."""
        return self.I0/self.h


class HarmonicCavity:
    """."""

    def __init__(self, params=None):
        """."""
        self.params = params or Params()

        self._form_factor = 1 + 1j*0
        self._wr = 0
        self._psih_harmonic = 0
        self._harmonic_phase = 0

    @property
    def k_harmonic_flat_potential(self):
        """."""
        nharm = self.params.nharm
        U0 = self.params.U0
        Vrf = self.params.Vrf
        kharm = 1/nharm**2 - ((U0/Vrf)**2)/(nharm**2-1)
        return kharm**(1/2)

    @property
    def phih_harmonic_flat_potential(self):
        """."""
        nharm = self.params.nharm
        U0 = self.params.U0
        Vrf = self.params.Vrf
        q = Vrf/U0
        numer = nharm/q
        denom = (nharm**2-1)**2 - (nharm**2/q)**2
        tan = -numer/(denom)**(1/2)
        return _np.arctan(tan)/nharm

    @property
    def harmonic_phase(self):
        """."""
        return self._harmonic_phase

    @harmonic_phase.setter
    def harmonic_phase(self, value):
        self._harmonic_phase = value

    @property
    def psih_harmonic_flat_potential(self):
        """."""
        return self.params.nharm * self.phih_harmonic_flat_potential - _pi/2

    @property
    def psih_harmonic(self):
        """."""
        return self._psih_harmonic

    @psih_harmonic.setter
    def psih_harmonic(self, value):
        self._psih_harmonic = value

    @property
    def perturbed_sync_phase(self):
        """."""
        nharm = self.params.nharm
        U0 = self.params.U0
        Vrf = self.params.Vrf
        phih = self.harmonic_phase
        kh = self.k_harmonic_flat_potential
        arg = U0/Vrf - kh*_np.sin(nharm*phih)
        return _pi - _np.arcsin(arg)

    @property
    def shunt_impedance(self):
        """."""
        kharm = self.k_harmonic_flat_potential
        psih = self.detune_angle
        Vrf = self.params.Vrf
        I0 = self.params.I0
        ib = 2*I0*abs(self.form_factor)
        return kharm * Vrf / (ib*abs(_np.cos(psih)))

    @property
    def wr(self):
        """."""
        return self._wr

    @wr.setter
    def wr(self, value):
        self._wr = value

    @property
    def form_factor(self):
        """."""
        return self._form_factor

    @form_factor.setter
    def form_factor(self, value):
        self._form_factor = value

    @property
    def detune_angle(self):
        """."""
        Q = self.params.Q
        nharm = self.params.nharm
        wrf = 2*_pi*self.params.frf
        wr = self.wr
        if wr == 0:
            raise Exception('wr cannot be zero!')
        if wrf == 0:
            raise Exception('wrf cannot be zero!')
        tan = Q * (wr/(nharm*wrf) - nharm*wrf/wr)
        angle = _np.arctan2(tan, 1)
        return _pi + angle

    @detune_angle.setter
    def detune_angle(self, value):
        Q = self.params.Q
        nharm = self.params.nharm
        wrf = 2*_pi*self.params.frf
        alpha = _np.tan(value)/Q
        self.wr = nharm*wrf/2
        self.wr *= (alpha + (alpha**2 + 4)**(1/2))

    def detune_passive_cavity(self, Rs, k=None):
        """."""
        kh = k or self.k_harmonic_flat_potential
        I0 = self.params.I0
        Vrf = self.params.Vrf
        ib = 2*I0*abs(self.form_factor)
        arg = kh*Vrf/(ib*Rs)
        return _np.arccos(arg)

    def print_flat_potential(self):
        """."""
        st = ''
        st += 'form factor                   : {:+08.3f} \n'
        st += 'k flat potential              : {:+08.3f} \n'
        st += 'harmonic phase [deg]          : {:+08.3f} \n'
        st += 'harmonic tuning angle [deg]   : {:+08.3f} \n'
        st += 'detuning frequency    [kHz]   : {:+08.3f} \n'
        st += 'shunt impedance flat [M.ohm]  : {:+08.4f} \n'
        st += 'unperturbed sync. phase [deg] : {:+08.3f} \n'
        st += 'perturbed sync. phase [deg]   : {:+08.3f} \n'
        rad2deg = 180/_pi

        kh = self.k_harmonic_flat_potential
        phih = self.phih_harmonic_flat_potential
        psih = self.psih_harmonic_flat_potential
        self.detune_angle = psih
        fr = self.wr/2/_pi
        df = (fr-self.params.nharm*self.params.frf)*1e-3
        Rs_fp = self.shunt_impedance*1e-6
        phis = self.params.sync_phase
        new_phis = self.perturbed_sync_phase
        ffactor = self.form_factor
        print(
            st.format(
                ffactor, kh, phih*rad2deg, psih*rad2deg, df, Rs_fp,
                phis*rad2deg, new_phis*rad2deg))

    def integrated_potential(self, z, voltage=None, harmonic=True):
        """."""
        alpha = self.params.alpha
        if voltage is None:
            wrf = 2*_pi*self.params.frf
            sigmae = self.params.espread
            phis = self.params.sync_phase
            h = self.params.h
            tunes = self.params.tunes
            new_phis = self.perturbed_sync_phase
            nh = self.params.nharm
            phih = self.harmonic_phase
            kh = 0
            if harmonic:
                kh = self.k_harmonic_flat_potential
            phase = wrf*z/_c
            pot = (alpha**2 * sigmae**2)
            pot /= _np.cos(phis)*(h*alpha*sigmae/tunes)**2
            t1 = _np.cos(new_phis)
            t2 = _np.cos(phase + new_phis)
            t3 = kh/nh * (_np.cos(nh*phih) - _np.cos(nh*phase + nh*phih))
            t4 = (_np.sin(new_phis) + kh*_np.sin(nh*phih))*phase
            pot *= (t1 - t2 + t3 - t4)
        else:
            dpot = (voltage - self.params.U0)/self.params.E0
            pot = -scy_int.cumtrapz(dpot, z, initial=0)
            pot -= _np.min(pot)
            pot *= alpha/self.params.L0
        return pot

    def calc_distribution(self, z, voltage=None, harmonic=True):
        """."""
        alpha = self.params.alpha
        espread = self.params.espread
        pot = self.integrated_potential(
            z=z, voltage=voltage, harmonic=harmonic)
        dist = _np.exp(-pot/(alpha**2 * espread**2))
        dist /= _np.trapz(dist, z)
        return dist.ravel()

    @staticmethod
    def calc_sync_phase(z, dist):
        """."""
        return _np.trapz(z*dist, z)

    @staticmethod
    def calc_bunch_length(z, dist):
        """."""
        zm = HarmonicCavity.calc_sync_phase(z, dist)
        z2 = _np.trapz(z**2 * dist, z)
        return _np.sqrt(z2 - zm**2)

    @staticmethod
    def complex_form_factor(z, w, rho):
        """."""
        return _np.trapz(
            rho*_np.exp(1j*w*z/_c), z)/_np.trapz(rho, z)

    def loop_form_factor(
            self, z, dist0, kh=None, include_phase=True, update_detune=False):
        """."""
        wrf = self.params.wrf
        nharm = self.params.nharm
        Rs = self.params.Rs
        Vrf = self.params.Vrf
        self.form_factor = self.complex_form_factor(
            z, nharm*wrf, dist0)

        if update_detune:
            k = kh or self.k_harmonic_flat_potential
            detune = self.detune_passive_cavity(self.params.Rs, k)
            self.detune_angle = detune

        ffang = _np.angle(self.form_factor)
        if not include_phase:
            ffang = 0

        self.harmonic_phase = (self.detune_angle + _pi/2)/nharm

        vharm0 = self.calc_passive_voltage(
            z=0, Rs=Rs, detune_phase=self.detune_angle,
            form_factor_phase=-ffang)
        arg = (self.params.U0 - vharm0)/Vrf
        phis_pert = _pi - _np.arcsin(arg)
        vmain_pert = Vrf * _np.sin(wrf*z/_c + phis_pert)

        vharm = self.calc_passive_voltage(
            z=z, Rs=Rs, detune_phase=self.detune_angle,
            form_factor_phase=-ffang)
        dist = self.calc_distribution(z=z, voltage=vmain_pert+vharm)
        return vmain_pert, vharm, dist

    def convergence_distribution(
            self, z, rho0, niter, tol, beta, method='type1',
            kh=None, include_phase=True, update_detune=False,
            print_iter=False):
        """."""
        err_hist = []
        entropy_hist = []
        rho_n = rho0
        X_k = []
        G_k = []

        def haissinski(
                rho0, z=z, kh=kh,
                include_phase=include_phase,
                update_detune=update_detune):
            _, _, rho_new = self.loop_form_factor(
                z, rho0, kh=kh,
                include_phase=include_phase, update_detune=update_detune)
            return rho_new

        xs_ = [rho0, haissinski(rho0)]
        gs_ = [haissinski(rhoval) - rhoval for rhoval in xs_]
        G_k = [gs_[1] - gs_[0], ]
        X_k = [xs_[1] - xs_[0], ]
        m = 3
        for nit in range(niter):
            rho = haissinski(rho0)
            rho_n_1 = rho
            diff = rho_n_1 - rho_n
            err = _np.sqrt(_np.trapz(diff**2, z))
            err *= (1+beta)
            sint = rho*_np.log(rho)
            filt = _np.logical_not(_np.isnan(sint))
            entropy = -_np.trapz(sint[filt], z[filt])
            err_hist.append(err)
            entropy_hist.append(entropy)
            if print_iter:
                print('{:03d}: {:.2e}'.format(nit, err))
            if err < tol:
                break
            if method == 'type1':
                rho_n_2 = rho_n
                rho_n = (rho_n_1 + beta*rho_n_2)/(1+beta)
            elif method == 'type2':
                _, _, rho_n = haissinski(rho_n)
            elif method == 'anderson_acc':
                mk_ = min(m, nit)
                Xkmat = _np.array(X_k, ndmin=3)
                Gkmat = _np.array(G_k, ndmin=3)
                gamma_k = []
                xk_ = []
                for grid, _ in enumerate(z):
                    # print(Gkmat[:, :, grid])
                    # print(gs_[nit][grid])
                    gamma_ki = _np.linalg.lstsq(
                        Gkmat[:, :, grid], [gs_[nit][grid]], rcond=None)[0]
                    gamma_k.append(gamma_ki)
                    xki_p1 = xs_[nit][grid] + gs_[nit][grid]
                    mat = Xkmat[:, :, grid] + Gkmat[:, :, grid]
                    xki_p1 -= mat @ gamma_ki
                    xki_p1 *= beta
                    xki_p2 = xs_[nit][grid] - Xkmat[:, :, grid] @ gamma_ki
                    xki_p2 *= (1-beta)
                    xk_.append((xki_p1 + xki_p2)[0])
                xs_.append(xk_)
                gs_.append(haissinski(xs_[-1]) - xs_[-1])
                G_k.append(gs_[nit+1] - gs_[nit])
                X_k.append(_np.array(xs_[nit+1]) - _np.array(xs_[nit]))
                rho = xs_[-1]
                rho_n = rho
                if len(X_k) > mk_:
                    X_k = X_k[-m:]
                    G_k = G_k[-m:]
                # print(err)
                if abs(err) < tol:
                    break
            rho0 = rho_n
        return rho, _np.array(err_hist), _np.array(entropy_hist), nit

    def robinson_growth_rate(self, w, wr, approx=False):
        """."""
        alpha = self.params.alpha
        h = self.params.h
        I0 = self.params.I0
        E0 = self.params.E0
        Rs = self.params.Rs
        Q = self.params.Q
        w0 = 2*_pi*self.params.frf/h
        ws = self.params.tunes * w0
        wp = w + ws
        wn = w - ws
        const = I0*alpha*w0/(4*_pi*ws*E0)
        if approx:
            x = w/wr
            const_approx = const*4*ws
            growth = const_approx
            growth *= Rs*Q**2
            growth *= (1-x**2)*(1 + x**2)
            growth /= x**4 * (1 + Q**2 * (1/x - x)**2)**2
        else:
            Zlp = imp.longitudinal_resonator(Rs=Rs, Q=Q, wr=wr, w=wp)
            Zln = imp.longitudinal_resonator(Rs=Rs, Q=Q, wr=wr, w=wn)
            growth = const*(wp*Zlp.real - wn*Zln.real)
        return growth

    def tuneshifts_cbi(self, w, wr, m=1, nbun_fill=None, radiation=False):
        """."""
        ring = si.create_ring()
        ring.nus = self.params.tunes
        ring.w0 = 2*_pi*self.params.frf/self.params.h
        ring.mom_cmpct = self.params.alpha
        ring.E = self.params.E0
        ring.nbun = nbun_fill if nbun_fill is not None else self.params.h
        ring.nom_cur = self.params.I0
        ring.bunlen = self.params.bunlen
        ring.dampte = self.params.dampte
        if not radiation:
            ring.dampte = _np.inf
        Rs = self.params.Rs
        Q = self.params.Q
        Zl = imp.longitudinal_resonator(Rs=Rs, Q=Q, wr=wr, w=w)
        deltaw, wp, interpol_Z, spectrum = ring.longitudinal_cbi(
            w=w, Zl=Zl, bunlen=ring.bunlen(), m=m, inverse=False)

        # Relative tune-shifts must be multiplied by ws
        deltaw *= self.params.tunes * ring.w0
        return deltaw, Zl, wp, interpol_Z, spectrum

    def calc_voltages(self, z):
        """."""
        Vrf = self.params.Vrf
        wrf = 2*_pi*self.params.frf
        phis0 = self.params.sync_phase
        phis_pert = self.perturbed_sync_phase
        kh = self.k_harmonic_flat_potential
        nh = self.params.nharm
        phih = self.harmonic_phase

        phase = wrf*z/_c
        Vmain0 = Vrf*_np.sin(phase + phis0)
        Vmain_pert = Vrf*_np.sin(phase + phis_pert)
        Vharm = Vrf*kh*_np.sin(nh*phase + nh*phih)
        return Vmain0, Vmain_pert, Vharm

    def calc_passive_voltage(self, z, Rs, detune_phase, form_factor_phase=0):
        """."""
        I0 = self.params.I0
        wrf = 2*_pi*self.params.frf
        phase = wrf*z/_c
        nh = self.params.nharm
        ib = 2*I0*abs(self.form_factor)
        volt = -ib*Rs*_np.cos(detune_phase)
        volt *= _np.cos(nh*phase+detune_phase+form_factor_phase)
        return volt

    def plot_voltages(self, z, vmain, vharm, vtotal):
        """."""
        fig = plt.figure(figsize=(6, 4))
        gs = mpl_gs.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0, 0])

        ax1.plot(z*100, vmain * 1e-6, label=r'$V_{rf}$ Main', color='C0')
        ax1.plot(z*100, vharm * 1e-6, label=r'$V_{3h}$ 3HC', color='C1')
        ax1.plot(
            z*100, vtotal * 1e-6, label=r'$V_{rf} + V_{3h}$ Total', color='C3')
        ax1.axhline(
            y=self.params.U0*1e-6, color='tab:gray', ls='--', label=r'$U_0$')

        ax1.set_xlabel(r'$z$ [cm]')
        ax1.set_ylabel('RF voltage [MV]')
        ax1.legend(loc='upper right')
        ax1.grid(ls='--', alpha=0.5)
        return fig
