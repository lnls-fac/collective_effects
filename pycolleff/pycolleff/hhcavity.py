"""."""

from functools import partial as _partial

import numpy as _np
import scipy.integrate as scy_int
import matplotlib.pyplot as plt
import matplotlib.gridspec as mpl_gs

import mathphys

import pycolleff.impedances as imp
import pycolleff.sirius as si

_LSPEED = mathphys.constants.light_speed
_PI = _np.pi


class HarmonicCavity:
    """."""

    def __init__(self, ring=None):
        """."""
        self.ring = ring or si.create_ring()

        self.cavity_harm_num = 3
        self.cavity_ang_freq = 0.0
        self.cavity_Q = 0.0
        self.cavity_shunt_impedance = 0.0
        self._form_factor = 1.0 + 0.0j
        self._psih_harmonic = 0.0
        self._harmonic_phase = 0.0

    @property
    def flat_potential_k_harmonic(self):
        """."""
        nharm = self.cavity_harm_num
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        kharm = 1/nharm**2 - ((U0/Vrf)**2)/(nharm**2-1)
        return kharm**(1/2)

    @property
    def flat_potential_phih_harmonic(self):
        """."""
        nharm = self.cavity_harm_num
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        q = Vrf/U0
        numer = nharm/q
        denom = (nharm**2-1)**2 - (nharm**2/q)**2
        tan = -numer/(denom)**(1/2)
        return _np.arctan(tan)/nharm

    @property
    def flat_potential_psih_harmonic(self):
        """."""
        return self.cavity_harm_num * self.flat_potential_phih_harmonic - _PI/2

    @property
    def flat_potential_shunt_impedance(self):
        """."""
        kharm = self.flat_potential_k_harmonic
        psih = self.detune_angle
        Vrf = self.ring.gap_voltage
        I0 = self.ring.total_current
        ib = 2*I0*abs(self.form_factor)
        return kharm * Vrf / (ib*abs(_np.cos(psih)))

    @property
    def harmonic_phase(self):
        """."""
        return self._harmonic_phase

    @harmonic_phase.setter
    def harmonic_phase(self, value):
        self._harmonic_phase = value

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
        nharm = self.cavity_harm_num
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        phih = self.harmonic_phase
        kh = self.flat_potential_k_harmonic
        arg = U0/Vrf - kh*_np.sin(nharm*phih)
        return _PI - _np.arcsin(arg)

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
        Q = self.cavity_Q
        nharm = self.cavity_harm_num
        wrf = 2*_PI*self.ring.rf_freq
        wr = self.cavity_ang_freq
        if wr == 0:
            raise Exception('wr cannot be zero!')
        if wrf == 0:
            raise Exception('wrf cannot be zero!')
        tan = Q * (wr/(nharm*wrf) - nharm*wrf/wr)
        angle = _np.arctan2(tan, 1)
        return _PI + angle

    @detune_angle.setter
    def detune_angle(self, value):
        Q = self.cavity_Q
        nharm = self.cavity_harm_num
        wrf = 2*_PI*self.ring.rf_freq
        alpha = _np.tan(value)/Q
        self.cavity_ang_freq = nharm*wrf/2
        self.cavity_ang_freq *= (alpha + (alpha**2 + 4)**(1/2))

    def sirius_params(self):
        """."""
        self.ring.en_lost_rad = 871e3
        self.ring.gap_voltage = 3e6
        self.ring.rf_freq = 499.663824e6
        self.ring.harm_num = 864
        self.ring.total_current = 350e-3
        self.ring.energy = 3e9
        self.ring.mom_comp = 1.6446e-4
        self.ring.sync_tune = 4.6520358017681e-3
        self.ring.bunlen = 2.46e-3
        self.ring.espread = 8.4359e-4
        self.ring.dampte = 7.02985e-3

        self.cavity_harm_num = 3

    def maxiv_params(self):
        """."""
        self.ring.en_lost_rad = 856e3
        self.ring.gap_voltage = 1.63e6
        self.ring.rf_freq = 99.93e6
        self.ring.harm_num = 176
        self.ring.total_current = 500e-3
        self.ring.energy = 3e9
        self.ring.mom_comp = 3.07e-4
        self.ring.sync_tune = 1.994e-3
        self.ring.bunlen = 10.1e-3
        self.ring.espread = 7.82e-4
        self.ring.dampte = 25.6e-3

        self.cavity_Q = 21600
        self.cavity_harm_num = 3

    def detune_passive_cavity(self, Rs, k=None):
        """."""
        if k is None:
            kh = self.flat_potential_k_harmonic
        else:
            kh = k
        I0 = self.ring.total_current
        Vrf = self.ring.gap_voltage
        ib = 2*I0*abs(self.form_factor)
        arg = kh*Vrf/(ib*Rs)
        return _np.arccos(arg)

    def print_flat_potential(self):
        """."""
        st = ''
        st += 'form factor                   : {:+8.03f} \n'
        st += 'k flat potential              : {:+08.3f} \n'
        st += 'harmonic phase [deg]          : {:+08.3f} \n'
        st += 'harmonic tuning angle [deg]   : {:+08.3f} \n'
        st += 'detuning frequency    [kHz]   : {:+08.3f} \n'
        st += 'shunt impedance flat [M.ohm]  : {:+08.4f} \n'
        st += 'unperturbed sync. phase [deg] : {:+08.3f} \n'
        st += 'perturbed sync. phase [deg]   : {:+08.3f} \n'
        rad2deg = 180/_PI

        kh = self.flat_potential_k_harmonic
        phih = self.flat_potential_phih_harmonic
        psih = self.flat_potential_psih_harmonic

        # TODO: Since this is a print method we should not set any property
        # here that changes the object state. We need to fix this
        self.detune_angle = psih

        fr = self.cavity_ang_freq/2/_PI
        df = (fr-self.cavity_harm_num*self.ring.rf_freq)*1e-3
        Rs_fp = self.flat_potential_shunt_impedance*1e-6
        phis = self.ring.sync_phase
        new_phis = self.perturbed_sync_phase
        ffactor = self.form_factor
        print(
            st.format(
                ffactor, kh, phih*rad2deg, psih*rad2deg, df, Rs_fp,
                phis*rad2deg, new_phis*rad2deg))

    def integrated_potential(self, z, voltage=None, harmonic=True):
        """."""
        alpha = self.ring.mom_comp
        if voltage is None:
            wrf = 2*_PI*self.ring.rf_freq
            sigmae = self.ring.espread
            phis = self.ring.sync_phase
            h = self.ring.harm_num
            tunes = self.ring.sync_tune
            new_phis = self.perturbed_sync_phase
            nh = self.cavity_harm_num
            phih = self.harmonic_phase
            kh = 0
            if harmonic:
                kh = self.flat_potential_k_harmonic
            phase = wrf*z/_LSPEED
            pot = (alpha**2 * sigmae**2)
            pot /= _np.cos(phis)*(h*alpha*sigmae/tunes)**2
            t1 = _np.cos(new_phis)
            t2 = _np.cos(phase + new_phis)
            t3 = kh/nh * (_np.cos(nh*phih) - _np.cos(nh*phase + nh*phih))
            t4 = (_np.sin(new_phis) + kh*_np.sin(nh*phih))*phase
            pot *= (t1 - t2 + t3 - t4)
        else:
            dpot = (voltage - self.ring.en_lost_rad)/self.ring.energy
            pot = -scy_int.cumtrapz(dpot, z, initial=0)
            pot -= _np.min(pot)
            pot *= alpha/self.ring.circum
        return pot

    def calc_distribution(self, z, voltage=None, harmonic=True):
        """."""
        alpha = self.ring.mom_comp
        espread = self.ring.espread
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
    def calc_complex_form_factor(z, w, rho):
        """."""
        return _np.trapz(
            rho*_np.exp(1j*w*z/_LSPEED), z)/_np.trapz(rho, z)

    def calc_voltages(self, z):
        """."""
        Vrf = self.ring.gap_voltage
        wrf = 2*_PI*self.ring.rf_freq
        phis0 = self.ring.sync_phase
        phis_pert = self.perturbed_sync_phase
        kh = self.flat_potential_k_harmonic
        nh = self.cavity_harm_num
        phih = self.harmonic_phase

        phase = wrf*z/_LSPEED
        Vmain0 = Vrf*_np.sin(phase + phis0)
        Vmain_pert = Vrf*_np.sin(phase + phis_pert)
        Vharm = Vrf*kh*_np.sin(nh*phase + nh*phih)
        return Vmain0, Vmain_pert, Vharm

    def calc_passive_voltage(self, z, Rs, detune_phase, form_factor_phase=0):
        """."""
        I0 = self.ring.total_current
        wrf = 2*_PI*self.ring.rf_freq
        phase = wrf*z/_LSPEED
        nh = self.cavity_harm_num
        ib = 2*I0*abs(self.form_factor)
        volt = -ib*Rs*_np.cos(detune_phase)
        volt *= _np.cos(nh*phase+detune_phase+form_factor_phase)
        return volt

    def loop_form_factor(
            self, z, dist0, kh=None, update_detune=False, full=True):
        """."""
        wrf = self.ring.rf_freq * 2*_PI
        nharm = self.cavity_harm_num
        Rs = self.cavity_shunt_impedance
        Vrf = self.ring.gap_voltage
        self.form_factor = self.calc_complex_form_factor(
            z, nharm*wrf, dist0)

        if update_detune:
            if kh is None:
                k = self.flat_potential_k_harmonic
            else:
                k = kh
            detune = self.detune_passive_cavity(self.cavity_shunt_impedance, k)
            self.detune_angle = detune

        ffang = _np.angle(self.form_factor)
        self.harmonic_phase = (self.detune_angle + _PI/2)/nharm

        vharm0 = self.calc_passive_voltage(
            z=0, Rs=Rs, detune_phase=self.detune_angle,
            form_factor_phase=-ffang)
        arg = (self.ring.en_lost_rad - vharm0)/Vrf
        phis_pert = _PI - _np.arcsin(arg)
        vmain_pert = Vrf * _np.sin(wrf*z/_LSPEED + phis_pert)

        vharm = self.calc_passive_voltage(
            z=z, Rs=Rs, detune_phase=self.detune_angle,
            form_factor_phase=-ffang)
        dist = self.calc_distribution(z=z, voltage=vmain_pert+vharm)

        if full:
            return vmain_pert, vharm, dist
        return dist

    def convergence_distribution(
            self, z, rho0, niter, tol, beta=1, method='anderson_acc',
            kh=None, update_detune=False, print_iter=False):
        """."""
        err_hist = []
        entropy_hist = []
        rho_n = rho0
        X_k = []
        G_k = []

        haissinski = _partial(
            self.loop_form_factor, z, kh=kh, update_detune=update_detune,
            full=False)

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
            idx_nzeros = (rho > 1e-12)
            sint = rho[idx_nzeros]*_np.log(rho[idx_nzeros])
            filt = _np.logical_not(_np.isnan(sint))
            entropy = -_np.trapz(sint[filt], z[idx_nzeros][filt])
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
                if abs(err) < tol:
                    break
            rho0 = rho_n
        return rho, _np.array(err_hist), _np.array(entropy_hist), nit

    def calc_robinson_growth_rate(self, w, wr, approx=False):
        """."""
        alpha = self.ring.mom_comp
        I0 = self.ring.total_current
        E0 = self.ring.energy
        Rs = self.cavity_shunt_impedance
        Q = self.cavity_Q
        w0 = self.ring.rev_ang_freq
        ws = self.ring.sync_tune * w0
        wp = w + ws
        wn = w - ws
        const = I0*alpha*w0/(4*_PI*ws*E0)
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
        return growth - 1/self.ring.dampte

    def calc_tuneshifts_cbi(self, w, wr, m=1, nbun_fill=None, radiation=False):
        """."""
        ring = self.ring
        num_bun = ring.num_bun
        dampte = ring.dampte

        ring.num_bun = nbun_fill if nbun_fill is not None else ring.harm_num
        if not radiation:
            ring.dampte = _np.inf

        Rs = self.cavity_shunt_impedance
        Q = self.cavity_Q
        Zl = imp.longitudinal_resonator(Rs=Rs, Q=Q, wr=wr, w=w)

        deltaw, wp, interpol_Z, spectrum = ring.longitudinal_cbi(
            w=w, Zl=Zl, m=m, inverse=False, full=True)

        # Relative tune-shifts must be multiplied by ws
        deltaw *= (ring.sync_tune * ring.rev_ang_freq)

        ring.num_bun = num_bun
        ring.dampte = dampte

        return deltaw, Zl, wp, interpol_Z, spectrum

    # ---------------------- Make Some Plots ---------------------------------
    def plot_voltages(self, z, vmain, vharm, vtotal):
        """."""
        fig = plt.figure(figsize=(8, 6))
        gs = mpl_gs.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0, 0])

        ax1.plot(z*100, vmain * 1e-6, label=r'$V_{rf}$ Main', color='C0')
        ax1.plot(z*100, vharm * 1e-6, label=r'$V_{3h}$ 3HC', color='C1')
        ax1.plot(
            z*100, vtotal * 1e-6, label=r'$V_{rf} + V_{3h}$ Total', color='C3')
        ax1.axhline(
            y=self.ring.en_lost_rad*1e-6, color='tab:gray',
            ls='--', label=r'$U_0$')

        ax1.set_xlabel(r'$z$ [cm]')
        ax1.set_ylabel('RF voltage [MV]')
        ax1.legend(loc='upper right', fontsize='x-small')
        ax1.grid(True, ls='--', alpha=0.5)
        fig.tight_layout()
        return fig
