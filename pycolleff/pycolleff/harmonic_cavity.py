"""."""
import time as _time
import numpy as _np
import scipy.integrate as _scy_int

import mathphys as _mathphys

import pycolleff.impedances as _imp

_LSPEED = _mathphys.constants.light_speed
_PI = _np.pi


class HarmonicCavity:
    """."""

    def __init__(self):
        """."""
        self.cavity_harm_num = None
        self.cavity_Q = None
        self.cavity_shunt_impedance = None

        # SSRF 3HC
        # self.cavity_Q = 4.0e8
        # self.cavity_shunt_impedance = 90*self.cavity_Q
        # self.rf_freq = rf_freq or si.create_ring().rf_freq

        # MAX-IV 3HC
        # self.cavity_Q = 20800
        # self.cavity_shunt_impedance = 8.25e6
        # self.rf_freq = rf_freq or maxiv.create_ring().rf_freq

        self.wrf = None
        self.cavity_ang_freq = None

    def longitudinal_impedance(self, w):
        """."""
        return _imp.longitudinal_resonator(
                Rs=self.cavity_shunt_impedance,
                Q=self.cavity_Q,
                wr=self.cavity_ang_freq, w=w)

    @property
    def RoverQ(self):
        """."""
        return self.cavity_shunt_impedance/self.cavity_Q

    @property
    def detune_w(self):
        """."""
        return self.cavity_ang_freq-self.cavity_harm_num*self.wrf

    @detune_w.setter
    def detune_w(self, value):
        """."""
        wrf = 2*_PI*self.rf_freq
        wr = self.cavity_harm_num*wrf + value
        self.cavity_ang_freq = wr

    @property
    def detune_angle(self):
        """."""
        Q = self.cavity_Q
        nharm = self.cavity_harm_num
        wrf = 2*_PI*self.rf_freq
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
        wrf = 2*_PI*self.rf_freq

        delta = _np.tan(value)/2/Q
        self.cavity_ang_freq = nharm*wrf*(delta + (1+delta**2)**(1/2))


class LongitudinalEquilibrium:
    """."""

    def __init__(
            self, ring,
            resonators, fillpattern=_np.array([])):
        """."""
        self.ring = ring
        self.resonators = resonators
        self._zgrid = None
        self._dist = None
        self._fillpattern = None
        self._main_voltage = None
        self.max_mode = 10*self.ring.harm_num
        self.min_mode0_ratio = 1e-9
        self.fillpattern = fillpattern
        self.zgrid = self.create_zgrid()
        self.main_voltage = self.voltage_main_cavity()

    @property
    def max_mode(self):
        """."""
        return self._max_mode

    @max_mode.setter
    def max_mode(self, value):
        if value % self.ring.harm_num:
            raise Exception('mode must be a multiple of harmonic number!')
        else:
            self._max_mode = value

    @property
    def zgrid(self):
        """."""
        return self._zgrid

    @zgrid.setter
    def zgrid(self, value):
        self._zgrid = value
        main_voltage = self.voltage_main_cavity()
        self.main_voltage = main_voltage
        self.dist = self.calc_distribution_from_voltage(main_voltage)

    @property
    def main_voltage(self):
        """."""
        return self._main_voltage

    @main_voltage.setter
    def main_voltage(self, value):
        self._main_voltage = value

    def to_dict(self):
        """."""
        return dict(
            energy=self.ring.energy,
            en_lost_rad=self.ring.en_lost_rad,
            harm_num=self.ring.harm_num,
            rf_freq=self.ring.rf_freq,
            espread=self.ring.espread,
            mom_comp=self.ring.mom_comp,
            gap_voltage=self.ring.gap_voltage,
            sync_tune=self.ring.sync_tune,
            )

    def from_dict(self, dic):
        """."""
        self.energy = dic['energy']
        self.en_lost_rad = dic['en_lost_rad']
        self.harm_num = dic['harm_num']
        self.rf_freq = dic['rf_freq']
        self.espread = dic['espread']
        self.mom_comp = dic['mom_comp']
        self.gap_voltage = dic['gap_voltage']
        self.sync_tune = dic['sync_tune']

    def create_zgrid(self, nr_points=1001, sigmas=30):
        """."""
        return sigmas*self.ring.bunlen*_np.linspace(-1, 1, nr_points)

    def calc_first_moment(self, dist=None):
        """."""
        if dist is None:
            dist = self.dist
        return _np.trapz(self.zgrid[None, :]*dist, self.zgrid, axis=1)

    def calc_second_moment(self, dist=None):
        """."""
        if dist is None:
            dist = self.dist
        zm = self.calc_first_moment(dist)
        zgrid2 = self.zgrid**2
        z2 = _np.trapz(zgrid2[None, :]*dist, self.zgrid, axis=1)
        return _np.sqrt(z2 - zm**2)

    def gaussian_distribution(self, sigmaz, z0=0):
        """."""
        arg = (self.zgrid - z0)/sigmaz
        dist = _np.exp(-arg**2/2)
        norm = _np.trapz(dist, self.zgrid)
        dist = dist/norm
        dist = _np.tile(dist, (self.ring.harm_num, 1))
        return dist

    @property
    def fillpattern(self):
        """."""
        return self._fillpattern

    @fillpattern.setter
    def fillpattern(self, value):
        self._fillpattern = value

    @property
    def filled_buckets(self):
        """."""
        idx = _np.where(self.fillpattern != 0)[0]
        return idx

    @property
    def dist(self):
        """."""
        return self._dist

    @dist.setter
    def dist(self, value):
        self._dist = value

    def harmonic_voltage_for_flat_potential(self):
        """."""
        nharm = self.resonators[0].cavity_harm_num
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        kharm = 1/nharm**2 - ((U0/Vrf)**2)/(nharm**2-1)
        return kharm**(1/2)

    def detune_for_fixed_harmonic_voltage(self, peak_harm_volt):
        """."""
        Rs = self.resonators[0].cavity_shunt_impedance
        I0 = _np.sum(self.fillpattern)
        # <<< WARNING >>>
        # This way of including the form factor is temporary
        ib = 2*I0*_np.mean(_np.abs(self.form_factor))
        arg = peak_harm_volt/ib/Rs
        return _np.arccos(arg)

    def harmonic_voltage_for_fixed_detune(self, detune):
        """."""
        Rs = self.resonators[0].cavity_shunt_impedance
        I0 = _np.sum(self.fillpattern)
        # <<< WARNING >>>
        # This way of including the form factor is temporary
        ib = 2*I0*_np.mean(_np.abs(self.form_factor))
        peak_harm_volt = Rs*ib*_np.cos(detune)
        return _np.abs(peak_harm_volt)

    @property
    def form_factor(self):
        """."""
        wrf = 2*_PI*self.ring.rf_freq
        w = self.resonators[0].cavity_harm_num*wrf
        arg = self.dist*_np.exp(1j*w*self.zgrid/_LSPEED)[None, :]
        factor = _np.trapz(arg, self.zgrid, axis=1)
        norm = _np.trapz(self.dist, self.zgrid, axis=1)
        return factor/norm

    def calc_distribution_from_voltage(self, voltage):
        """."""
        flag = False
        if len(voltage.shape) < 2:
            flag = True
            # voltage must be (h, zgrid) or (1, zgrid)
            voltage = voltage[None, :]

        # subtract U0
        U0 = self.ring.en_lost_rad
        pot = -_scy_int.cumtrapz(
            voltage-U0, self.zgrid, initial=0)

        # subtract minimum value for all bunches
        pot -= _np.min(pot, axis=1)[:, None]

        const = self.ring.espread**2
        const *= self.ring.mom_comp
        const *= self.ring.circum
        # normalize by E0
        const *= self.ring.energy

        dist = _np.exp(-pot/const)
        norm = _np.trapz(dist, self.zgrid, axis=1)
        if flag:
            dist = _np.tile(dist, (self.ring.harm_num, 1))
        # distribution must be normalized
        return dist/norm[:, None]

    def distribution_fourier_transform(self, w, dist=None):
        """."""
        if dist is None:
            dist = self.dist
        arg = dist*_np.exp(1j*w*self.zgrid/_LSPEED)[None, :]
        return _np.trapz(arg, self.zgrid, axis=1)

    def voltage_main_cavity(self, z=None, sync_phase=None):
        """."""
        wrf = 2*_PI*self.ring.rf_freq
        phase0 = sync_phase or self.ring.sync_phase
        if z is None:
            z = self.zgrid
        phase = wrf*z/_LSPEED
        phase += phase0
        voltage = self.ring.gap_voltage*_np.sin(phase)
        return voltage

    def get_harmonics_impedance_and_filling(self, w=None):
        """."""
        if w is None:
            w = self._create_freqs()
        h = self.ring.harm_num
        max_harm = self.max_mode//h
        zl_wp = self.get_total_impedance(w=w)
        fill_fft = _np.abs(_np.fft.fft(self.fillpattern))
        fill_fft = _np.tile(fill_fft, (max_harm, 1)).ravel()
        zl_fill = (zl_wp * fill_fft).real
        modes = _np.where(
            zl_fill > zl_fill.max()*self.min_mode0_ratio)[0]
        return modes, zl_wp[modes]

    def get_total_impedance(self, w=None):
        """."""
        if w is None:
            w = self._create_freqs()
        total_zl = _np.zeros(w.shape, dtype=_np.complex)
        for reson in self.resonators:
            total_zl += reson.longitudinal_impedance(w=w)
        return total_zl

    def _create_freqs(self):
        w0 = self.ring.rev_ang_freq
        p = _np.arange(0, self.max_mode)
        return p*w0

    def voltage_harmonic_cavity_impedance(self, dist=None):
        """."""
        if dist is None:
            dist = self.dist
        w0 = self.ring.rev_ang_freq
        voltage = _np.zeros(
            (self.ring.harm_num, self.zgrid.size), dtype=_np.complex)

        ind = _np.arange(self.ring.harm_num)
        zn_ph = 2*_PI/self.ring.harm_num*ind
        z_ph = w0*self.zgrid/_LSPEED

        ps, zl_wps = self.get_harmonics_impedance_and_filling()
        t0a = _time.time()
        for idx, p in enumerate(ps):
            wp = p * w0
            dist_fourier = self.distribution_fourier_transform(w=wp, dist=dist)

            exp_phase = _np.exp(-1j*p*zn_ph)
            beam_part = exp_phase * self.fillpattern * dist_fourier.conj()
            beam_part = _np.sum(beam_part)
            beam_part = beam_part/exp_phase
            zl_wp = zl_wps[idx].conj() * _np.exp(1j*p*z_ph)
            # sum over positive frequencies only -> factor 2
            voltage += 2 * zl_wp[None, :] * beam_part[:, None]
        t0b = _time.time() - t0a
        print(f'     E.T. to sum {ps.size:02d} harmonics: {t0b:.3f}s ')
        return -voltage.real

    def calc_longitudinal_equilibrium(self, niter=100, tol=1e-10, beta=1, m=3):
        """."""
        dists = [self.dist, ]
        dists = self.anderson_acceleration(dists, niter, tol, beta=beta, m=m)
        dists = [self._reshape_dist(rho) for rho in dists]
        self.dist = dists[-1]
        return dists

    def anderson_acceleration(self, dists, niter, tol, m=3, beta=1):
        """."""
        if beta < 0:
            raise Exception('relaxation parameter beta must be positive.')

        xold = dists[-1].ravel()
        xnew = self._ffunc(xold)
        dists.append(xnew)
        g = _np.array([xnew - xold, self._ffunc(xnew) - xnew])
        G_k = _np.array(g[-1] - g[-2], ndmin=2).T
        X_k = _np.array(g[0], ndmin=2).T

        for k in range(niter):
            t0 = _time.time()
            m_k = min(k+1, m)
            gamma_k = _np.linalg.lstsq(G_k, g[-1], rcond=None)[0]
            mat = X_k + G_k
            xold = xnew
            xnew = beta * (xold + g[-1] - mat @ gamma_k)
            if beta != 1:
                xnew += (1-beta) * (xold - X_k @ gamma_k)
            dists.append(xnew)
            gnew = self._ffunc(xnew)-xnew
            g = _np.vstack([g[-1], gnew])
            G_k = _np.hstack([G_k, (g[-1]-g[-2])[:, None]])
            X_k = _np.hstack([X_k, (xnew-xold)[:, None]])

            if X_k.shape[1] > m_k:
                # only m previous iterations are needed
                X_k = X_k[:, -m:]
                G_k = G_k[:, -m:]

            diff = self._reshape_dist(gnew)
            diff = _np.trapz(_np.abs(diff), self.zgrid, axis=1)
            idx = _np.argmax(diff)
            tf = _time.time() - t0
            print(
                f'Iter.: {k+1:03d}, Dist. Diff.: {diff[idx]:.3e}' +
                f' (bucket {idx:03d}), E.T.: {tf:.3f}s')
            if diff[idx] < tol:
                print('distribution ok!')
                break
        return dists

    def _ffunc(self, xk):
        """Haissinski operator."""
        xk = self._reshape_dist(xk)
        hvolt = self.voltage_harmonic_cavity_impedance(dist=xk)
        tvolt = self.main_voltage[None, :] + hvolt
        fxk = self.calc_distribution_from_voltage(tvolt)
        return fxk.ravel()

    def _reshape_dist(self, dist):
        return dist.reshape((self.ring.harm_num, self.zgrid.size))

    def calc_robinson_growth_rate(self, w, wr, approx=False):
        """."""
        alpha = self.ring.mom_comp
        I0 = self.ring.total_current
        E0 = self.ring.energy
        Rs = self.resonators[0].cavity_shunt_impedance
        Q = self.resonators[0].cavity_Q
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
            Zlp = self.resonators[0].longitudinal_impedance(w=wp)
            Zln = self.resonators[0].longitudinal_impedance(w=wn)
            growth = const*(wp*Zlp.real - wn*Zln.real)
        return growth - 1/self.ring.dampte

    def calc_tuneshifts_cbi(self, w, m=1, nbun_fill=None, radiation=False):
        """."""
        ring = self.ring
        num_bun = ring.num_bun
        dampte = ring.dampte

        ring.num_bun = nbun_fill if nbun_fill is not None else ring.harm_num
        if not radiation:
            ring.dampte = _np.inf

        Zl = self.resonators[0].longitudinal_impedance(w=w)

        deltaw, wp, interpol_Z, spectrum = ring.longitudinal_cbi(
            w=w, Zl=Zl, m=m, inverse=False, full=True)

        # Relative tune-shifts must be multiplied by ws
        deltaw *= (ring.sync_tune * ring.rev_ang_freq)

        ring.num_bun = num_bun
        ring.dampte = dampte

        return deltaw, Zl, wp, interpol_Z, spectrum
