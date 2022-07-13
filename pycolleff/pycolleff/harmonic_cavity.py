"""."""
import time as _time
import numpy as _np
import scipy.integrate as scy_int

import mathphys

import pycolleff.impedances as imp
import pycolleff.sirius as si

_LSPEED = mathphys.constants.light_speed
_PI = _np.pi


class HarmonicCavity:
    """."""

    def __init__(self, rf_freq=None):
        """."""
        self.cavity_harm_num = 3
        self.cavity_Q = 4.0e8
        self.cavity_shunt_impedance = 90*self.cavity_Q
        self.rf_freq = rf_freq or si.create_ring().rf_freq
        self.wrf = 2*_PI*self.rf_freq
        self.cavity_ang_freq = self.cavity_harm_num*self.wrf

    def longitudinal_impedance(self, w):
        """."""
        return imp.longitudinal_resonator(
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
            self, ring=None,
            resonators=None, fillpattern=_np.array([])):
        """."""
        self.ring = ring or si.create_ring()
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
        pot = -scy_int.cumtrapz(
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
            w0 = self.ring.rev_ang_freq
            p = _np.arange(0, self.max_mode)
            w = p*w0

        h = self.ring.harm_num
        max_harm = self.max_mode//h
        # zl_wp = self.harmcav.longitudinal_impedance(w=p*w0)
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
            w0 = self.ring.rev_ang_freq
            p = _np.arange(0, self.max_mode)
            w = p*w0
        total_zl = _np.zeros(w.shape, dtype=_np.complex)
        for reson in self.resonators:
            total_zl += reson.longitudinal_impedance(w=w)
        return total_zl

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
            voltage += 2 * zl_wp[None, :] * beam_part[:, None]
            # impedance has dim (zgrid, )
            # beam part has dim (h, )
            # voltage must have dim (h, zgrid)
            # sum over positive frequencies only
        t0b = _time.time() - t0a
        print(f'     E.T. to sum {ps.size:02d} harmonics: {t0b:.3f}s ')
        return -voltage.real

    def calc_longitudinal_equilibrium(self, niter=100, tol=1e-10, beta=1):
        """."""
        dists = [self.dist, ]
        dists = self.anderson_acceleration(dists, niter, tol, beta=beta)
        dists = [self._reshape_dist(rho) for rho in dists]
        self.dist = dists[-1]
        return dists

    def anderson_acceleration(self, dists, niter, tol, m=3, beta=1):
        """."""
        x0 = dists[-1].ravel()
        x = _np.array([x0, self._ffunc(x0)])
        g = _np.array([x[1] - x[0], self._ffunc(x[1]) - x[1]])
        G_k = _np.array(g[1] - g[0], ndmin=2).T
        X_k = _np.array(x[1] - x[0], ndmin=2).T

        for k in range(niter):
            t0 = _time.time()
            m_k = min(k+1, m)
            xnew = []
            gamma_k = _np.linalg.lstsq(
                G_k, g[k], rcond=None)[0]
            mat = X_k + G_k

            xnewi_p1 = beta*(x[k] + g[k] - mat @ gamma_k)
            if beta != 1:
                xnewi_p2 = (1-beta)*(x[k] - X_k @ gamma_k)
                xnew.append(xnewi_p1 + xnewi_p2)
            else:
                xnew.append(xnewi_p1)
            xnew = _np.array(xnew)
            x = _np.r_[x, xnew]
            gnew = _np.array([self._ffunc(x[-1])-x[-1]])
            g = _np.r_[g, gnew]
            G_k = _np.hstack([G_k, (g[-1]-g[-2])[:, None]])
            X_k = _np.hstack([X_k, (x[-1]-x[-2])[:, None]])

            self.dist = self._reshape_dist(x[-1])

            if X_k.shape[1] > m_k:
                X_k = X_k[:, -m:]
                G_k = G_k[:, -m:]

            diff = x[-1]-x[-2]
            diff = self._reshape_dist(diff)
            diff = _np.sqrt(_np.trapz(diff*diff, self.zgrid, axis=1))
            diff = _np.max(diff)
            tf = _time.time() - t0
            print(
                f'Iter.: {k+1:03d}, Dist. Diff.: {diff:.3e}, E.T.: {tf:.3f}s')
            if diff < tol:
                print('distribution ok!')
                break
        return x

    def _ffunc(self, xk):
        """Haissinski operator."""
        xk = self._reshape_dist(xk)
        hvolt = self.voltage_harmonic_cavity_impedance(dist=xk)
        tvolt = self.main_voltage[None, :] + hvolt
        fxk = self.calc_distribution_from_voltage(tvolt)
        return fxk.ravel()

    def _reshape_dist(self, dist):
        return dist.reshape((self.ring.harm_num, self.zgrid.size))
