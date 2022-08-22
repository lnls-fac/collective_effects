"""."""
import time as _time

import numpy as _np
import scipy.integrate as _scy_int

from mathphys.constants import light_speed as _LSPEED

from . import impedances as _imp
from .colleff import Ring as _Ring

_PI = _np.pi


class Resonator:
    """."""

    def __init__(self, Rs=0, Q=0, ang_freq=0, harm_rf=3):
        """."""
        self.ang_freq = ang_freq
        self.Q = Q
        self.shunt_impedance = Rs

        self.harm_rf = harm_rf
        self.ang_freq_rf = 0

    def get_impedance(self, w):
        """."""
        return _imp.longitudinal_resonator(
            Rs=self.shunt_impedance, Q=self.Q, wr=self.ang_freq, w=w)

    @property
    def RoverQ(self):
        """."""
        return self.shunt_impedance/self.Q

    @property
    def detune_w(self):
        """."""
        return self.ang_freq - self.harm_rf*self.ang_freq_rf

    @detune_w.setter
    def detune_w(self, value):
        """."""
        wr = self.harm_rf*self.ang_freq_rf + value
        self.ang_freq = wr

    @property
    def detune_angle(self):
        """."""
        Q = self.Q
        nharm = self.harm_rf
        wrf = self.ang_freq_rf
        wr = self.ang_freq
        if wr == 0:
            raise Exception('wr cannot be zero!')
        if wrf == 0:
            raise Exception('wrf cannot be zero!')
        tan = Q * (wr/(nharm*wrf) - nharm*wrf/wr)
        angle = _np.arctan2(tan, 1)
        return _PI + angle

    @detune_angle.setter
    def detune_angle(self, value):
        Q = self.Q
        nharm = self.harm_rf
        wrf = self.ang_freq_rf

        delta = _np.tan(value)/2/Q
        self.ang_freq = nharm*wrf*(delta + (1+delta**2)**(1/2))

    def to_dict(self):
        """Save state to dictionary."""
        return dict(
            ang_freq=self.ang_freq,
            Q=self.Q,
            shunt_impedance=self.shunt_impedance,
            harm_rf=self.harm_rf,
            ang_freq_rf=self.ang_freq_rf)

    def from_dict(self, dic):
        """Load state from dictionary."""
        self.ang_freq = dic.get('ang_freq', self.ang_freq)
        self.Q = dic.get('Q', self.Q)
        self.shunt_impedance = dic.get('shunt_impedance', self.shunt_impedance)
        self.harm_rf = dic.get('harm_rf', self.harm_rf)
        self.ang_freq_rf = dic.get('ang_freq_rf', self.ang_freq_rf)


class LongitudinalEquilibrium:
    """."""

    def __init__(
            self, ring: _Ring, resonators: list, fillpattern=None):
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
        self.main_voltage = self.ring.get_voltage_waveform(self._zgrid)

    @property
    def main_voltage(self):
        """."""
        return self._main_voltage

    @main_voltage.setter
    def main_voltage(self, value):
        if value.shape[-1] != self._zgrid.shape[0]:
            raise ValueError('Wrong shape for voltage.')
        self._main_voltage = value
        self.distributions = self.calc_distributions_from_voltage(
            self._main_voltage)

    @property
    def fillpattern(self):
        """."""
        return self._fillpattern

    @fillpattern.setter
    def fillpattern(self, value):
        if value.size != self.ring.harm_num:
            raise ValueError('Wrong size for fillparttern.')
        self._fillpattern = value

    @property
    def filled_buckets(self):
        """."""
        idx = _np.where(self.fillpattern != 0)[0]
        return idx

    @property
    def distributions(self):
        """."""
        return self._dist

    @distributions.setter
    def distributions(self, value):
        if value.ndim != 2:
            raise ValueError('Distributions must have 2 dimensions.')
        elif value.shape[0] != self.ring.harm_num:
            raise ValueError('First dimension must be equal ring.harm_num.')
        elif value.shape[1] != self._zgrid.size:
            raise ValueError('Second dimension must be equal zgrid.size.')
        self._dist = value

    def to_dict(self):
        """Save state to dictionary."""
        return dict(
            ring=self.ring.to_dict(),
            resonators=[res.to_dict() for res in self.resonators],
            zgrid=self._zgrid,
            dist=self._dist,
            fillpatern=self._fillpattern,
            main_voltage=self._main_voltage,
            max_mode=self.max_mode,
            min_mode0_ratio=self.min_mode0_ratio)

    def from_dict(self, dic):
        """Load state from dictionary."""
        self.ring.from_dict(dic.get('ring', dict()))
        resons = []
        for res in dic.get('resonators', self.resonators):
            re = Resonator()
            re.from_dict(res)
            resons.append(re)
        self.resonators = resons
        self._zgrid = dic.get('zgrid', self._zgrid)
        self._dist = dic.get('dist', self._dist)
        self._fillpattern = dic.get('fillpattern', self._fillpattern)
        self._main_voltage = dic.get('main_voltage', self._main_voltage)
        self.max_mode = dic.get('max_mode', self.max_mode)
        self.min_mode0_ratio = dic.get('min_mode0_ratio', self.min_mode0_ratio)

    def create_zgrid(self, nr_points=1001, sigmas=30):
        """."""
        return sigmas*self.ring.bunlen*_np.linspace(-1, 1, nr_points)

    def calc_moments(self, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
        zm = _np.trapz(self.zgrid[None, :]*dist, self.zgrid, axis=1)
        zgrid2 = self.zgrid**2
        z2 = _np.trapz(zgrid2[None, :]*dist, self.zgrid, axis=1)
        return zm, _np.sqrt(z2 - zm**2)

    def get_gaussian_distributions(self, sigmaz, z0=0):
        """."""
        arg = (self.zgrid - z0)/sigmaz
        dist = _np.exp(-arg**2/2)
        norm = _np.trapz(dist, self.zgrid)
        dist = dist/norm
        dist = _np.tile(dist, (self.ring.harm_num, 1))
        return dist

    def calc_harmonic_voltage_for_flat_potential(self, harm_rf=3):
        """."""
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        kharm = 1/harm_rf**2 - ((U0/Vrf)**2)/(harm_rf**2-1)
        return kharm**(1/2)

    def calc_detune_for_fixed_harmonic_voltage(
            self, peak_harm_volt, harm_rf=3, Rs=0):
        """."""
        I0 = _np.sum(self.fillpattern)
        # TODO: This way of including the form factor is temporary. Fix it.
        wr = 2 * _PI * self.ring.rf_freq * harm_rf
        form_factor = self.calc_fourier_transform(wr)
        ib = 2 * I0 * _np.abs(form_factor).mean()
        arg = peak_harm_volt/ib/Rs
        return _np.arccos(arg)

    def calc_harmonic_voltage_for_fixed_detune(self, detune, harm_rf=3, Rs=0):
        """."""
        I0 = _np.sum(self.fillpattern)
        # TODO: This way of including the form factor is temporary. Fix it.
        wr = 2 * _PI * self.ring.rf_freq * harm_rf
        form_factor = self.calc_fourier_transform(wr)
        ib = 2 * I0 * _np.abs(form_factor).mean()
        peak_harm_volt = Rs * ib * _np.cos(detune)
        return _np.abs(peak_harm_volt)

    def calc_distributions_from_voltage(self, voltage):
        """."""
        flag = False
        if len(voltage.shape) < 2:
            flag = True
            # voltage must be (h, zgrid) or (1, zgrid)
            voltage = voltage[None, :]

        # subtract U0
        U0 = self.ring.en_lost_rad
        pot = -_scy_int.cumtrapz(voltage-U0, self.zgrid, initial=0)

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

    def calc_fourier_transform(self, w, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
        arg = dist*_np.exp(1j*w*self.zgrid/_LSPEED)[None, :]
        return _np.trapz(arg, self.zgrid, axis=1)

    def get_impedance(self, w=None):
        """."""
        if w is None:
            w = self._create_freqs()
        total_zl = _np.zeros(w.shape, dtype=_np.complex)
        for reson in self.resonators:
            total_zl += reson.get_impedance(w=w)
        return total_zl

    def get_harmonics_impedance_and_filling(self, w=None):
        """."""
        if w is None:
            w = self._create_freqs()
        h = self.ring.harm_num
        max_harm = self.max_mode//h
        diff = self.max_mode - max_harm*h
        if diff > 0:
            max_harm = max_harm + 1
        zl_wp = self.get_impedance(w=w)
        fill_fft = _np.fft.fft(self.fillpattern)
        fill_fft = _np.tile(fill_fft, (max_harm, 1)).ravel()
        if diff > 0:
            fill_fft = fill_fft[:-(h-diff)]
        zl_fill = _np.abs(zl_wp * fill_fft)
        modes = _np.where(
            zl_fill > zl_fill.max()*self.min_mode0_ratio)[0]
        return modes, zl_wp[modes]

    def calc_voltage_harmonic_cavity_impedance(self, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
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
            dist_fourier = self.calc_fourier_transform(w=wp, dist=dist)

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
        dists = [self.distributions, ]
        dists = self._apply_anderson_acceleration(
            dists, niter, tol, beta=beta, m=m)
        dists = [self._reshape_dist(rho) for rho in dists]
        self.distributions = dists[-1]
        return dists

    def calc_robinson_growth_rate(
            self, w, approx=False, wr=None, Rs=None, Q=None):
        """."""
        alpha = self.ring.mom_comp
        I0 = self.ring.total_current
        E0 = self.ring.energy
        w0 = self.ring.rev_ang_freq
        ws = self.ring.sync_tune * w0
        wp = w + ws
        wn = w - ws
        const = I0*alpha*w0/(4*_PI*ws*E0)
        if approx and None not in {wr, Rs, Q}:
            x = w/wr
            const_approx = const*4*ws
            growth = const_approx
            growth *= Rs*Q**2
            growth *= (1-x**2)*(1 + x**2)
            growth /= x**4 * (1 + Q**2 * (1/x - x)**2)**2
        else:
            Zlp = self.get_impedance(w=wp)
            Zln = self.get_impedance(w=wn)
            growth = const*(wp*Zlp.real - wn*Zln.real)
        return growth

    def calc_tuneshifts_cbi(self, w, m=1, nbun_fill=None, radiation=False):
        """."""
        ring = self.ring
        num_bun = ring.num_bun
        dampte = ring.dampte

        ring.num_bun = nbun_fill if nbun_fill is not None else ring.harm_num
        if not radiation:
            ring.dampte = _np.inf

        total_zl = self.get_impedance(w=w)

        deltaw, wp, interpol_Z, spectrum = ring.longitudinal_cbi(
            w=w, Zl=total_zl, m=m, inverse=False, full=True)

        # Relative tune-shifts must be multiplied by ws
        deltaw *= (ring.sync_tune * ring.rev_ang_freq)

        ring.num_bun = num_bun
        ring.dampte = dampte
        return deltaw, total_zl, wp, interpol_Z, spectrum

    # -------------------- auxiliary methods ----------------------------------
    def _apply_anderson_acceleration(self, dists, niter, tol, m=3, beta=1):
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

    def _create_freqs(self):
        w0 = self.ring.rev_ang_freq
        p = _np.arange(0, self.max_mode)
        return p*w0

    def _ffunc(self, xk):
        """Haissinski operator."""
        xk = self._reshape_dist(xk)
        hvolt = self.calc_voltage_harmonic_cavity_impedance(dist=xk)
        tvolt = self.main_voltage[None, :] + hvolt
        fxk = self.calc_distributions_from_voltage(tvolt)
        return fxk.ravel()

    def _reshape_dist(self, dist):
        return dist.reshape((self.ring.harm_num, self.zgrid.size))
