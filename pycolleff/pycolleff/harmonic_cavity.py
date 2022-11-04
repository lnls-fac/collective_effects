"""."""
import time as _time
import numpy as _np
import numexpr as _ne

from mathphys.constants import light_speed as _LSPEED
from mathphys.functions import get_namedtuple as _get_namedtuple

from scipy.optimize import least_squares as _least_squares
from . import impedances as _imp
from .colleff import Ring as _Ring

_PI = _np.pi


def _mytrapz(y, dx, cumul=False):
    """Perform trapezoidal integration along last axis of array.

    Args:
        y (numpy.ndarray, (..., N)): array where integration is performed.
        dx (float): step size.
        cumul (bool, optional): Whether or not to return cummulative integral.
            Defaults to False.

    Returns:
        numpy.ndarray: if cumul is True, then the shape matches y, else the
            number of dimensions is reduced.

    """
    y1 = y[..., :-1]
    y2 = y[..., 1:]
    if cumul:
        intg = _np.zeros_like(y)
        intg[..., 1:] = _ne.evaluate('(y1 + y2)*dx/2.0')
        res = _np.cumsum(intg, axis=-1)
        res.shape = intg.shape
        return res
    else:
        return _ne.evaluate('(y1 + y2)*dx/2.0').sum(axis=-1)


class ImpedanceSource:
    """."""

    Methods = _get_namedtuple('Methods', ['Impedance', 'Wake'])
    ActivePassive = _get_namedtuple('ActivePassive', ['Active', 'Passive'])

    def __init__(
            self, Rs=0, Q=0, ang_freq=0, harm_rf=3,
            calc_method=Methods.Wake, active_passive=ActivePassive.Passive):
        """."""
        self._calc_method = None
        self._active_passive = None

        self.ang_freq = ang_freq
        self.Q = Q
        self.shunt_impedance = Rs

        self.harm_rf = harm_rf
        self.ang_freq_rf = 0
        self._loop_ctrl_freq = 0
        self.calc_method = calc_method
        self.active_passive = active_passive

    @property
    def calc_method_str(self):
        """."""
        return self.Methods._fields[self._calc_method]

    @property
    def calc_method(self):
        """."""
        return self._calc_method

    @calc_method.setter
    def calc_method(self, value):
        if value is None:
            return
        if isinstance(value, str):
            self._calc_method = self.Methods._fields.index(value)
        elif int(value) in self.Methods:
            self._calc_method = int(value)
        else:
            raise ValueError('Wrong value for calc_method.')

    @property
    def active_passive_str(self):
        """."""
        return self.ActivePassive._fields[self._active_passive]

    @property
    def active_passive(self):
        """."""
        return self._active_passive

    @active_passive.setter
    def active_passive(self, value):
        if value is None:
            return
        if isinstance(value, str):
            self._active_passive = self.ActivePassive._fields.index(value)
        elif int(value) in self.ActivePassive:
            self._active_passive = int(value)
        else:
            raise ValueError('Wrong value for active_passive.')

    def get_impedance(self, w):
        """."""
        imp = _imp.longitudinal_resonator(
            Rs=self.shunt_impedance, Q=self.Q, wr=self.ang_freq, w=w)
        return imp

    @property
    def loop_ctrl_freq(self):
        """."""
        return self._loop_ctrl_freq

    @loop_ctrl_freq.setter
    def loop_ctrl_freq(self, value):
        """."""
        self._loop_ctrl_freq = value

    @property
    def loop_ctrl_ang_freq(self):
        """."""
        return 2*_PI*self._loop_ctrl_freq

    @loop_ctrl_ang_freq.setter
    def loop_ctrl_ang_freq(self, value):
        """."""
        self._loop_ctrl_freq = value/2/_PI

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
    def detune_freq(self):
        """."""
        return self.detune_w/2/_PI

    @property
    def alpha(self):
        """."""
        return self.ang_freq/2/self.Q

    @property
    def ang_freq_bar(self):
        """."""
        wr_ = self.ang_freq
        alpha = self.alpha
        return (wr_*wr_-alpha*alpha)**0.5

    @property
    def beta(self):
        """."""
        return (self.alpha - 1j*self.ang_freq_bar)/_LSPEED

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

    def __str__(self):
        """."""
        stmp = '{0:20s}: {1:}  {2:s}\n'.format
        ftmp = '{0:20s}: {1:3.2f}  {2:s}\n'.format
        etmp = '{0:20s}: {1:.2e}  {2:s}\n'.format
        mega = 1e-6
        stg = stmp('calc_method', self.calc_method_str, '')
        stg += stmp('active_passive', self.active_passive_str, '')
        stg += ftmp('ang_freq_rf', self.ang_freq_rf*mega, '[MHz]')
        stg += ftmp('ang_freq', self.ang_freq*mega, '[Mrad/s]')
        stg += ftmp('shunt_impedance', self.shunt_impedance*mega, '[MOhm]')
        stg += etmp('Q', self.Q, '')
        stg += ftmp('RoverQ', self.RoverQ, '[Ohm]')
        stg += ftmp('harm_rf', self.harm_rf, '')
        stg += ftmp('detune_angle', self.detune_angle, '[rad]')
        stg += ftmp('detune_freq', self.detune_freq/1e3, '[kHz]')
        stg += ftmp('detune_w', self.detune_w, '[rad/s]')
        stg += ftmp('alpha', self.alpha, '[rad/s]')
        stg += ftmp('ang_freq_bar', self.ang_freq_bar*mega, '[Mrad/s]')
        return stg


class LongitudinalEquilibrium:
    """."""

    def __init__(self, ring: _Ring, impedance_sources: list, fillpattern=None):
        """."""
        self._zgrid = None
        self._dist = None
        self._fillpattern = None
        self._main_voltage = None
        self._calc_fun = None
        self._calc_method = None
        self._print_flag = False
        self._exp_z = None
        self._wake_matrix = None
        self.beamload_active = None

        self.ring = ring
        self.impedance_sources = impedance_sources
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
        self._max_mode = value

    @property
    def zgrid(self):
        """."""
        return self._zgrid

    @zgrid.setter
    def zgrid(self, value):
        self._zgrid = value
        self.main_voltage = self.ring.get_voltage_waveform(self._zgrid)
        self._exp_z = None

    @property
    def main_voltage(self):
        """."""
        return self._main_voltage

    @main_voltage.setter
    def main_voltage(self, value):
        if value.shape[-1] != self._zgrid.shape[0]:
            raise ValueError('Wrong shape for voltage.')
        self._main_voltage = value
        self.distributions, _ = self.calc_distributions_from_voltage(
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
        self._wake_matrix = None

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

    @property
    def print_flag(self):
        """."""
        return self._print_flag

    @print_flag.setter
    def print_flag(self, value):
        self._print_flag = value

    def to_dict(self):
        """Save state to dictionary."""
        return dict(
            ring=self.ring.to_dict(),
            impedance_sources=[
                imp.to_dict() for imp in self.impedance_sources],
            zgrid=self._zgrid,
            dist=self._dist,
            fillpatern=self._fillpattern,
            main_voltage=self._main_voltage,
            max_mode=self.max_mode,
            min_mode0_ratio=self.min_mode0_ratio,
            calc_method=self.calc_method_str)

    def from_dict(self, dic):
        """Load state from dictionary."""
        self.ring.from_dict(dic.get('ring', dict()))
        imps = []
        for imp in dic.get('impedance_sources', self.impedance_sources):
            _imp = ImpedanceSource()
            _imp.from_dict(imp)
            imps.append(_imp)
        self.impedance_sources = imps
        self._zgrid = dic.get('zgrid', self._zgrid)
        self._dist = dic.get('dist', self._dist)
        self._fillpattern = dic.get('fillpattern', self._fillpattern)
        self._main_voltage = dic.get('main_voltage', self._main_voltage)
        self.max_mode = dic.get('max_mode', self.max_mode)
        self.min_mode0_ratio = dic.get('min_mode0_ratio', self.min_mode0_ratio)
        self.calc_method = dic.get('calc_method', self.calc_method)

    def create_zgrid(self, nr_points=1001, sigmas=30):
        """."""
        return sigmas*self.ring.bunlen*_np.linspace(-1, 1, nr_points)

    @staticmethod
    def calc_moments(zgrid, dist):
        """."""
        dz = zgrid[1] - zgrid[0]
        zm = _mytrapz(zgrid[None, :]*dist, dz)

        zgrid2 = zgrid*zgrid
        z2 = _mytrapz(zgrid2[None, :]*dist, dz)
        return zm, _np.sqrt(z2 - zm**2)

    def get_gaussian_distributions(self, sigmaz, z0=0):
        """."""
        dz = self.zgrid[1] - self.zgrid[0]

        arg = (self.zgrid - z0)/sigmaz
        dist = _np.exp(-arg**2/2)
        dist /= _mytrapz(dist, dz)
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
        form_factor = self.calc_fourier_transform(wr)[self.filled_buckets]
        ib = 2 * I0 * _np.abs(form_factor).mean()
        arg = peak_harm_volt/ib/Rs
        return _np.arccos(arg)

    def calc_harmonic_voltage_for_fixed_detune(self, detune, harm_rf=3, Rs=0):
        """."""
        I0 = _np.sum(self.fillpattern)
        # TODO: This way of including the form factor is temporary. Fix it.
        wr = 2 * _PI * self.ring.rf_freq * harm_rf
        form_factor = self.calc_fourier_transform(wr)[self.filled_buckets]
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

        dz = self.zgrid[1] - self.zgrid[0]

        # subtract U0
        U0 = self.ring.en_lost_rad
        pot = -_mytrapz(voltage-U0, dz, cumul=True)

        # subtract minimum value for all bunches
        pot -= _np.min(pot, axis=1)[:, None]

        const = self.ring.espread**2
        const *= self.ring.mom_comp
        const *= self.ring.circum
        # normalize by E0
        const *= self.ring.energy

        dist = _ne.evaluate('exp(-pot/const)')
        dist /= _mytrapz(dist, dz)[:, None]
        if flag:
            dist = _np.tile(dist, (self.ring.harm_num, 1))
        # distribution must be normalized
        return dist, pot

    def calc_fourier_transform(self, w, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
        arg = _np.exp((1j*w/_LSPEED)*self.zgrid)[None, :]
        arg = _ne.evaluate('dist * arg')
        dz = self.zgrid[1] - self.zgrid[0]
        return _mytrapz(arg, dz)

    def get_impedance(self, w=None, apply_filter=False):
        """."""
        if w is None:
            w = self._create_freqs()
        total_zl = _np.zeros(w.shape, dtype=_np.complex)
        imp_idx = self._get_impedance_types_idx()
        loop_freq_idx = []
        for iidx in imp_idx:
            imp = self.impedance_sources[iidx]
            cond = imp.active_passive == ImpedanceSource.ActivePassive.Active
            cond |= apply_filter
            if cond:
                # this part is broken for instabilities calculations
                w_loop = imp.loop_ctrl_ang_freq
                if w_loop not in w:
                    w = _np.sort(_np.r_[w, w_loop])
                loop_freq_idx = _np.where(w == w_loop)[0][0]
            _zl = imp.get_impedance(w=w)
            # low-level control loop makes the impedance
            # go to zero at the loop actuation frequency
            # for active impedance sources
            _zl[loop_freq_idx] = 0 + 0j
            total_zl += _zl
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
        zl_wp = self.get_impedance(w=w, apply_filter=True)
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
        h = self.ring.harm_num
        w0 = self.ring.rev_ang_freq

        if dist is None:
            dist = self.distributions
        fillpattern = self.fillpattern
        zgrid = self.zgrid

        zn_ph = (2j*_PI/h) * _np.arange(h)[None, :]
        z_ph = (1j*w0/_LSPEED) * zgrid[None, :]

        ps, zl_wps = self.get_harmonics_impedance_and_filling()
        ps = ps[:, None]
        expphz = _ne.evaluate('exp(ps*z_ph)')
        expphn = _ne.evaluate('exp(ps*zn_ph)')
        ps = ps.ravel()

        # remove last point to do not overlap domains
        dist_beam = (fillpattern[:, None]*dist[:, :-1]).ravel()
        dist_dft = _np.fft.rfft(dist_beam) * (zgrid[1] - zgrid[0])
        dist_dft = dist_dft[ps]
        # shift phase due to difference of DFT reference frame
        # (initial point z=0) and our reference (initial point z=-lambda_rf/2)
        dist_dft *= _np.exp(1j*ps*_PI/h)
        dist_dft *= zl_wps.conj()
        # sum over positive frequencies only -> factor 2
        expphn *= dist_dft[:, None]
        harm_volt = -2*_np.dot(expphn.T, expphz)
        return harm_volt.real

    def calc_voltage_harmonic_cavity_wake(
            self, wake_source, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
        fillpattern = self.fillpattern[:, None]
        zgrid = self.zgrid

        h = self.ring.harm_num
        circum = self.ring.circum
        rev_time = self.ring.rev_time

        alpha = wake_source.alpha
        beta = wake_source.beta
        wrbar = wake_source.ang_freq_bar
        rsh = wake_source.shunt_impedance

        if self._exp_z is None:
            self._exp_z = _ne.evaluate('exp(beta*zgrid)')[None, :]

        dist_exp_z = _np.zeros(dist.shape, dtype=_np.complex)
        dist_exp_z += dist
        dist_exp_z *= fillpattern
        dist_exp_z *= self._exp_z
        dz = zgrid[1] - zgrid[0]
        dist_fourier = _mytrapz(dist_exp_z, dz)

        # NOTE: Alternative implementation without matrix multiplication. This
        # calculation did not reduce the evaluation time too much, then the
        # original implementation was kept for readability.
        # ind = _np.arange(h)
        # exp_betac0 = _np.exp(beta*circum)
        # exp_ind = _ne.evaluate('exp(beta*circum*ind/h)')
        # vec = exp_ind * dist_fourier
        # cum_sum = _np.r_[0, _np.cumsum(vec)]

        # V = exp_betac0*cum_sum[:-1]
        # V += cum_sum[-1]
        # V -= cum_sum[:-1]
        # V /= exp_ind
        # V /= exp_betac0 - 1

        if self._wake_matrix is None:
            exp_betac0 = _np.exp(-beta*circum)
            Ll = 1/(1-exp_betac0)  # buckets ahead current one (l<n)
            Gl = Ll*exp_betac0     # buckets behind current one (l>=n)
            wmat = Ll*_np.tri(h, h, -1,)
            wmat += Gl*_np.tri(h, h).T

            ind = _np.arange(h)
            diff = ind[:, None] - ind[None, :]
            wmat *= _ne.evaluate('exp(-beta*circum*diff/h)')
            self._wake_matrix = wmat
        V = _np.dot(self._wake_matrix, dist_fourier)

        Vt = _mytrapz(dist_exp_z, dz, cumul=True)
        Vt += V[:, None]
        Vt /= self._exp_z

        harm_volt = Vt.real
        harm_volt -= alpha / wrbar * Vt.imag
        harm_volt *= -2*alpha * rsh * rev_time
        return harm_volt

    def calc_longitudinal_equilibrium(
            self, niter=100, tol=1e-10, beta=1, m=3, print_flag=True):
        """."""
        self.print_flag = print_flag
        dists = [self.distributions, ]
        dists = self._apply_anderson_acceleration(
            dists, niter, tol, beta=beta, m=m)
        dists = [self._reshape_dist(rho) for rho in dists]
        self.distributions = dists[-1]
        # Flush pre-calculated data
        self._wake_matrix = None
        self._exp_z = None
        return dists

    def get_generator_voltage(self, feedback=True, method='phasor'):
        """."""
        if feedback:
            if self.beamload_active is None:
                raise ValueError(
                    'Feedback is on but beam loading voltage is None!')
            if method == 'phasor':
                # Phasor compensation
                _vg = self._feedback_phasor()
            elif method == 'lstqr':
                # Least-squares minimization
                _vg = self._feedback_least_squares()
            else:
                raise ValueError(
                    "Wrong feedback method: must be 'phasor' or 'lstqr'")
        else:
            _vg = self.ring.get_voltage_waveform(self.zgrid)[None, :]
        return _vg

    # -------------------- instabilities calculations -------------------------
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
            Zlp = self.get_impedance(w=wp, apply_filter=False)
            Zln = self.get_impedance(w=wn, apply_filter=False)
            growth = const*(wp*Zlp.real - wn*Zln.real)
        return growth

    def calc_tuneshifts_cbi(
            self, w, m=1, nbun_fill=None, radiation=False):
        """."""
        ring = self.ring
        num_bun = ring.num_bun
        dampte = ring.dampte

        ring.num_bun = nbun_fill if nbun_fill is not None else ring.harm_num
        if not radiation:
            ring.dampte = _np.inf

        if _np.array(w).size == 2:
            # Automatically sample the impedance at revolution harmonics,
            # only the min and max frequencies are required: w = [w_min, w_max]
            Zl = self.get_impedance
        else:
            Zl = self.get_impedance(w=w, apply_filter=False)

        deltaw, wp, interpol_Z, spectrum = ring.longitudinal_cbi(
            w=w, Zl=Zl, m=m, inverse=False, full=True)

        # Relative tune-shifts must be multiplied by ws
        deltaw *= (ring.sync_tune * ring.rev_ang_freq)

        ring.num_bun = num_bun
        ring.dampte = dampte
        return deltaw, Zl, wp, interpol_Z, spectrum

    def calc_mode_coupling(
            self, w, cbmode, max_azi=10, max_rad=12, nbun_fill=None,
            modecoup_matrix=None, fokker_matrix=None):
        """."""
        ring = self.ring
        num_bun = ring.num_bun
        dampte = ring.dampte

        ring.num_bun = nbun_fill if nbun_fill is not None else ring.harm_num

        if _np.array(w).size == 2:
            Zl = self.get_impedance
        else:
            Zl = self.get_impedance(w=w, apply_filter=False)

        eigenfreq, modecoup_matrix, fokker_matrix = \
            ring.longitudinal_mode_coupling(
                w=w, Zl=Zl, cbmode=cbmode, max_azi=max_azi, max_rad=max_rad,
                modecoup_matrix=modecoup_matrix, fokker_matrix=fokker_matrix)

        # Relative tune-shifts must be multiplied by ws
        eigenfreq *= (ring.sync_tune * ring.rev_ang_freq)

        ring.num_bun = num_bun
        ring.dampte = dampte
        return eigenfreq, modecoup_matrix, fokker_matrix

    # -------------------- auxiliary methods ----------------------------------
    def _apply_anderson_acceleration(self, dists, niter, tol, m=3, beta=1):
        """."""
        if beta < 0:
            raise Exception('relaxation parameter beta must be positive.')

        xold = dists[-1].ravel()
        xnew = self._ffunc(xold)
        dists.append(xnew)

        where = 0
        G_k = _np.zeros((xnew.size, m), dtype=float)
        X_k = _np.zeros((xnew.size, m), dtype=float)
        mat = _np.zeros((xnew.size, m), dtype=float)

        gold = xnew - xold
        gnew = self._ffunc(xnew) - xnew
        G_k[:, where] = gnew - gold
        X_k[:, where] = gold
        mat[:, where] = gnew
        where += 1
        where %= m

        for k in range(niter):
            t0 = _time.time()
            gamma_k = _np.linalg.lstsq(G_k, gnew, rcond=None)[0]
            xold = xnew
            xnew = xold + gnew
            xnew -= mat @ gamma_k
            xnew *= beta
            if beta != 1:
                xnew += (1-beta) * (xold - X_k @ gamma_k)
            dists.append(xnew)

            gold = gnew
            gnew = self._ffunc(xnew) - xnew
            G_k[:, where] = gnew - gold
            X_k[:, where] = xnew - xold
            mat[:, where] = G_k[:, where] + X_k[:, where]
            where += 1
            where %= m

            diff = self._reshape_dist(gnew)
            dz = self.zgrid[1] - self.zgrid[0]
            diff = _mytrapz(_np.abs(diff), dz)
            idx = _np.argmax(diff)
            tf = _time.time() - t0
            if self.print_flag:
                print(
                    f'Iter.: {k+1:03d}, Dist. Diff.: {diff[idx]:.3e}' +
                    f' (bucket {idx:03d}), E.T.: {tf:.3f}s')
            if diff[idx] < tol:
                if self.print_flag:
                    print('distribution ok!')
                break
        return dists

    def _create_freqs(self):
        w0 = self.ring.rev_ang_freq
        p = _np.arange(0, self.max_mode)
        return p*w0

    def _ffunc(self, xk, feedback=True):
        """Haissinski operator."""
        xk = self._reshape_dist(xk)
        total_volt = _np.zeros(xk.shape)
        self.beamload_active = _np.zeros(xk.shape)
        idx_wake = self._get_wake_types_idx()
        if idx_wake:
            wake_sources = [self.impedance_sources[idx] for idx in idx_wake]
            # Wakes need to be evaluated one by one
            _func = self.calc_voltage_harmonic_cavity_wake
            for wake in wake_sources:
                if wake.active_passive == ImpedanceSource.ActivePassive.Active:
                    self.beamload_active = _func(wake_source=wake, dist=xk)
                    self._wake_matrix = None
                    self._exp_z = None
                total_volt += _func(wake_source=wake, dist=xk)
                self._wake_matrix = None
                self._exp_z = None
        else:
            # Impedances can be summed and calculated once
            _func = self.calc_voltage_harmonic_cavity_impedance
            total_volt += _func(dist=xk)

        total_volt += self.get_generator_voltage(
            feedback=feedback)
        fxk, _ = self.calc_distributions_from_voltage(total_volt)
        return fxk.ravel()

    def _get_impedance_types_idx(self):
        """."""
        imp_idx = []
        for idx, imp in enumerate(self.impedance_sources):
            if imp.calc_method == ImpedanceSource.Methods.Impedance:
                imp_idx.append(idx)
        return imp_idx

    def _get_wake_types_idx(self):
        """."""
        wake_idx = []
        for idx, imp in enumerate(self.impedance_sources):
            if imp.calc_method == ImpedanceSource.Methods.Wake:
                wake_idx.append(idx)
        return wake_idx

    def _feedback_phasor(self):
        vgap = self.ring.gap_voltage
        wrf = 2*_PI*self.ring.rf_freq
        phase = wrf*self.zgrid/_LSPEED
        dz = self.zgrid[1] - self.zgrid[0]

        vref_phasor = vgap*_np.exp(
            1j*(_PI/2 - self.ring.sync_phase))
        vbeamload_phasor = _np.mean(_mytrapz(
            self.beamload_active*_np.exp(1j*phase)[None, :], dz))
        vbeamload_phasor *= 2/(self.zgrid[-1]-self.zgrid[0])
        vg_phasor = vref_phasor - vbeamload_phasor
        vg = _np.real(vg_phasor*_np.exp(-1j*phase))
        return vg[None, :]

    def _feedback_least_squares(self):
        x0 = [self.ring.gap_voltage, 0]
        wrf = 2*_PI*self.ring.rf_freq
        phase = wrf*self.zgrid/_LSPEED
        dz = self.zgrid[1] - self.zgrid[0]

        vref = self.ring.get_voltage_waveform(self.zgrid)
        res = _least_squares(
            fun=self._feedback_err, x0=x0,
            args=(phase, dz, self.beamload_active, vref),
            method='lm')
        vg = LongitudinalEquilibrium._generator_model(
            phase, res.x[0], res.x[1])
        return vg[None, :]

    @staticmethod
    def _feedback_err(x, *args):
        phase, dz, vbeamload, vref = args
        vgen = LongitudinalEquilibrium._generator_model(phase, x[0], x[1])
        err = (vgen[None, :] + vbeamload) - vref[None, :]
        err = _mytrapz(err*err, dz)
        return err

    @staticmethod
    def _generator_model(phase, a, b):
        return a*_np.sin(phase) + b*_np.cos(phase)

    def _reshape_dist(self, dist):
        return dist.reshape((self.ring.harm_num, self.zgrid.size))
