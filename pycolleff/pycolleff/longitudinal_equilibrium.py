"""."""

import time as _time

import numexpr as _ne
import numpy as _np
from mathphys.constants import light_speed as _c
from mathphys.functions import get_namedtuple as _get_namedtuple
from scipy.fft import fft as _fft, irfft as _irfft, rfft as _rfft
from scipy.integrate import quad as _quad
from scipy.interpolate import interp1d as _interp1d
from scipy.optimize import least_squares as _least_squares, root as _root
from scipy.special import gamma as _gammafunc

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
        intg[..., 1:] = _ne.evaluate("(y1 + y2)*dx/2.0")
        res = _np.cumsum(intg, axis=-1)
        res.shape = intg.shape
        return res
    else:
        return _ne.evaluate("(y1 + y2)*dx/2.0").sum(axis=-1)


class ImpedanceSource:
    """."""

    Methods = _get_namedtuple(
        "Methods",
        ["ImpedanceDFT", "ImpedanceModeSel", "Wake", "UniformFillAnalytic"],
    )
    ActivePassive = _get_namedtuple("ActivePassive", ["Active", "Passive"])

    def __init__(
        self,
        Rs=0,
        Q=0,
        ang_freq=0,
        harm_rf=3,
        calc_method=Methods.ImpedanceDFT,
        active_passive=ActivePassive.Passive,
    ):
        """."""
        self._calc_method = None
        self._active_passive = None

        self.ang_freq = ang_freq
        self.Q = Q
        self.shunt_impedance = Rs

        self.harm_rf = harm_rf
        self.ang_freq_rf = 0
        self._loop_ctrl_freq = 0
        self._loop_ctrl_transfer = 0
        self._zl_table = None
        self._ang_freq_table = None
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
            raise ValueError("Wrong value for calc_method.")

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
            raise ValueError("Wrong value for active_passive.")

    def get_impedance(self, w):
        """."""
        if self.zl_table is None:
            imp = _imp.longitudinal_resonator(
                Rs=self.shunt_impedance, Q=self.Q, wr=self.ang_freq, w=w
            )
        else:
            w_tab = self.ang_freq_table
            zl_tab = self.zl_table
            imp = _np.interp(w, w_tab, zl_tab.imag) * 1j
            imp += _np.interp(w, w_tab, zl_tab.real)
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
        return 2 * _PI * self._loop_ctrl_freq

    @loop_ctrl_ang_freq.setter
    def loop_ctrl_ang_freq(self, value):
        """."""
        self._loop_ctrl_freq = value / 2 / _PI

    @property
    def loop_ctrl_transfer(self):
        """."""
        return self._loop_ctrl_transfer

    @loop_ctrl_transfer.setter
    def loop_ctrl_transfer(self, func):
        """."""
        self._loop_ctrl_transfer = func

    @property
    def RoverQ(self):
        """."""
        return self.shunt_impedance / self.Q

    @property
    def detune_w(self):
        """."""
        return self.ang_freq - self.harm_rf * self.ang_freq_rf

    @detune_w.setter
    def detune_w(self, value):
        """."""
        wr = self.harm_rf * self.ang_freq_rf + value
        self.ang_freq = wr

    @property
    def detune_freq(self):
        """."""
        return self.detune_w / 2 / _PI

    @property
    def alpha(self):
        """."""
        return self.ang_freq / 2 / self.Q

    @property
    def ang_freq_bar(self):
        """."""
        wr_ = self.ang_freq
        alpha = self.alpha
        return (wr_ * wr_ - alpha * alpha) ** 0.5

    @property
    def beta(self):
        """."""
        return (self.alpha - 1j * self.ang_freq_bar) / _c

    @property
    def detune_angle(self):
        """."""
        Q = self.Q
        nharm = self.harm_rf
        wrf = self.ang_freq_rf
        wr = self.ang_freq
        if wr == 0:
            raise Exception("wr cannot be zero!")
        if wrf == 0:
            raise Exception("wrf cannot be zero!")
        tan = Q * (wr / (nharm * wrf) - nharm * wrf / wr)
        angle = _np.arctan2(tan, 1)
        return angle

    @detune_angle.setter
    def detune_angle(self, value):
        Q = self.Q
        nharm = self.harm_rf
        wrf = self.ang_freq_rf

        delta = _np.tan(value) / 2 / Q
        self.ang_freq = nharm * wrf * (delta + (1 + delta**2) ** (1 / 2))

    @property
    def zl_table(self):
        """."""
        return self._zl_table

    @zl_table.setter
    def zl_table(self, value):
        self._zl_table = value

    @property
    def ang_freq_table(self):
        """."""
        return self._ang_freq_table

    @ang_freq_table.setter
    def ang_freq_table(self, value):
        self._ang_freq_table = value

    def to_dict(self):
        """Save state to dictionary."""
        return dict(
            ang_freq=self.ang_freq,
            Q=self.Q,
            shunt_impedance=self.shunt_impedance,
            harm_rf=self.harm_rf,
            ang_freq_rf=self.ang_freq_rf,
        )

    def from_dict(self, dic):
        """Load state from dictionary."""
        self.ang_freq = dic.get("ang_freq", self.ang_freq)
        self.Q = dic.get("Q", self.Q)
        self.shunt_impedance = dic.get("shunt_impedance", self.shunt_impedance)
        self.harm_rf = dic.get("harm_rf", self.harm_rf)
        self.ang_freq_rf = dic.get("ang_freq_rf", self.ang_freq_rf)

    def __str__(self):
        """."""
        stmp = "{0:20s}: {1:}  {2:s}\n".format
        ftmp = "{0:20s}: {1:3.2f}  {2:s}\n".format
        etmp = "{0:20s}: {1:.2e}  {2:s}\n".format
        mega = 1e-6
        stg = stmp("calc_method", self.calc_method_str, "")
        stg += stmp("active_passive", self.active_passive_str, "")
        stg += ftmp("ang_freq_rf", self.ang_freq_rf * mega, "[Mrad/s]")
        stg += ftmp("ang_freq", self.ang_freq * mega, "[Mrad/s]")
        stg += ftmp("shunt_impedance", self.shunt_impedance * mega, "[MOhm]")
        stg += etmp("Q", self.Q, "")
        stg += ftmp("RoverQ", self.RoverQ, "[Ohm]")
        stg += ftmp("harm_rf", self.harm_rf, "")
        stg += ftmp("detune_angle", self.detune_angle, "[rad]")
        stg += ftmp("detune_freq", self.detune_freq / 1e3, "[kHz]")
        stg += ftmp("detune_w", self.detune_w, "[rad/s]")
        stg += ftmp("alpha", self.alpha, "[rad/s]")
        stg += ftmp("ang_freq_bar", self.ang_freq_bar * mega, "[Mrad/s]")
        return stg


class LongitudinalEquilibrium:
    """."""

    FeedbackMethod = _get_namedtuple(
        "FeedbackMethod", ["Phasor", "LeastSquares"]
    )

    def __init__(self, ring: _Ring, impedance_sources: list, fillpattern=None):
        """."""
        self._zgrid = None
        self._dist = None
        self._fillpattern = None
        self._main_voltage = None
        self._calc_fun = None
        self._calc_method = None
        self._print_flag = False
        self._wake_matrix = None
        self.beamload_active = None
        self.total_voltage = None

        self.ring = ring
        self.impedance_sources = impedance_sources
        self.fillpattern = fillpattern
        self._max_mode = 10 * self.ring.harm_num
        self.min_mode0_ratio = 1e-9

        self.feedback_method = self.FeedbackMethod.LeastSquares
        self.feedback_on = False

        self.main_ref_phase_offset = 0.0  # [radian]

        # If feedback_on = True, this will not be used
        self.main_ref_amp = self.ring.gap_voltage
        self.main_ref_phase = self.ring.sync_phase

        # If feedback_on = True, this will be updated
        self.main_gen_amp_mon = None
        self.main_gen_phase_mon = None

    @property
    def feedback_method_str(self):
        """."""
        return self.FeedbackMethod._fields[self._feedback_method]

    @property
    def feedback_method(self):
        """."""
        return self._feedback_method

    @feedback_method.setter
    def feedback_method(self, value):
        if value is None:
            return
        if isinstance(value, str):
            self._feedback_method = self.FeedbackMethod._fields.index(value)
        elif int(value) in self.FeedbackMethod:
            self._feedback_method = int(value)
        else:
            raise ValueError("Wrong value for feedback_method.")

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
            raise ValueError("Wrong shape for voltage.")
        self._main_voltage = value
        self.distributions, _ = self.calc_distributions_from_voltage(
            self._main_voltage
        )

    @property
    def fillpattern(self):
        """."""
        return self._fillpattern

    @fillpattern.setter
    def fillpattern(self, value):
        if value.size != self.ring.harm_num:
            raise ValueError("Wrong size for fillparttern.")
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
            raise ValueError("Distributions must have 2 dimensions.")
        elif value.shape[0] != self.ring.harm_num:
            raise ValueError("First dimension must be equal ring.harm_num.")
        elif value.shape[1] != self._zgrid.size:
            raise ValueError("Second dimension must be equal zgrid.size.")
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
                imp.to_dict() for imp in self.impedance_sources
            ],
            zgrid=self._zgrid,
            dist=self._dist,
            fillpatern=self._fillpattern,
            main_voltage=self._main_voltage,
            max_mode=self.max_mode,
            min_mode0_ratio=self.min_mode0_ratio,
            calc_method=self.calc_method_str,
        )

    def from_dict(self, dic):
        """Load state from dictionary."""
        self.ring.from_dict(dic.get("ring", dict()))
        imps = []
        for imp in dic.get("impedance_sources", self.impedance_sources):
            _imp = ImpedanceSource()
            _imp.from_dict(imp)
            imps.append(_imp)
        self.impedance_sources = imps
        self._zgrid = dic.get("zgrid", self._zgrid)
        self._dist = dic.get("dist", self._dist)
        self._fillpattern = dic.get("fillpattern", self._fillpattern)
        self._main_voltage = dic.get("main_voltage", self._main_voltage)
        self.max_mode = dic.get("max_mode", self.max_mode)
        self.min_mode0_ratio = dic.get("min_mode0_ratio", self.min_mode0_ratio)
        self.calc_method = dic.get("calc_method", self.calc_method)

    def create_zgrid(self, nr_points=1001, sigmas=30):
        """."""
        return sigmas * self.ring.bunlen * _np.linspace(-1, 1, nr_points)

    @staticmethod
    def calc_moments(zgrid, dist):
        """."""
        dz = zgrid[1] - zgrid[0]
        zm = _mytrapz(zgrid[None, :] * dist, dz)

        zgrid2 = zgrid * zgrid
        z2 = _mytrapz(zgrid2[None, :] * dist, dz)
        return zm, _np.sqrt(z2 - zm**2)

    def get_gaussian_distributions(self, sigmaz, z0=0):
        """."""
        dz = self.zgrid[1] - self.zgrid[0]

        arg = (self.zgrid - z0) / sigmaz
        dist = _np.exp(-(arg**2) / 2)
        dist /= _mytrapz(dist, dz)
        dist = _np.tile(dist, (self.ring.harm_num, 1))
        return dist

    def calc_harmonic_voltage_for_flat_potential(self, harm_rf=3):
        """."""
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        n2 = harm_rf**2
        kharm = 1 / n2 - ((U0 / Vrf) ** 2) / (n2 - 1)
        return kharm ** (1 / 2)

    def calc_detune_for_fixed_harmonic_voltage(
        self, peak_harm_volt, harm_rf=3, Rs=0, form_factor=None
    ):
        """."""
        I0 = _np.sum(self.fillpattern)
        # TODO: This way of including the form factor is temporary. Fix it.
        wr = 2 * _PI * self.ring.rf_freq * harm_rf
        if form_factor is None:
            form_factor = self.calc_fourier_transform(wr)[self.filled_buckets]
        ib = 2 * I0 * _np.abs(form_factor).mean()
        arg = peak_harm_volt / ib / Rs
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
        pot = -_mytrapz(voltage - U0, dz, cumul=True)

        # subtract minimum value for all bunches
        pot -= _np.min(pot, axis=1)[:, None]

        const = self.ring.espread**2
        const *= self.ring.mom_comp
        const *= self.ring.circum
        # normalize by E0
        const *= self.ring.energy

        dist = _ne.evaluate("exp(-pot/const)")
        # distribution must be normalized
        dist /= _mytrapz(dist, dz)[:, None]
        if flag:
            dist = _np.tile(dist, (self.ring.harm_num, 1))
        return dist, pot / const

    def calc_fourier_transform(self, w, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
        arg = _np.exp((1j * w / _c) * self.zgrid)[None, :]
        arg = _ne.evaluate("dist * arg")
        dz = self.zgrid[1] - self.zgrid[0]
        return _mytrapz(arg, dz)

    def get_impedance(self, w=None, apply_filter=False):
        """."""
        if w is None:
            w = self._create_freqs()
        total_zl = _np.zeros(w.shape, dtype=complex)
        imp_idx = self._get_impedance_types_idx()
        for iidx in imp_idx:
            imp = self.impedance_sources[iidx]
            cond = imp.active_passive == ImpedanceSource.ActivePassive.Active
            cond &= apply_filter
            _zl0 = imp.get_impedance(w=w)
            if cond:
                # closed-loop impedance
                transf = imp.loop_ctrl_transfer(w, imp.loop_ctrl_ang_freq)
                _zl = _zl0 / (1 + transf * _zl0)
            else:
                # open-loop impedance
                _zl = _zl0
            total_zl += _zl
        return total_zl

    def get_harmonics_impedance_and_filling(self, w=None):
        """."""
        if w is None:
            w = self._create_freqs()
        h = self.ring.harm_num
        zl_wp = self.get_impedance(w=w, apply_filter=True)
        fill_fft = _fft(self.fillpattern)
        fill_fft = _np.tile(fill_fft, (zl_wp.size // h, 1)).ravel()
        zl_fill = _np.abs(zl_wp * fill_fft)

        # # select modes based on max peak neighbors
        # peak = _np.argmax(zl_fill)
        # nr_modes = 0
        # if self.max_mode is not None:
        #     nr_modes = (self.max_mode // 2)
        # modes = _np.arange(nr_modes + 1)
        # modes = _np.r_[-modes[:0:-1], modes] + peak
        # out = modes, zl_wp[modes], zl_fill

        # select modes based sorted imp * fill spectrum
        modes = _np.where(zl_fill >= zl_fill.max() * self.min_mode0_ratio)[0]

        idx_sort = _np.argsort(_np.abs(zl_fill[modes]))[::-1]
        if self.max_mode is not None:
            idx_sort = idx_sort[: self.max_mode]
        out = modes[idx_sort], zl_wp[modes][idx_sort], zl_fill
        return out

    def calc_induced_voltage_uniform_filling(self, wake_source, dist=None):
        """."""
        if dist is None:
            dist = self.distributions
        wr = wake_source.harm_rf * wake_source.ang_freq_rf
        form = self.calc_fourier_transform(wr, dist=dist)
        F0 = _np.abs(form)[0]
        Phi0 = _np.angle(form)[0]

        It = _np.sum(self.fillpattern)
        ang = wake_source.detune_angle
        Rs = wake_source.shunt_impedance

        volt = -2 * It * F0 * Rs * _np.cos(ang)
        volt *= _np.cos(wr * self.zgrid / _c + ang - Phi0)
        return _np.tile(volt, (self.ring.harm_num, 1))

    def calc_induced_voltage_impedance_mode_selection(self, dist=None):
        """."""
        h = self.ring.harm_num
        w0 = self.ring.rev_ang_freq

        if dist is None:
            dist = self.distributions
        fillpattern = self.fillpattern
        zgrid = self.zgrid

        zn_ph = (2j * _PI / h) * _np.arange(h)[None, :]
        z_ph = (1j * w0 / _c) * zgrid[None, :]

        ps, zl_wps, _ = self.get_harmonics_impedance_and_filling()
        ps = ps[:, None]
        zl_wp = _ne.evaluate("exp(ps*z_ph)")
        zl_wp *= zl_wps[:, None].conj()

        expph = _ne.evaluate("exp(-ps*zn_ph)")
        harm_volt = _np.zeros((h, zgrid.size), dtype=complex)
        for idx, p in enumerate(ps):
            dist_fourier = self.calc_fourier_transform(w=p * w0, dist=dist)

            exp_phase = expph[idx]
            beam_part = _np.einsum(
                "i,i,i", exp_phase, fillpattern, dist_fourier.conj()
            )
            beam_part = beam_part / exp_phase

            # sum over positive frequencies only -> factor 2
            harm_volt += -2 * zl_wp[idx] * beam_part[:, None]
        return harm_volt.real

    def calc_induced_voltage_impedance_dft(self, dist=None):
        """."""
        if dist is None:
            dist = self.distributions

        ps, zl_wps, _ = self.get_harmonics_impedance_and_filling()

        did_zero_pad = False
        rf_lamb = self.ring.rf_lamb
        if self.zgrid[0] != -rf_lamb / 2 or self.zgrid[-1] != rf_lamb / 2:
            dist, idx_ini = self._do_zero_padding(dist)
            did_zero_pad = True

        # remove last point to do not overlap domains
        dist_beam = (self.fillpattern[:, None] * dist[:, :-1]).ravel()
        dist_dft_ = _rfft(dist_beam)

        # calculate with DFT
        dist_dft = _np.zeros(dist_dft_.size, dtype=complex)
        dist_dft[ps] = dist_dft_[ps]
        dist_dft[ps] *= zl_wps.conj()
        _harm_volt = (-self.ring.circum) * _irfft(dist_dft)
        harm_volt = _np.zeros_like(dist, dtype=complex)
        harm_volt[:, :-1] = _harm_volt.reshape((dist.shape[0], -1))
        harm_volt[:-1, -1] = harm_volt[1:, 0]
        harm_volt[-1, -1] = harm_volt[0, 0]
        if did_zero_pad:
            harm_volt = harm_volt[:, idx_ini : idx_ini + self.zgrid.size]
        return harm_volt.real

    def calc_induced_voltage_wake(self, wake_source, dist=None):
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
            self._exp_z = _ne.evaluate("exp(beta*zgrid)")[None, :]

        dist_exp_z = _np.zeros(dist.shape, dtype=complex)
        dist_exp_z += dist
        dist_exp_z *= fillpattern
        dist_exp_z *= self._exp_z
        dz = zgrid[1] - zgrid[0]
        Sn = _mytrapz(dist_exp_z, dz, cumul=True)
        dist_laplace = Sn[:, -1]

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
            exp_betac0 = _np.exp(-beta * circum)
            # buckets ahead current one (l<n)
            log_Ll = -_np.log(1 - exp_betac0)
            # buckets behind current one (l>=n)
            log_Gl = log_Ll - beta * circum
            log_wmat = log_Ll * _np.tri(h, h, -1)
            log_wmat += log_Gl * _np.tri(h, h).T
            ind = _np.arange(h)
            diff = ind[:, None] - ind[None, :]
            log_wmat += -beta * circum * diff / h
            self._wake_matrix = _ne.evaluate("exp(log_wmat)")
        V = _np.dot(self._wake_matrix, dist_laplace)
        Vt = (Sn + V[:, None]) / self._exp_z

        harm_volt = Vt.real
        harm_volt -= alpha / wrbar * Vt.imag
        harm_volt *= -2 * alpha * rsh * rev_time
        return harm_volt

    def calc_longitudinal_equilibrium(
        self, niter=100, tol=1e-10, beta=1, m=3, print_flag=True
    ):
        """."""
        self.print_flag = print_flag
        dists = [
            self.distributions,
        ]
        dists, converged = self._apply_anderson_acceleration(
            dists, niter, tol, beta=beta, m=m
        )

        # dists = self._apply_random_convergence(dists, niter, tol)
        dists = [self._reshape_dist(rho) for rho in dists]
        self.distributions = dists[-1]
        # Flush pre-calculated data
        self._wake_matrix = None
        self._exp_z = None
        return dists, converged

    def _apply_random_convergence(self, dists, niter, tol):
        xold = dists[-1].ravel()
        for k in range(niter):
            xnew = self._ffunc(xold)
            dists.append(xnew)
            diff = self._reshape_dist(xnew - xold)
            dz = self.zgrid[1] - self.zgrid[0]
            diff = _mytrapz(_np.abs(diff), dz)
            idx = _np.argmax(diff)
            if self.print_flag:
                print(
                    f"Iter.: {k+1:03d}, Dist. Diff.: {diff[idx]:.3e}"
                    + f" (bucket {idx:03d})"
                )
                print("-" * 20)
            if diff[idx] < tol:
                if self.print_flag:
                    print("distribution ok!")
                break
            r = _np.random.randn() / 2
            xold = (1 - r) * xnew + r * xold
        return dists

    def get_generator_voltage(self):
        """."""
        if self.feedback_on:
            err = "Feedback is on but there is no active beam loading voltage!"
            if self.beamload_active is not None:
                val = _np.sum(self.beamload_active)
                if not val:
                    raise ValueError(err)
            else:
                raise ValueError(err)
            if self.feedback_method == self.FeedbackMethod.Phasor:
                # Phasor compensation
                _vg = self._feedback_phasor()
            elif self.feedback_method == self.FeedbackMethod.LeastSquares:
                # Least-squares minimization
                _vg = self._feedback_least_squares()
            else:
                raise ValueError(
                    "Wrong feedback method: must be"
                    + "'Phasor' or 'LeastSquares'"
                )
        else:
            amp = self.main_ref_amp
            phase = self.main_ref_phase
            phase += self.main_ref_phase_offset
            _vg = self.ring.get_voltage_waveform(
                self.zgrid, amplitude=amp, phase=phase
            )[None, :]
        return _vg

    def calc_synchrotron_frequency(
        self, total_voltage, method="action", max_amp=5, nrpts=100
    ):
        """Calculate synchrotron frequencies for given total voltage."""
        # _warnings.filterwarnings("error")

        lambda0, phiz = self.calc_distributions_from_voltage(total_voltage)
        zgrid = self.zgrid
        ring = self.ring
        phiz = phiz[0, :]

        z0, sigmaz0 = self.calc_moments(zgrid, lambda0)
        z0, sigmaz0 = z0[0], sigmaz0[0]

        zmin = zgrid[_np.argmin(phiz)]
        zgrid -= zmin

        alpha = ring.mom_comp
        sigmae0 = ring.espread
        factor = alpha * _c * sigmae0**2

        if method == "action":
            phiz_interp = _interp1d(zgrid, phiz * factor, kind="cubic")

            dz = max_amp * sigmaz0 / nrpts
            # avoid zero amplitude
            zamps = _np.arange(1, nrpts) * dz

            actions, periods, hamiltonian, zamps_ = [], [], [], []

            def energy_deviation(z):
                return hamiltonian0i - phiz_interp(z)

            def intg(z):
                return _np.sqrt(
                    (2 / alpha / _c) * _np.abs(hamiltonian0i - phiz_interp(z))
                )

            def iintg(z):
                return 1 / intg(z)

            for zamp in zamps:
                zri = z0 * 0 + zamp
                zli = z0 * 0 - zamp
                hamiltonian0i = phiz_interp(zri)

                turn_pts = _root(energy_deviation, x0=zli)
                if turn_pts.success:
                    zli = turn_pts.x[0]
                else:
                    raise Exception("Problem in finding turning points.")
                zli, zri = (zli, zri) if zli <= zri else (zri, zli)

                action, _ = _quad(intg, zli, zri)
                period, _ = _quad(iintg, zli, zri)

                actions.append(action / _PI)
                periods.append(period * 2 / alpha / _c)
                hamiltonian.append(hamiltonian0i)
                zamps_.append(zamp)

            actions, periods, hamiltonian, zamps = (
                _np.array(actions),
                _np.array(periods),
                _np.array(hamiltonian),
                _np.array(zamps_),
            )

            freqs_deriv = _np.gradient(hamiltonian, actions) / 2 / _PI
            freqs = 1 / periods
            nan_idx = ~(
                _np.isnan(actions) | _np.isnan(freqs) | _np.isnan(freqs_deriv)
            )
            diverge_idx1 = (_np.abs(freqs) < ring.rev_freq) & (freqs >= 0)
            diverge_idx2 = (_np.abs(freqs_deriv) < ring.rev_freq) & (
                freqs_deriv >= 0
            )
            filter_idx = nan_idx & diverge_idx1 & diverge_idx2
            actions, freqs, freqs_deriv, hamiltonian, zamps = (
                actions[filter_idx],
                freqs[filter_idx],
                freqs_deriv[filter_idx],
                hamiltonian[filter_idx],
                zamps[filter_idx],
            )

            lambda0_ = _np.exp(-hamiltonian / factor)
            lambda0_ /= _np.trapz(lambda0_, actions)

            fs_avg = _np.trapz(freqs * lambda0_, actions)
            fs_std = _np.trapz(freqs * freqs * lambda0_, actions)
            fs_std = _np.sqrt(fs_std - fs_avg**2)

            out = dict()
            out["sync_freq"] = freqs
            out["sync_freq_numeric_derivative"] = freqs_deriv
            out["avg_sync_freq"] = fs_avg
            out["std_sync_freq"] = fs_std
            out["action_distribution"] = lambda0_
            out["action"] = actions
            out["hamiltonian"] = hamiltonian
            out["amplitude"] = zamps
            return out

        elif method == "derivative":
            wrf_c = ring.rf_ang_freq / _c
            factor = _np.sqrt(
                alpha * ring.harm_num / (2 * _PI * ring.energy) / wrf_c
            )

            fil = _np.abs(zgrid) < max_amp * sigmaz0
            zgrid_ = zgrid[fil]

            dv = -_np.gradient(total_voltage[0, fil], zgrid_)
            remove_neg = dv > 0
            zgrid_, lambda0_, dv = (
                zgrid_[remove_neg],
                lambda0[0, fil][remove_neg],
                dv[remove_neg],
            )

            lambda0_ /= _np.trapz(lambda0_, zgrid_)
            freqs = factor * _np.sqrt(dv) * ring.rev_freq
            fs_avg = _np.trapz(freqs * lambda0_, zgrid_)
            fs_std = _np.sqrt(
                _np.trapz((freqs - fs_avg) ** 2 * lambda0_, zgrid_)
            )

            out = dict()
            out["sync_freq"] = freqs
            out["avg_sync_freq"] = fs_avg
            out["std_sync_freq"] = fs_std
            return out

    def calc_canonical_transformation(self, res, total_voltage):
        """."""
        npts = 100
        twopi = 2 * _np.pi
        U0 = self.ring.en_lost_rad
        E0 = self.ring.energy
        T0 = self.ring.rev_time
        vtotal = (total_voltage - U0) / E0 / T0
        zgrid = self.zgrid

        zj = []
        phi = []
        for _, fs in enumerate(res["sync_freq"]):
            t = _np.linspace(0, 1, npts) * (1 / fs)
            dt = 1 / fs / npts
            z0 = res["amplitude"]
            z, _ = self._solve_eom(dt, z0, npts, zgrid, vtotal)
            zj.append(z)
            phi.append(twopi * fs * t)
        return _np.array(zj), _np.array(phi)

    def _solve_eom(self, dt, z0, npts, zgrid, vtotal):
        z = [z0]
        p = [0]
        alpha = self.ring.mom_comp
        for _ in range(npts):
            dp = _np.interp(z[-1], zgrid, vtotal) * dt
            p.append(p[-1] + dp)
            dz = -alpha * _c * p * dt
            z.append(z[-1] + dz)
        return _np.array(z), _np.array(p)

    def calc_synchrotron_frequency_quadratic_potential(self):
        """."""
        ring = self.ring
        nus0 = ring.mom_comp * ring.harm_num
        nus0 *= -ring.gap_voltage * _np.cos(ring.sync_phase)
        nus0 /= 2 * _PI * ring.energy
        nus0 = _np.sqrt(nus0)
        return nus0 * ring.rev_freq

    def calc_synchrotron_frequency_quartic_potential(self, bunch_length):
        """."""
        ring = self.ring
        fs_avg = 2 ** (3 / 4) / _gammafunc(1 / 4) ** 2
        fs_avg *= ring.mom_comp * _c * ring.espread / bunch_length
        fs_std = _np.sqrt((_PI - 2 ** (3 / 2))) / 2 ** (3 / 4) * fs_avg
        return fs_avg, fs_std

    # -------------------- instabilities calculations -------------------------
    def calc_robinson_growth_rate(
        self, w, approx=False, wr=None, Rs=None, Q=None
    ):
        """."""
        alpha = self.ring.mom_comp
        I0 = self.ring.total_current
        E0 = self.ring.energy
        w0 = self.ring.rev_ang_freq
        ws = self.ring.sync_tune * w0
        wp = w + ws
        wn = w - ws
        const = I0 * alpha * w0 / (4 * _PI * ws * E0)
        if approx and None not in {wr, Rs, Q}:
            x = w / wr
            const_approx = const * 4 * ws
            growth = const_approx
            growth *= Rs * Q**2
            growth *= (1 - x**2) * (1 + x**2)
            growth /= x**4 * (1 + Q**2 * (1 / x - x) ** 2) ** 2
        else:
            Zlp = self.get_impedance(w=wp, apply_filter=False)
            Zln = self.get_impedance(w=wn, apply_filter=False)
            growth = const * (wp * Zlp.real - wn * Zln.real)
        return growth

    def calc_tuneshifts_cbi(self, w, m=1, nbun_fill=None, radiation=False):
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
            w=w, Zl=Zl, m=m, inverse=False, full=True
        )

        # Relative tune-shifts must be multiplied by ws
        deltaw *= ring.sync_tune * ring.rev_ang_freq

        ring.num_bun = num_bun
        ring.dampte = dampte
        return deltaw, Zl, wp, interpol_Z, spectrum

    def calc_mode_coupling(
        self,
        w,
        cbmode,
        max_azi=10,
        max_rad=12,
        nbun_fill=None,
        modecoup_matrix=None,
        fokker_matrix=None,
        use_fokker=True,
    ):
        """."""
        ring = self.ring
        num_bun = ring.num_bun
        dampte = ring.dampte

        ring.num_bun = nbun_fill if nbun_fill is not None else ring.harm_num

        if _np.array(w).size == 2:
            Zl = self.get_impedance
        else:
            Zl = self.get_impedance(w=w, apply_filter=False)

        (
            eigenfreq,
            modecoup_matrix,
            fokker_matrix,
        ) = ring.longitudinal_mode_coupling(
            w=w,
            Zl=Zl,
            cbmode=cbmode,
            max_azi=max_azi,
            max_rad=max_rad,
            modecoup_matrix=modecoup_matrix,
            fokker_matrix=fokker_matrix,
            use_fokker=use_fokker,
        )

        # Relative tune-shifts must be multiplied by ws
        eigenfreq *= ring.sync_tune * ring.rev_ang_freq

        ring.num_bun = num_bun
        ring.dampte = dampte
        return eigenfreq, modecoup_matrix, fokker_matrix

    # -------------------- auxiliary methods ----------------------------------
    def _apply_anderson_acceleration(self, dists, niter, tol, m=None, beta=1):
        """."""
        if beta < 0:
            raise Exception("relaxation parameter beta must be positive.")

        xold = dists[-1].ravel()
        xnew = self._ffunc(xold)
        dists.append(xnew)

        m = m or niter
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

        converged = False

        for k in range(niter):
            t0 = _time.time()
            gamma_k = _np.linalg.lstsq(G_k, gnew, rcond=None)[0]
            # tf1 = _time.time()
            # print(f'AndersonLeastSquares: {tf1-t0:.3f}s')
            xold = xnew
            xnew = xold + gnew
            xnew -= mat @ gamma_k
            xnew *= beta
            if beta != 1:
                xnew += (1 - beta) * (xold - X_k @ gamma_k)
            dists.append(xnew)

            gold = gnew
            # tf2 = _time.time()
            # print(f'MatrixMul: {tf2-tf1:.3f}s')
            gnew = self._ffunc(xnew) - xnew
            # tf3 = _time.time()
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
            # print(f'Trapz: {tf-tf3:.3f}s')
            if self.print_flag:
                print(
                    f"Iter.: {k+1:03d}, Dist. Diff.: {diff[idx]:.3e}"
                    + f" (bucket {idx:03d}), E.T.: {tf:.3f}s"
                )
                # print(f"Iter.: {k+1:03d}, E.T.: {tf-t0:.3f}s")
                print("-" * 20)
            if diff[idx] < tol:
                converged = True
                if self.print_flag:
                    print("distribution ok!")
                break
        return dists, converged

    def _create_freqs(self):
        w0 = self.ring.rev_ang_freq
        p = _np.arange(0, self.max_mode)
        return p * w0

    def _ffunc(self, xk):
        """Haissinski operator."""
        # t0 = _time.time()
        xk = self._reshape_dist(xk)
        total_volt = _np.zeros(xk.shape)
        self.beamload_active = _np.zeros(xk.shape)
        idx_wake = self._get_wake_types_idx()
        if idx_wake:
            wake_sources = [self.impedance_sources[idx] for idx in idx_wake]
            # Wakes need to be evaluated one by one
            _func = self.calc_induced_voltage_wake
            for wake in wake_sources:
                if wake.active_passive == ImpedanceSource.ActivePassive.Active:
                    self.beamload_active += _func(wake_source=wake, dist=xk)
                    self._wake_matrix = None
                    self._exp_z = None
                total_volt += _func(wake_source=wake, dist=xk)
                self._wake_matrix = None
                self._exp_z = None

        idx_imp = self._get_impedance_types_idx()
        if idx_imp:
            # Impedances can be summed and calculated once
            mthd = self.impedance_sources[idx_imp[0]].calc_method
            if mthd == ImpedanceSource.Methods.ImpedanceDFT:
                _func = self.calc_induced_voltage_impedance_dft
            elif mthd == ImpedanceSource.Methods.ImpedanceModeSel:
                _func = self.calc_induced_voltage_impedance_mode_selection
            else:
                raise Exception(
                    "Methods must be ImpedanceDFT or ImpedanceModeSel."
                )
            total_volt += _func(dist=xk)

        if not idx_imp and not idx_wake:
            wksrc = self.impedance_sources[0]
            mthd = wksrc.calc_method
            if mthd == ImpedanceSource.Methods.UniformFillAnalytic:
                _func = self.calc_induced_voltage_uniform_filling
            total_volt += _func(wake_source=wksrc, dist=xk)
        # tf1 = _time.time()
        # print(f'CalcIndVoltage: {tf1-t0:.3f}s')
        # tf2 = _time.time()
        total_volt += self.get_generator_voltage()
        # print(f'CalcGenVoltage: {tf2-tf1:.3f}s')
        self.total_voltage = total_volt
        # tf3 = _time.time()
        fxk, _ = self.calc_distributions_from_voltage(total_volt)
        # print(f'CalcDist: {tf3-tf2:.3f}s')
        return fxk.ravel()

    def _get_impedance_types_idx(self):
        """."""
        imp_idx = []
        for idx, imp in enumerate(self.impedance_sources):
            if "impedance" in imp.calc_method_str.lower():
                imp_idx.append(idx)
        return imp_idx

    def _get_wake_types_idx(self):
        """."""
        wake_idx = []
        for idx, imp in enumerate(self.impedance_sources):
            if "wake" in imp.calc_method_str.lower():
                wake_idx.append(idx)
        return wake_idx

    def _feedback_phasor(self):
        ref_amp = self.main_ref_amp
        ref_phase = self.main_ref_phase
        ref_phase += self.main_ref_phase_offset
        wrf = 2 * _PI * self.ring.rf_freq
        phase = wrf * self.zgrid / _c
        dz = _np.diff(self.zgrid)[0]
        vref_phasor = ref_amp * _np.exp(1j * (_PI / 2 - ref_phase))
        vbeamload_phasor = _np.mean(
            _mytrapz(self.beamload_active * _np.exp(1j * phase)[None, :], dz)
        )
        vbeamload_phasor *= 2 / (self.zgrid[-1] - self.zgrid[0])
        vg_phasor = vref_phasor - vbeamload_phasor
        vg = _np.real(vg_phasor * _np.exp(-1j * phase))
        self.main_gen_amp_mon = _np.abs(vg_phasor)
        self.main_gen_phase_mon = _np.angle(vg_phasor)
        return vg[None, :]

    def _feedback_least_squares(self):
        ref_amp = self.main_ref_amp
        ref_phase = self.main_ref_phase
        ref_phase += self.main_ref_phase_offset
        x0 = [ref_amp, ref_phase]
        wrf = 2 * _PI * self.ring.rf_freq
        phase = wrf * self.zgrid / _c
        dz = self.zgrid[1] - self.zgrid[0]

        vref = self.ring.get_voltage_waveform(
            self.zgrid, amplitude=ref_amp, phase=ref_phase
        )
        res = _least_squares(
            fun=self._feedback_err,
            x0=x0,
            args=(phase, dz, self.beamload_active, vref),
            method="lm",
        )
        gen_amp = _np.sqrt(res.x[0] ** 2 + res.x[1] ** 2)
        gen_phase = _np.arctan2(res.x[1], res.x[0])

        self.main_gen_amp_mon = gen_amp
        self.main_gen_phase_mon = gen_phase
        vg = self.ring.get_voltage_waveform(
            self.zgrid, amplitude=gen_amp, phase=gen_phase
        )
        return vg[None, :]

    def _do_zero_padding(self, dist):
        rf_lamb = self.ring.rf_lamb
        dz = _np.diff(self.zgrid)[0]
        # zero-padding
        nr_pts = int(rf_lamb / dz) + 1
        if not nr_pts % 2:
            nr_pts -= 1
        zgrid_full = _np.linspace(-1, 1, nr_pts) * rf_lamb / 2
        dist_new = _np.zeros((dist.shape[0], nr_pts))
        idx_ini = _np.searchsorted(zgrid_full, self.zgrid[0])
        dist_new[:, idx_ini : idx_ini + self.zgrid.size] = dist
        dist = dist_new
        return dist, idx_ini

    @staticmethod
    def _feedback_err(x, *args):
        phase, dz, vbeamload, vref = args
        vgen = LongitudinalEquilibrium._generator_model(phase, x[0], x[1])
        err = (vgen[None, :] + vbeamload) - vref[None, :]
        err = _mytrapz(err * err, dz)
        return err

    @staticmethod
    def _generator_model(phase, a, b):
        return a * _np.sin(phase) + b * _np.cos(phase)

    def _reshape_dist(self, dist):
        return dist.reshape((self.ring.harm_num, self.zgrid.size))
