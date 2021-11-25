#!/usr/bin/env python3

import os as _os
import re as _re
import sh as _sh
import gzip as _gzip
import pickle as _pickle
import logging as _log
import json as _json

import numpy as _np
from scipy import integrate as _scy_int
import matplotlib.pyplot as _plt
from matplotlib import rc as _rc

from mathphys.constants import light_speed as c

from . import sirius as _sirius

try:
    from pyaccel import naff as _naff
    bool_pyaccel = True
except Exception:
    bool_pyaccel = False

_rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# for Palatino and other serif fonts use:
# _rc('font',**{'family':'serif','serif':['Palatino']})
_rc('text', usetex=True)

# # Troughout the code I am assuming:
# s positive means particle behind source -->  Wl, Wt = 0 s < 0
# Wl(s) = -c/Q * int El(ct-s,t) dt
# Wx(s) = - int_-inf^s dWl/dx ds'
# Zl =   int exp(i*w*s) Wl(s) ds
# Zx = i*int exp(i*w*s) Wx(s) ds

_jnpth = _os.path.sep.join
_si = _sirius.create_ring()

DEFAULT_FNAME_SAVE = 'SimulData.pickle'
FNAME_ECHOZ1 = r"wake.dat"
FNAME_ECHOZ2 = r"wake[LT]{1}.dat"
# the older .dat files are not treated:
FNAME_ECHOZR2D = r"(wakeL_([0-9]{2}).txt)"
FNAME_ECHO3D = r'wake3Dindirect.bin'
FNAME_GDFIDL = r"[\w-]+W[YXq]{1}_AT_XY.[0-9]{4}"

ANALYSIS_TYPES = {
    'dx',  # horizontal impedance
    'dy',  # vertical impedance
    'db',  # both planes are symmetric
    'll'   # longitudinal and transverse quadrupolar impedances
    }

PLANES = ('ll', 'dx', 'dy', 'qx', 'qy')
TITLES = {
    'll': 'Longitudinal',
    'dx': 'Dipolar Horizontal',
    'dy': 'Dipolar Vertical',
    'qx': 'Quadrupolar Horizontal',
    'qy': 'Quadrupolar Vertical'}
WAKE_YLABELS = {
    'll': r'$W_l$ [V/pC]',
    'dx': r'$W_{{D_x}}$ [V/pC/m]',
    'dy': r'$W_{{D_y}}$ [V/pC/m]',
    'qx': r'$W_{{Q_x}}$ [V/pC/m]',
    'qy': r'$W_{{Q_y}}$ [V/pC/m]'}
IMPS_YLABELS = {
    'll': r'$Z_l$ [$\Omega$]',
    'dx': r'$Z_{{D_x}}$ [$\Omega$/m]',
    'dy': r'$Z_{{D_y}}$ [$\Omega$/m]',
    'qx': r'$Z_{{Q_x}}$ [$\Omega$/m]',
    'qy': r'$Z_{{Q_y}}$ [$\Omega$/m]'}


# ########################## Public Interface ##########################
class EMSimulData:
    """."""

    def __init__(self, code=None):
        """."""
        # CST, ACE3P, GdfidL, ECHOz1 ECHOz2, ECHOzR, ECHO2D, ECHO3D...
        self.code = code
        # Bunch Length Used in simulation [m]
        self.bunlen = 0.0
        # positions where bunch is defined [m]
        self.sbun = _np.array([], dtype=float)
        # bunch profile used in the simulation [As/m]
        self.bun = _np.array([], dtype=float)
        # axis: distance from following to drive bunch [m]
        self.s = _np.array([], dtype=float)
        # Longitudinal Wakepotential [V/C]
        self.Wll = _np.array([], dtype=float)
        # Dipolar Horizontal Wakepotential [V/C/m]
        self.Wdx = _np.array([], dtype=float)
        # Dipolar Vertical Wakepotential [V/C/m]
        self.Wdy = _np.array([], dtype=float)
        # Quadrupolar Horizontal Wakepotential [V/C/m]
        self.Wqx = _np.array([], dtype=float)
        # Quadrupolar Vertical Wakepotential [V/C/m]
        self.Wqy = _np.array([], dtype=float)
        # axis: frequency obtained from FFT [GHz]
        self.freq = _np.array([], dtype=float)
        # Longitudinal Impedance [Ohm]
        self.Zll = _np.array([], dtype=complex)
        # Dipolar Horizontal Impedance [Ohm]
        self.Zdx = _np.array([], dtype=complex)
        # Dipolar Vertical Impedance [Ohm]
        self.Zdy = _np.array([], dtype=complex)
        # Quadrupolar Horizontal Impedance [Ohm]
        self.Zqx = _np.array([], dtype=complex)
        # Quadrupolar Vertical Impedance [Ohm]
        self.Zqy = _np.array([], dtype=complex)
        self._klossW = None
        self._kckdxW = None
        self._kckdyW = None
        self._kckqxW = None
        self._kckqyW = None

    def copy(self):
        """."""
        other = EMSimulData()
        other.code = self.code
        other.bunlen = self.bunlen
        other.sbun = self.sbun.copy()
        other.bun = self.bun.copy()
        other.s = self.s.copy()
        other.Wll = self.Wll.copy()
        other.Wdx = self.Wdx.copy()
        other.Wdy = self.Wdy.copy()
        other.Wqx = self.Wqx.copy()
        other.Wqy = self.Wqy.copy()
        other.freq = self.freq.copy()
        other.Zll = self.Zll.copy()
        other.Zdx = self.Zdx.copy()
        other.Zdy = self.Zdy.copy()
        other.Zqx = self.Zqx.copy()
        other.Zqy = self.Zqy.copy()
        other._klossW = self._klossW
        other._kckdxW = self._kckdxW
        other._kckdyW = self._kckdyW
        other._kckqxW = self._kckqxW
        other._kckqyW = self._kckqyW
        return other

    def get_klossW(self):
        """."""
        if self._klossW:
            return self._klossW
        sigs, spos = self.bunlen, self.s
        wake = self.Wll
        if wake is None or _np.allclose(wake, 0, atol=0):
            return None
        rhos = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW = _np.trapz(wake*rhos, x=spos)
        self._klossW = kW
        return kW

    def get_PlossW(self, rev_time=None, harm_num=None, Iavg=500e-3):
        """."""
        rev_time = rev_time or _si.T0
        harm_num = harm_num or _si.harm_num
        kW = self.klossW()
        Ploss = kW * Iavg**2 * rev_time * 1e12 / harm_num
        return Ploss

    def get_kick_factorW(self, pl='dy'):
        """."""
        kick = getattr(self, '_kck'+pl+'W')
        if kick:
            return kick
        sigs, spos = self.bunlen, self.s
        wake = getattr(self, 'W'+pl)
        if wake is None or _np.all(wake == 0):
            return None
        rhos = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW = _np.trapz(wake*rhos, x=spos)
        setattr(self, '_kck'+pl+'W', kW)
        return kW

    def get_klossZ(self, bunlen=2.65e-3, nbunches=1):
        """."""
        _si.nbun = nbunches
        klossZ, *_ = _si.loss_factor(
            w=self.freq*2*_np.pi, Zl=self.Zll, bunlen=bunlen)
        return klossZ

    def get_kick_factorZ(self, pl='dy', bunlen=2.65e-3, nbunches=1):
        """."""
        _si.nbun = nbunches
        Z = getattr(self, 'Z'+pl)
        if Z is None or _np.allclose(Z, 0, atol=0):
            return None
        kckZ, *_ = _si.kick_factor(
            w=self.freq*2*_np.pi, Z=Z, bunlen=bunlen)
        return kckZ

    def get_PlossZ(self, bunlen=2.65e-3, nbunches=1):
        """."""
        _si.nbun = nbunches
        _, PlossZ, *_ = _si.loss_factor(
            w=self.freq*2*_np.pi, Zl=self.Zll, bunlen=bunlen)
        return PlossZ


def load_raw_data(
        simul_data=None, code=None, path=None, anal_pl=None, silent=False):
    """."""
    if not simul_data:
        simul_data = EMSimulData()

    _log.basicConfig(level=_log.CRITICAL if silent else _log.INFO)

    if path is None:
        path = _os.path.abspath('.')

    _log.info('#'*60 + '\nLoading Simulation Data')

    # First try to guess the code used in simulation, if not supplied:
    if code is None:
        _log.info('Simulation Code not supplied.')
        code = _get_code(path)
    _log.info(code)
    simul_data.code = code

    # Now try to guess the plane of the analysis:
    if anal_pl is None:
        _log.info('Plane of Analysis not supplied.')
        anal_pl = _get_plane_of_analysis(path, code)
    _log.info(anal_pl)

    # changes in simul_data are made implicitly
    CODES[code](simul_data, path=path, anal_pl=anal_pl)

    _log.info('#'*60+'\n')
    return simul_data


def calc_impedance(
        simul_data, use_win='phase', pl=None, cutoff=2, s_min=None,
        s_max=None, silent=False):
    """."""
    def choose_fft_length(n, max_factor=1000):
        for p in range(n):
            n2, i = n-p, 2
            while (i * i <= n2 and i < max_factor):
                if n2 % i:
                    i += 1
                else:
                    n2 //= i
                if n2 < max_factor:
                    return p, n-p

    def _get_impedance(spos, wake, sigt, cutoff):
        # frequency scale (Hz):
        dt = (spos[-1]-spos[0]) / (spos.shape[0]-1) / c
        # fft == \int exp(-i*2pi*f*t/n) G(t) dt:
        VHat = _np.fft.fft(wake, wake.shape[0]) * dt
        freq = _np.fft.fftfreq(wake.shape[0], d=dt)
        VHat = _np.fft.fftshift(VHat)  # shift the negative frequencies
        freq = _np.fft.fftshift(freq)  # to the center of the spectrum
        # Longitudinal position shift to match center of the bunch with zero z:
        w = 2*_np.pi*freq
        VHat *= _np.exp(-1j*w*(spos[0])/c)
        # find the maximum useable frequency
        wmax = cutoff/sigt
        indcs = _np.abs(w) <= wmax
        # Deconvolve the Transform with a gaussian bunch:
        Jwlist = _np.exp(-(w*sigt)**2/2)
        Z = VHat[indcs]/Jwlist[indcs]
        return freq[indcs], Z

    _log.basicConfig(level=_log.CRITICAL if silent else _log.INFO)
    _log.info('#'*60 + '\n' + 'Calculating Impedances')

    planes = PLANES
    if pl is not None:
        planes = [pl]

    # Extracts Needed Variables
    sigt = simul_data.bunlen / c  # bunch time-length
    spos = simul_data.s.copy()

    if s_min is None:
        s_min = spos[0]
    if s_max is None:
        s_max = spos[-1]
    inds = _np.logical_and(spos >= s_min, spos <= s_max)
    spos = spos[inds]
    # numpy's fft algorithm is slow for large primes:
    p, n = choose_fft_length(spos.shape[0])
    spos = spos[:n]

    if p > 0:
        _log.info(
            'Last {0:d} point{1:s} of wake {2:s} '.format(
                p, *(('s', 'were') if p > 1 else ('', 'was'))) +
            'not considered to gain performance in FFT.')

    if use_win is True:
        _log.info('Using Half-Hanning Window')
        # Half Hanning window to zero the end of the signal
        window = _np.hanning(2*spos.shape[0])[spos.shape[0]:]
    elif isinstance(use_win, str) and use_win.lower().startswith('phase'):
        _log.info('Using Half-Hanning Window to correct the phases')
        # Half Hanning window to smooth the final of the signal
        window = _np.hanning(2*spos.shape[0])[spos.shape[0]:]
    else:
        _log.info('Not using Window')
        window = _np.ones(spos.shape[0])

    _log.info(f'Cutoff frequency w = {cutoff:d}/sigmat')

    for pl in planes:
        _log.info(f'Performing FFT on W{pl:s}: ')
        Wpl = getattr(simul_data, 'W'+pl).copy()
        if Wpl is None or _np.all(Wpl == 0):
            _log.info('No Data found.')
            continue
        _log.info('Data found. ')
        Wpl = Wpl[inds]
        Wpl = Wpl[:n]

        simul_data.freq, Zpl = _get_impedance(spos, Wpl*window, sigt, cutoff)
        if isinstance(use_win, str) and use_win.lower().startswith('phase'):
            _, Zpl2 = _get_impedance(spos, Wpl, sigt, cutoff)
            Zpl = _np.abs(Zpl2)*_np.exp(1j*_np.angle(Zpl))

        if pl == 'll':
            # I have to take the conjugate of the fft because:
            # fftt == \int exp(-i*2pi*f*t/n) G(t) dt
            # while impedance, according to Chao and Ng, is given by:
            # Z == \int exp(i*2pi*f*t/n) G(t) dt
            simul_data.Zll = Zpl.conj()
        else:
            # the Transverse impedance, according to Chao and Ng, is given by:
            # Z == i\int exp(i*2pi*f*t/n) G(t) dt
            setattr(simul_data, 'Z'+pl, 1j*Zpl.conj())
        _log.info('Impedance Calculated.')

    _log.info('#'*60 + '\n')


def calc_impedance_naff(
        simul_data, pl='ll', s_min=None, s_max=None, win=1, nr_ff=20):
    """."""
    if pl not in PLANES:
        raise Exception(
            'Value of variable pl not accepted. Must be one of these: ' +
            ', '.join(PLANES))

    # Extracts Needed Variables
    sigt = simul_data.bunlen / c  # bunch time-length
    spos = simul_data.s.copy()
    W = getattr(simul_data, 'W' + pl).copy()
    sizeW = len(W)
    if W is None or not len(W) or _np.all(W == 0):
        raise Exception('No Data found.')

    if s_min is None:
        s_min = spos[0]
    if s_max is None:
        s_max = spos[-1]

    inds = _np.logical_and(spos >= s_min, spos <= s_max)
    spos = spos[inds]
    W = W[inds]
    dt = (spos[1]-spos[0])/c
    leng = len(W) - (len(W)-1) % 6
    spos = spos[-leng:]
    W = W[-leng:]
    if 0.49 < win < 0.51:
        W *= _np.hanning(2*spos.shape[0])[spos.shape[0]:]
        tu, a = _naff.naff_general(
            W, use_win=0, is_real=False, nr_ff=nr_ff)
    elif isinstance(win, int):
        tu, a = _naff.naff_general(
            W, use_win=win, is_real=False, nr_ff=nr_ff)
    else:
        raise Exception(
            'Win must be 1/2 for half-hanning window or an integer '
            'for other windows(0 --> no window).')

    freq = tu/dt
    w = 2*_np.pi*freq
    # Longitudinal position shift to match center of the bunch with zero z:
    a *= _np.exp(-1j*w*(spos[0])/c)

    # Reconstruct the signal
    S = _np.zeros(leng, dtype=complex)
    for wi, ai in zip(w, a):
        S += ai*_np.exp(1j*wi*spos/c)
    S = S.real

    # Deconvolve the Transform with a gaussian bunch:
    a /= _np.exp(-(w*sigt)**2/2)

    # Must multiply by the vector length due to difference in the meaning
    # of the amplitune in the NAFF transform and the fourier transform
    Z = a*dt*sizeW

    if pl == 'll':
        # I have to take the conjugate of the fft because:
        # fftt == \int exp(-i*2pi*f*t/n) G(t) dt
        # while impedance, according to Chao and Ng, is given by:
        # Z == \int exp(i*2pi*f*t/n) G(t) dt
        Z = Z.conj()
    else:
        # the Transverse impedance, according to Chao and Ng, is given by:
        # Z == i\int exp(i*2pi*f*t/n) G(t) dt
        Z = 1j*Z.conj()
    return freq, Z, leng, S


def save_processed_data(simul_data, silent=False, pth2sv=None):
    """."""
    _log.basicConfig(level=_log.CRITICAL if silent else _log.INFO)
    _log.info('#'*60 + '\nSaving Processed data:')
    spos = simul_data.s
    freq = simul_data.freq

    if pth2sv is None:
        _log.info('Saving in the same folder of the raw data')
        pth2sv = _os.path.abspath('.')
    elif type(pth2sv) is str:
        _log.info('Saving to subfolder: ' + pth2sv)
        if not _os.path.isdir(pth2sv):
            _log.info('Folder does not exist. Creating it...')
            _os.mkdir(pth2sv)
    else:
        msg = 'pth2sv must be a string or None object'
        _log.info(msg)
        raise Exception(msg)

    # Save wakes
    for par in PLANES:
        unit = 'V/C' if par == 'll' else 'V/C/m'
        header = '{0:30s} {1:30s}'.format(
            's [m]', 'W{0:s} [{1:s}]'.format(par, unit))
        fname = _jnpth([pth2sv, 'W'+par+'.gz'])
        wake = getattr(simul_data, 'W'+par)
        if wake is None or _np.all(wake == 0):
            continue
        _log.info('Saving W' + par + ' data to .gz file')
        _np.savetxt(
            fname, _np.array([spos, wake]).transpose(),
            fmt=['%30.16g', '%30.16g'], header=header)

    # Save Impedances
    for par in PLANES:
        unit = 'Ohm' if par == 'll' else 'Ohm/m'
        header = '{0:30s} {1:30s} {2:30s}'.format(
            'Frequency [GHz]',
            'ReZ{0:s} [{1:s}]'.format(par, unit),
            'ImZ{0:s} [{1:s}]'.format(par, unit))
        fname = _jnpth([pth2sv, 'Z'+par+'.gz'])
        Z = getattr(simul_data, 'Z'+par)
        if Z is None or _np.allclose(Z, 0, atol=0):
            continue
        _log.info('Saving Z' + par + ' data to .gz file')
        _np.savetxt(
            fname, _np.array([freq/1e9, Z.real, Z.imag]).transpose(),
            fmt=['%30.16g', '%30.16g', '%30.16g'], header=header)

    _log.info('Saving the Complete EMSimulData structure to a .pickle file.')
    with _gzip.open(_jnpth((pth2sv, DEFAULT_FNAME_SAVE)), 'wb') as f:
        _pickle.dump(simul_data, f, _pickle.HIGHEST_PROTOCOL)

    _log.info('All Data Saved\n' + '#'*60)


def load_processed_data(filename):
    """."""
    with _gzip.open(filename, 'rb') as fh:
        simul_data = _pickle.load(fh)
    return simul_data


def plot_wakes(simul_data, save_figs=False, pth2sv=None, show=False, pls=None):
    """."""
    sbun = simul_data.sbun
    sigs = simul_data.bunlen
    spos = simul_data.s
    if pls is None:
        pls = PLANES
    for pl in pls:
        wake = getattr(simul_data, 'W'+pl)*1e-12  # V/C -> V/pC
        if wake is None or _np.allclose(wake, 0, atol=0):
            continue
        max_wake = wake[_np.abs(wake).argmax()]
        bunchshape = simul_data.bun * (max_wake/simul_data.bun.max())

        f, axs = _plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(14, 6))
        ax = axs[0]
        b = ax.get_position()
        b.x0, b.x1 = 0.05, 0.45
        ax.set_position(b)
        ax.plot(sbun*1000, bunchshape, 'b', linewidth=2, label='Bunch Shape')
        ax.plot(spos*1000, wake, 'r', linewidth=2, label='Wake Potential')
        ax.grid(True)
        ax.set_ylabel(WAKE_YLABELS[pl], fontsize=13)
        ax.set_xlim([spos[0]*1000, 8000*sigs])
        ax.set_ylim([wake.min()*1.1, wake.max()*1.1])
        ax.legend(loc='best')

        ax = axs[1]
        b = ax.get_position()
        b.x0, b.x1 = 0.45, 0.95
        ax.set_position(b)
        ax.plot(spos*1000, wake, 'r', linewidth=2)
        ax.grid(True)
        tit = ax.set_title(TITLES[pl], fontsize=13)
        tit.set_x(0.1)
        xl = ax.set_xlabel('s [mm]', fontsize=13)
        xl.set_x(0.1)
        ax.set_xlim([8000*sigs, spos[-1]*1000])
        ax.set_ylim([wake.min()*1.1, wake.max()*1.1])

        if save_figs:
            f.savefig(_jnpth((pth2sv, 'W'+pl+'.svg')))
    if show:
        _plt.show()


def plot_impedances(
        simul_data, save_figs=False, pth2sv=None, show=False, pls=None):
    """."""
    freq = simul_data.freq
    if pls is None:
        pls = PLANES
    for pl in pls:
        Z = getattr(simul_data, 'Z'+pl)
        if Z is None or _np.allclose(Z, 0, atol=0):
            continue

        _plt.figure()
        _plt.plot(freq/1e9, Z.real, 'r', linewidth=2, label='Re')
        _plt.plot(freq/1e9, Z.imag, 'b--', linewidth=2, label='Im')
        _plt.xlabel('Frequency [GHz]', fontsize=13)
        _plt.grid(True)
        _plt.title(TITLES[pl], fontsize=13)
        _plt.ylabel(IMPS_YLABELS[pl], fontsize=13)
        _plt.legend(loc='best')
        _plt.xlim(freq[[0, -1]]/1e9)
        if save_figs:
            _plt.savefig(_jnpth((pth2sv, 'Z'+pl+'.svg')))
    if show:
        _plt.show()


def plot_losskick_factors(
        simul_data, save_figs=False, pth2sv=None, show=False, pls=None):
    """."""
    # Extracts and Initialize Needed Variables:
    _si.nom_cur = 500e-3
    # bunch length scenarios:
    sigvec = _np.array([2.65, 5, 8, 10, 15], dtype=float)*1e-3
    # current scenarios
    Ivec = _np.linspace(10e-3, _si.nom_cur, num=50)

    bunlen = simul_data.bunlen
    sigi = _np.linspace(bunlen, 18e-3, num=50)
    fill_pat = _np.array([1, 864, 864/2, 864/4], dtype=int)
    if pls is None:
        pls = PLANES
    pls2 = []
    for pl in pls:
        W = getattr(simul_data, 'W'+pl)
        if W is None or _np.allclose(W, 0, atol=0):
            continue
        Z = getattr(simul_data, 'Z'+pl)
        if Z is None or _np.allclose(Z, 0, atol=0):
            continue
        pls2.append(pl)

    for pl in pls2:
        if pl == 'll':
            f, axs = _plt.subplots(
                nrows=1, ncols=2, figsize=(12, 6),
                gridspec_kw=dict(left=0.08, right=0.97))
            ax = axs[0]
            fname = 'Loss_factor'
            for i in range(fill_pat.shape[0]):
                kZi = _np.zeros(sigi.shape[0])
                for j in range(sigi.shape[0]):
                    kZi[j] = simul_data.get_klossZ(
                        bunlen=sigi[j], nbunches=fill_pat[i]) * 1e-12  # V/pC
                ax.semilogy(
                    sigi * 1e3, kZi * 1e3, 'o', markersize=4,
                    label=r'$n = {0:03d}$'.format(fill_pat[i]))
                if not i:
                    kZ = kZi[0]
            # Calculates klossW
            kW = simul_data.get_klossW() * 1e-12
            # Print loss factor calculated in both ways
            ax.semilogy(
                bunlen * 1e3, kW * 1e3, '*', markersize=7, color=[1, 0, 0],
                label=r'$K_L^W$')
            ax.set_title('Loss Factor for $n$ equally spaced bunches.')
            ax.set_xlabel(r'$\sigma$ [mm]')
            ax.set_ylabel(r'$K_L$ [mV/pC]')
            ax.legend(loc='best')
            ax.grid(True)
            ax.annotate(
                r'$K_L^W = {0:5.2f}$ mV/pC'.format(kW*1e3),
                xy=(bunlen*1.1e3, kW*1e3), fontsize=12)
            ax = axs[1]
            kZvec = _np.zeros(sigvec.shape[0])
            labels = []
            for i in range(sigvec.shape[0]):
                kZvec[i] = simul_data.get_klossZ(bunlen=sigvec[i])  # V/C
                labels.append(r'$\sigma = {0:05.2f}$ mm'.format(sigvec[i]*1e3))
            Plossvec = kZvec[None, :] * Ivec[:, None]**2
            Plossvec *= _si.T0 / _si.harm_num

            ax.semilogy(Ivec*1e3, Plossvec, markersize=4)
            ax.set_title(
                'Power Loss for ${0:d}$ equally spaced bunches.'.format(
                    _si.harm_num))
            ax.set_xlabel(r'$I_{{avg}}$ [mA]')
            ax.set_ylabel(r'Power [W]')
            ax.legend(labels, loc='best')
            ax.grid(True)
        else:
            f = _plt.figure(figsize=(6, 6))
            ax = _plt.axes()
            fname = 'Kck'+pl+'_factor'
            for i in range(fill_pat.shape[0]):
                kZi = _np.zeros(sigi.shape[0])
                for j in range(sigi.shape[0]):
                    kZi[j] = simul_data.get_kick_factorZ(
                        pl=pl, bunlen=sigi[j],
                        nbunches=fill_pat[i]) * 1e-12  # V/pC/m
                ax.plot(
                    sigi*1e3, kZi, 'o', markersize=4,
                    label=r'n = {0:03d}'.format(fill_pat[i]))
                if not i:
                    kickZ = kZi[0]

            # Calculates kickW:
            kickW = simul_data.get_kick_factorW(pl=pl) * 1e-12
            # Print loss factor calculated in both ways
            ax.plot(
                bunlen*1e3, kickW, '*', markersize=7, color=[1, 0, 0],
                label=r'$K_{{{0:s}_{1:s}}}^W$'.format(pl[0].upper(), pl[1]))
            ax.set_title('Kick Factor for $n$ equally spaced bunches.')
            ax.set_xlabel(r'$\sigma$ [mm]', fontsize=13)
            ax.set_ylabel(
                r'$K_{{{0:s}_{1:s}}}$ [V/pC/m]'.format(pl[0].upper(), pl[1]),
                fontsize=13)
            ax.legend(loc='best')
            ax.grid(True)
            stri = r'$K_{{{0:s}_{1:s}}}^W = {2:5.2f}$ V/pC/m'.format(
                pl[0].upper(), pl[1], kickW)
            ax.annotate(stri, xy=(bunlen*1.1e3, kickW), fontsize=13)
        if save_figs:
            _plt.savefig(_jnpth((pth2sv, fname+'.svg')))
    if show:
        _plt.show()


def create_make_fig_file(path=None):
    """."""
    if path is None:
        path = _os.path.abspath('.')
    fname = _jnpth([path, 'create_figs.py'])
    analysis = '#!/usr/bin/env python3\n\n'
    analysis += 'import os\n'
    analysis += 'import pycolleff.process_wakes as ems\n\n'
    analysis += 'opts = dict(save_figs=False,show=False)\n'
    analysis += 'path = os.path.abspath(__file__).rpartition(os.path.sep)[0]\n'
    analysis += "file_name = os.path.sep.join([path,'{0:s}'])\n".format(
        DEFAULT_FNAME_SAVE)
    analysis += 'simul_data = ems.load_processed_data(file_name)\n'
    analysis += 'ems.plot_wakes(simul_data,**opts)\n'
    analysis += 'ems.plot_impedances(simul_data,**opts)\n'
    analysis += 'ems.plot_losskick_factors(simul_data,**opts)\n'
    analysis += 'ems.show_now()\n'
    with open(fname, 'w') as f:
        f.writelines(analysis)
    _sh.chmod('+x', fname)


show_now = _plt.show


# ########################## Auxiliary Methods ##########################

def _get_code(path):
    # Split the path to try to guess other parameters:

    _log.info('Trying to guess from path: ')
    path_split = set(path.lower().split(_os.path.sep))
    code_guess = list(CODES.keys() & path_split)
    if code_guess:
        return code_guess[0]
    _log.info('could not be guessed by path.')

    _log.info('Trying to guess from files in folder: ')
    f_in_dir = ' '.join(_os.listdir(path))
    if len(_re.findall(FNAME_GDFIDL, f_in_dir)):
        return 'gdfidl'

    if len(_re.findall(FNAME_ECHOZ1, f_in_dir)):
        return 'echoz1'

    if len(_re.findall(FNAME_ECHOZ2, f_in_dir)):
        return 'echoz2'

    if len(_re.findall(FNAME_ECHO3D, f_in_dir)):
        return 'echo3d'

    f_mat = None
    if len(_re.findall(FNAME_ECHOZR2D, f_in_dir)):
        fol = path
        f_mat = _re.findall(FNAME_ECHOZR2D, f_in_dir)
    elif _os.path.isdir(_jnpth([path, 'elec'])):
        fol = _jnpth([path, 'elec'])
        f_in_fol = ' '.join(_os.listdir(fol))
        f_mat = _re.findall(FNAME_ECHOZR2D, f_in_fol)
    elif _os.path.isdir(_jnpth([path, 'magn'])):
        fol = _jnpth([path, 'magn'])
        f_in_fol = ' '.join(_os.listdir(fol))
        f_mat = _re.findall(FNAME_ECHOZR2D, f_in_fol)

    if f_mat is not None and _os.path.isfile(_jnpth([fol, f_mat[0][0]])):
        with open(_jnpth([fol, f_mat[0][0]])) as f:
            code = 'echozr'
            if f.readline().find('[cm]') <= 0:
                code = 'echo2d'
            return code

    msg = 'Simulation Code was not supplied and could not be guessed.'
    _log.info(msg)
    raise Exception(msg)


def _get_plane_of_analysis(path, code):
    # Split the path to try to guess other parameters:
    path_split = set(path.lower().split(_os.path.sep))

    _log.info('Trying to guess from path: ')
    anal_pl_guess = list(ANALYSIS_TYPES & path_split)
    if anal_pl_guess:
        return anal_pl_guess[0]
    _log.info('could not be guessed by path.')

    _log.info('Trying to guess from files in folder and code: ')
    if code == 'echoz1':
        return 'll'

    if code == 'echoz2':
        return 'dy' if _os.path.isfile('wakeT.dat') else 'll'

    if code == 'gdfidl':
        f_in_dir = ' '.join(_os.listdir(path))
        f_mat = _re.findall(
            r"[\w-]+W([YXq]{2})_AT_XY.[0-9]{4}", f_in_dir)
        if len(f_mat) > 0:
            # f_mat = _re.findall(
            #     r"[\w-]+W([YX]{1})_AT_XY.[0-9]{4}",f_in_dir)
            # countx = [x for x in f_mat if x=='X']
            # county = [y for y in f_mat if y=='Y']
            # anal_pl = 'dy' if len(county) >= len(county) else 'dx'
            anal_pl = 'd'+f_mat[0][0].lower()
        else:
            anal_pl = 'll'
        return anal_pl

    if code == 'echozr':
        if _os.path.isdir(_jnpth([path, 'magn'])):
            return 'll'
        if _os.path.isdir(_jnpth([path, 'elec'])):
            return 'dy'
        if _os.path.isfile('wakeL_01.txt'):
            w = _np.loadtxt(
                'wakeL_01.txt', skiprows=3, usecols=(1, ), unpack=True)
            anal_pl = 'll'
            if _np.allclose(w, 0, atol=0):
                anal_pl = 'dy'
            return anal_pl
        msg = 'Plane of analysis was not supplied '
        msg += 'and could not be guessed.'
        _log.info(msg)
        raise Exception(msg)

    if code == 'echo2d':
        if _os.path.isdir(_jnpth([path, 'magn'])):
            return 'll'
        if _os.path.isdir(_jnpth([path, 'elec'])):
            return 'dy'

        anal_pl = 'dy'
        if _os.path.isfile('wakeL_00.txt'):
            anal_pl = 'll'
        return anal_pl

    if code == 'echo3d':
        if _os.path.isdir(_jnpth([path, 'magn'])):
            return 'll'
        if _os.path.isdir(_jnpth([path, 'elec'])):
            return 'dy'

    msg = 'Plane of analysis was not supplied '
    msg += 'and could not be guessed.'
    _log.info(msg)
    raise Exception(msg)


def _ACE3P_load_data(simpar):
    raise NotImplementedError('This function was not tested yet.')
    nsigmas = 5
    headerL = 3
    if wdir.startswith(tardir):
        cpfile = False
    else:
        cpfile = True

    wakepath = _jnpth([wdir, 'wakefield.out'])
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile:
        _sh.cp(wakepath, tardir)

    spos = loadres[:, 0]
    # I know this is correct for ECHO (2015/08/27):
    if m == 0:
        wake = -loadres[:, 1]
    else:
        wake = loadres[:, 1]

    spos = spos - nsigmas * bunlen  # Performs displacement over s axis
    return spos, wake


def _CST_load_data(simpar):
    raise NotImplementedError('This function was not tested yet.')
    headerL = 2
    if wdir.startswith(tardir):
        cpfile = False
    else:
        cpfile = True

    wakepath = _jnpth([wdir, 'wake.txt'])
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile:
        _sh.cp(wakepath, tardir)

    spos = loadres[:, 0]
    wake = loadres[:, 1]
    # Adjust s-axis (rescale or shift)
    spos = spos/1000  # Rescaling mm to m
    if m > 0:
        wake = -wake
    return spos, wake


def _ECHOz1_load_data(simul_data, path, anal_pl):

    if anal_pl == 'll':
        _log.info('Loading longitudinal Wake file:')
        fname = _jnpth([path, FNAME_ECHOZ1])
        if _os.path.isfile(fname):
            _log.info('Data found.')
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            _log.info('Not found.')
            raise Exception('Longitudinal wake file not found.')
    else:
        msg = 'ECHOz1 only calculates longitudinal wake.'
        _log.info(msg)
        raise Exception(msg)

    simul_data.s = loadres[:, 0]/100  # Rescaling cm to m
    # V/C/m (the minus sign is due to convention):
    simul_data.Wll = -loadres[:, 1] * 1e12

    # loading driving bunch info
    _log.info('Loading bunch length from wake.dat')
    sbun = simul_data.s.copy()
    ds = sbun[1] - sbun[0]
    bunlen = abs(sbun[0]-ds/2) / 5
    a = _np.argmin(_np.abs(sbun + sbun[0])) + 1
    sbun = sbun[:a]
    simul_data.bunlen = bunlen
    simul_data.sbun = sbun
    simul_data.bun = _np.exp(-sbun**2/(2*bunlen**2))
    simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen
    _log.info(f'Bunch length of the driving bunch: {bunlen*1e3:7.3g} mm')
    _log.info('Data Loaded.')


def _ECHOz2_load_data(simul_data, path, anal_pl):

    anal_pl_ori = None
    if anal_pl == 'db':
        anal_pl_ori = 'db'
        anal_pl = 'dy'
        _log.info(
            'Even though there is symmetry, I am loading data to the Y plane.')

    if anal_pl == 'll':
        _log.info('Loading longitudinal Wake file:')
        fname = _jnpth([path, 'wakeL.dat'])
        if _os.path.isfile(fname):
            _log.info('Data found.')
            spos, wl = _np.loadtxt(
                fname, skiprows=0, usecols=(0, 1), unpack=True)
        else:
            _log.info('Not found.')
            Exception('Longitudinal wake file not found.')

        simul_data.s = spos/100  # Rescaling cm to m
        simul_data.Wll = -wl * 1e12  # V/C (minus sign is due to convention)
    elif anal_pl in {'dx', 'dy'}:
        fname = _jnpth([path, 'wakeL.dat'])
        if _os.path.isfile(fname):
            _log.info(
                'Calculating Transverse wake from longitudinal wake file:')
            _log.info('Data found.')
            spos, wl = _np.loadtxt(
                fname, skiprows=0, usecols=(0, 1), unpack=True)
            simul_data.s = spos/100  # Rescaling cm to m
            # one minus sign due to convention and
            # the other due to Panofsky-Wenzel:
            wt = -_scy_int.cumtrapz(-wl, x=spos/100, initial=0)
            setattr(simul_data, 'W'+anal_pl, wt * 1e12)  # V/C/m
        else:
            _log.info('File not found.')
            _log.info('Loading transverse wake from transverse wake file.:')
            fname = _jnpth([path, 'wakeT.dat'])
            if _os.path.isfile(fname):
                _log.info('Data found.')
                _log.info(
                    'Depending on the ECHOz2 program version this '
                    'may lead to inacurate results.')
                spos, wt = _np.loadtxt(
                    fname, skiprows=0, usecols=(0, 1), unpack=True)
            else:
                _log.info('Not found.')
                Exception('Transverse wake file not found.')
            simul_data.s = spos/100  # Rescaling cm to m
            # there is an error in the integration of echoz2.
            # It is needed to subtract the first value to correct and offset
            # wt = -_scy_int.cumtrapz(-wl, x=spos/100, initial=0)
            setattr(simul_data, 'W'+anal_pl, (wt - wt[0]) * 1e12)
            # V/C/m (minus sign is due to convention)
    else:
        msg = f'Plane of analysis {anal_pl:s} does not match any '
        msg += 'of the possible options'
        _log.info(msg)
        raise Exception(msg)

    # loading driving bunch info
    _log.info('Loading bunch length from wake file')
    sbun = simul_data.s.copy()
    ds = sbun[1]-sbun[0]
    bunlen = abs(sbun[0] - ds/2) / 5
    a = _np.argmin(_np.abs(sbun + sbun[0])) + 1
    sbun = sbun[:a]
    simul_data.bunlen = bunlen
    simul_data.sbun = sbun
    simul_data.bun = _np.exp(-sbun**2/(2*bunlen**2))
    simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen
    _log.info(
        f'Bunch length of the driving bunch: {simul_data.bunlen*1e3:7.3g} mm')
    _log.info('Data Loaded.')

    if anal_pl_ori:
        anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
        _log.info(
            'There is symmetry. Copying the data from the ' +
            f'{anal_pl[1].upper():s} plane to ' +
            f'the {anal_pl_comp[1].upper():s} plane')
        setattr(
            simul_data, 'W'+anal_pl_comp,
            getattr(simul_data, 'W' + anal_pl).copy())


def _ECHOzR_load_data(simul_data, path, anal_pl):
    _ECHO_rect_load_data(simul_data, 'echozr', path, anal_pl)


def _ECHO2D_load_data(simul_data, path, anal_pl):

    _log.info('Trying to find out the geometry type: ')

    if (_os.path.isdir(_jnpth([path, 'magn'])) or
            _os.path.isdir(_jnpth([path, 'elec']))):
        geo_type = 'rectangular'
    elif (_os.path.isfile(_jnpth([path, 'wakeL_00.txt'])) or
            _os.path.isfile(_jnpth([path, 'wakeL_01.txt']))):
        geo_type = 'round'
    else:
        msg = 'Could not find out the geometry type.'
        _log.info(msg)
        raise Exception(msg)
    _log.info(geo_type)

    if geo_type == 'rectangular':
        _ECHO_rect_load_data(simul_data, 'echo2d', path, anal_pl)
    else:
        anal_pl_ori = None
        if anal_pl == 'db':
            anal_pl_ori = 'db'
            anal_pl = 'dy'
            _log.info(
                'Even though there is symmetry, '
                'I am loading data to the Y plane.')

        if anal_pl == 'll':
            _log.info('Loading longitudinal Wake file:')
            fname = _jnpth([path, 'wakeL_00.txt'])
            if _os.path.isfile(fname):
                _log.info('Data found.')
                with open(fname) as f:
                    f.readline()
                    mstep, offset = _np.fromstring(f.readline(), sep='\t')
                    f.readline()
                    _, bunlen = _np.fromstring(f.readline(), sep='\t')
                spos, Wm = _np.loadtxt(fname, skiprows=6, unpack=True)
                simul_data.s = spos
                simul_data.Wll = -Wm  # V/C the minus sign is due to convention
            else:
                _log.info('Not found.')
                Exception('Longitudinal wake file not found.')
        elif anal_pl in {'dx', 'dy'}:
            _log.info('Loading Transverse Wake file:')
            fname = _jnpth([path, 'wakeL_01.txt'])
            if _os.path.isfile(fname):
                _log.info('Data found.')
                with open(fname) as f:
                    f.readline()
                    mstep, offset = _np.fromstring(f.readline(), sep='\t')
                    f.readline()
                    _, bunlen = _np.fromstring(f.readline(), sep='\t')
                # transverse wakes are calculated in the middle of the mesh
                y0 = mstep*(offset+0.5)
                # m and V/C/m^2
                spos, Wm = _np.loadtxt(fname, skiprows=6, unpack=True)
                simul_data.s = spos
                # V/C/m the minus sign is due to convention
                Wdm = -_scy_int.cumtrapz(-Wm/(y0*y0), x=spos, initial=0)
                setattr(simul_data, 'W'+anal_pl, Wdm)
            else:
                _log.info('File not found.')
                Exception('Transverse wake file not found.')
        else:
            _log.info(
                f'Plane of analysis {anal_pl:s} does not match '
                'any of the possible options')
            raise Exception(
                f'Plane of analysis {anal_pl:s} does not match '
                'any of the possible options')

        if anal_pl_ori:
            anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
            _log.info(
                'There is symmetry. Copying the data from the ' +
                f'{anal_pl[1].upper():s} plane to' +
                f' the {anal_pl_comp[1].upper():s} plane')
            setattr(
                simul_data, 'W'+anal_pl_comp,
                getattr(simul_data, 'W'+anal_pl).copy())

        a = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        simul_data.bunlen = bunlen
        simul_data.sbun = sbun
        simul_data.bun = _np.exp(-sbun**2/(2*bunlen**2))
        simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen
        _log.info(
            'Bunch length of the driving bunch: ' +
            f'{simul_data.bunlen*1e3:7.3g} mm')
        _log.info(f'Mesh size used in the simulation: {mstep*1e6:7.4g} um')
        _log.info('Data Loaded.')


def _ECHO3D_load_wakes(path):
    fname = _jnpth([path, 'Results', FNAME_ECHO3D])
    nz, nx, ny = _np.fromfile(fname, dtype=_np.int32, count=3)
    data = _np.fromfile(fname, dtype=_np.float64, offset=4*3)
    x = data[:nx] * 1e-2  # from cm to m
    y = data[nx:nx+ny] * 1e-2  # from cm to m
    wake = data[nx+ny:]
    wake = wake.reshape(nz, nx, ny)
    wake *= -1  # the minus sign is due to convention
    wake *= 1e12  # from V/pC to V/C

    # find the index with the origin
    origx = _np.isclose(x, 0.0).nonzero()[0]
    if not origx.size:
        _log.warning(
            'Attention, wake was not calculated at origin X. '
            'Results might not be accurate.')
        origx = _np.argmin(_np.abs(x))
    else:
        origx = origx[0]
    origy = _np.isclose(y, 0.0).nonzero()[0]
    if not origy.size:
        _log.warning(
            'Attention, wake was not calculated at origin Y. '
            'Results might not be accurate.')
        origy = _np.argmin(_np.abs(y))
    else:
        origy = origy[0]
    return wake, origx, origy


def _ECHO3D_load_info(path):
    with open(_jnpth([path, 'input.txt']), 'r') as fil:
        data = fil.readlines()
    dic = dict()
    for line in data:
        if '=' not in line:
            continue
        lin = line.split('%')[0]
        var, val = lin.split('=')
        var = var.replace(' ', '')
        val = val.strip()
        if "'" in val:
            dic[var] = val.replace("'", '')
        else:
            valo = ''
            while val != valo:
                valo = val
                val = valo.replace('  ', ' ')
            val = val.replace(' ', ',')
            dic[var] = _json.loads(val)
    return dic


def _ECHO3D_find_files(path):
    path_ = _jnpth([path, 'Results'])
    f_match = None
    if _os.path.isdir(path_):
        f_in_dir = ' '.join(_os.listdir(path_))
        f_match = _re.findall(FNAME_ECHO3D, f_in_dir)
    return f_match


def _ECHO3D_load_data(simul_data, path, anal_pl):
    # list all files that match the name pattern for wakefields
    f_match = _ECHO3D_find_files(path)

    if anal_pl == 'll':
        if not f_match:
            msg = 'No files found for longitudinal analysis.'
            _log.info(msg)
            raise Exception(msg)

        wake, origx, origy = _ECHO3D_load_wakes(path)
        nz = wake.shape[0]

        dic = _ECHO3D_load_info(path)
        bunlen = dic['BunchSigma']*1e-3  # in m
        mstepz = dic['Steps'][0]*1e-3  # in m
        mstepx = dic['Steps'][1]*1e-3  # in m
        mstepy = dic['Steps'][2]*1e-3  # in m
        bunx = dic['BunchPosition'][0]  # in step index
        buny = dic['BunchPosition'][1]  # in step index

        if bunx != origx:
            _log.warning(
                'Attention, Bunch was not passed at origin X for '
                'longitudinal wake simulation.'
                'Results might not be accurate.')
        if buny != origy:
            _log.warning(
                'Attention, Bunch was not passed at origin Y for '
                'longitudinal wake simulation.'
                'Results might not be accurate.')

        spos = mstepz*_np.arange(nz, dtype=float)
        spos -= 5 * bunlen
        a = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        bun = _np.exp(-sbun**2/(2*bunlen*bunlen))
        bun /= _np.sqrt(2*_np.pi)*bunlen

        simul_data.s = spos
        simul_data.bunlen = bunlen
        simul_data.sbun = sbun
        simul_data.bun = bun
        _log.info(
            'Bunch length of the driving bunch: ' +
            f'{simul_data.bunlen*1e3:7.3g} mm')
        _log.info(f'Mesh size used in the simulation: {mstepz*1e6:7.4g} um')

        # Load longitudinal Wake
        simul_data.Wll = wake[:, origx, origy]

        # Horizontal quadrupolar wake
        _log.info('Loading Horizontal Quadrupolar Wake file:')
        # Panofky-Wenzel:
        # dWx/dz = -dWl(x, y, z)/dx
        # Wx = -int(dWl(x, y, zp)/dx, zp=-inf, z)
        wakex = _np.gradient(wake[:, :, origy], mstepx, axis=1)
        wakex = -_scy_int.cumtrapz(wakex, x=spos, initial=0, axis=0)
        # Isolate the quadrupolar wake:
        # Wqx = dWx(z)/dx
        simul_data.Wqx = _np.gradient(wakex, mstepx, axis=1)[:, origx]  # V/C/m

        # Vertical quadrupolar wake
        _log.info('Loading Vertical Quadrupolar Wake file:')
        # Panofky-Wenzel:
        # dWy/dz = -dWl(x, y, z)/dy
        # Wy = -int(dWl(x, y, zp)/dy, zp=-inf, z)
        wakey = _np.gradient(wake[:, origx, :], mstepy, axis=1)
        wakey = -_scy_int.cumtrapz(wakey, x=spos, initial=0, axis=0)
        # Isolate the quadrupolar wake:
        # Wqy = dWy(z)/dy
        simul_data.Wqy = _np.gradient(wakey, mstepy, axis=1)[:, origy]  # V/C/m
        _log.info('Longitudinal Data Loaded.')
        return

    anal_pl_ori = None
    if anal_pl == 'db':
        anal_pl_ori = 'db'
        if not f_match:
            anal_pl = 'dx' if _os.path.isdir('dxdpl') else 'dy'
        else:
            msg = (
                "I know this simulation works for both planes, but I " +
                "couldn't find out which plane was in fact used.")
            _log.info(msg)
            raise Exception(msg)
        _log.info(
            'There is symmetry y=x, calculation performed in the ' +
            anal_pl[1].upper() + ' plane.')

    if anal_pl not in {'dx', 'dy'}:
        msg = f'Plane of analysis {anal_pl:s} does not match any of '
        msg += 'the possible options'
        _log.info(msg)
        raise Exception(msg)

    elec_symm = False
    if not f_match:
        _log.info('There is no wake files in this folder.')
        elec_fol = _jnpth([path, 'elec'])
        if _os.path.isdir(elec_fol):
            _log.info(
                ' I found a folder named "elec". I will assume the '
                'simulation has this symmetry.')
            f_match = _ECHO3D_find_files(elec_fol)
            elec_symm = True

    if not f_match:
        _log.info(' I will assume there is no symmetry.')

        sposs, wakes, sbuns, buns, bunlens, xd, yd = [], [], [], [], [], [], []
        origs, msteps = [], []
        for sub_fol in ['dpl', 'dmi']:
            ext_path = _jnpth([path, anal_pl+sub_fol])
            _log.info('Looking for '+anal_pl+sub_fol+' subfolder:')
            if not _os.path.isdir(ext_path):
                _log.info(
                    'For non-symmetric structures, there must '
                    f'be subfolders {anal_pl:s}dpl {anal_pl:s}dmi '
                    'with the data')
                raise Exception('Files not found')
            # list all files that match the pattern
            f_match = _ECHO3D_find_files(ext_path)
            if not f_match:
                msg = 'No files found for transverse analysis.'
                _log.info(msg)
                raise Exception(msg)

            _log.info('Loading wake file and calulating {0:s} wake:'.format(
                'Horizontal' if anal_pl == 'dx' else 'Vertical'))
            wake, origx, origy = _ECHO3D_load_wakes(ext_path)
            nz = wake.shape[0]

            dic = _ECHO3D_load_info(ext_path)
            bunlen = dic['BunchSigma']*1e-3  # in m
            mstepz = dic['Steps'][0]*1e-3  # in m
            mstepx = dic['Steps'][1]*1e-3  # in m
            mstepy = dic['Steps'][2]*1e-3  # in m
            bunx = dic['BunchPosition'][0]  # in step index
            buny = dic['BunchPosition'][1]  # in step index

            spos = mstepz*_np.arange(nz, dtype=float)
            spos -= 5 * bunlen
            a = _np.argmin(_np.abs(spos + spos[0])) + 1
            sbun = spos[:a]
            bun = _np.exp(-sbun**2/(2*bunlen*bunlen))
            bun /= _np.sqrt(2*_np.pi)*bunlen

            sposs.append(spos)
            buns.append(bun)
            sbuns.append(sbun)
            bunlens.append(bunlen)
            xd.append(bunx)
            yd.append(buny)

            # Panofky-Wenzel:
            # dWt/dz = -dWl(t, z)/dt
            # Wt = -int(dWl(t, zp)/dt, zp=-inf, z)
            # with t in {x, y}
            orig = origx if anal_pl == 'dx' else origy
            origo = origy if anal_pl == 'dx' else origx
            mstep = mstepx if anal_pl == 'dx' else mstepy
            bunp = bunx if anal_pl == 'dx' else buny
            origs.append(orig)
            msteps.append(mstep)

            waket = wake.take(origo, axis=2 if anal_pl == 'dx' else 1)
            waket = _np.gradient(waket, mstep, axis=1)
            waket = -_scy_int.cumtrapz(waket, x=spos, initial=0, axis=0)
            wakes.append(waket[:, orig] / (mstep*(bunp-orig)))

        # If the simulation is not ready yet the lenghts may differ.
        # This line is used to truncate the longer wake in the length of
        # the shorter:
        l1 = min(len(sposs[0]), len(sposs[1]))

        bunp = xd if anal_pl == 'dx' else yd
        ndel = yd if anal_pl == 'dx' else xd
        if not (_np.allclose(sposs[0][:l1], sposs[1][:l1], atol=0) and
                _np.allclose(sbuns[0], sbuns[1], atol=0) and
                _np.allclose(buns[0], buns[1], atol=0) and
                _np.allclose(ndel[0], ndel[1], atol=0)):
            msg = 'There is a mismatch between the paramenters of the'
            msg += f'simulation in the {anal_pl:s}dpl and {anal_pl:s}dmi '
            msg += 'folders.'
            _log.info(msg)
            raise Exception(msg)
        simul_data.s = sposs[0][:l1]
        simul_data.bun = buns[0]
        simul_data.sbun = sbuns[0]
        simul_data.bunlen = bunlens[0]

        waket = (wakes[0][:l1] + wakes[1][:l1]) / 2
        setattr(simul_data, 'W'+anal_pl, waket)  # V / pC / m
    else:
        if elec_symm:
            path = elec_fol

        wake, origx, origy = _ECHO3D_load_wakes(path)
        nz = wake.shape[0]

        dic = _ECHO3D_load_info(path)
        bunlen = dic['BunchSigma']*1e-3  # in m
        mstepz = dic['Steps'][0]*1e-3  # in m
        mstepx = dic['Steps'][1]*1e-3  # in m
        mstepy = dic['Steps'][2]*1e-3  # in m
        bunx = dic['BunchPosition'][0]  # in step index
        buny = dic['BunchPosition'][1]  # in step index

        spos = mstepz*_np.arange(nz, dtype=float)
        spos -= 5 * bunlen
        a = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        bun = _np.exp(-sbun**2/(2*bunlen*bunlen))
        bun /= _np.sqrt(2*_np.pi)*bunlen

        simul_data.s = spos
        simul_data.bun = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        # Panofsky-Wenzel:
        # dWt/dz = -dWl(t, z)/dt
        # Wt = -int(dWl(t, zp)/dt, zp=-inf, z)
        # with t in {x, y}
        orig = origx if anal_pl == 'dx' else origy
        origo = origy if anal_pl == 'dx' else origx
        mstep = mstepx if anal_pl == 'dx' else mstepy
        bunp = bunx if anal_pl == 'dx' else buny

        axis = 1 if anal_pl == 'dx' else 2

        waket = wake.take(origo, axis=axis)
        waket = _np.gradient(waket, mstep, axis=1)[:, orig]
        waket = -_scy_int.cumtrapz(waket, x=spos, initial=0, axis=0)

        mstep *= 2 if elec_symm else 1
        waket /= mstep * (bunp-orig)
        setattr(simul_data, 'W'+anal_pl, waket)  # V/C/m
    _log.info('Transverse Data Loaded.')

    if anal_pl_ori:
        anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
        _log.info(
            'There is symmetry. Copying the data from the ' +
            f'{anal_pl[1].upper():s} plane to the ' +
            f'{anal_pl_comp[1].upper():s} plane')
        setattr(
            simul_data, 'W'+anal_pl_comp,
            getattr(simul_data, 'W'+anal_pl).copy())


def _GdfidL_load_data(simul_data, path, anal_pl):
    # list all the files that match the name pattern for wakefields
    f_in_dir = ' '.join(_os.listdir(path))
    f_match = _re.findall(FNAME_GDFIDL, f_in_dir)

    if anal_pl == 'll':
        if not f_match:
            msg = 'No files found for longitudinal analysis.'
            _log.info(msg)
            raise Exception(msg)

        # Load longitudinal Wake
        spos, wake, sbun, bun, bunlen, xd, yd = _GdfidL_get_longitudinal_info(
            path, f_match, pl='ll')
        simul_data.Wll = wake
        simul_data.s = spos
        simul_data.bun = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        # And quadrupolar Wakes, if existent:
        _log.info('Loading Horizontal Quadrupolar Wake file:')
        wake = _GdfidL_get_transversal_info(
            path, f_match, pl='qx')  # V/C/m
        if wake is not None:
            simul_data.Wqx = wake

        _log.info('Loading Vertical Quadrupolar Wake file:')
        wake = _GdfidL_get_transversal_info(
            path, f_match, pl='qy')  # V/C/m
        if wake is not None:
            simul_data.Wqy = wake
        _log.info('Longitudinal Data Loaded.')
        return

    anal_pl_ori = None
    if anal_pl == 'db':
        anal_pl_ori = 'db'
        if not f_match:
            anal_pl = 'dx' if _os.path.isdir('dxdpl') else 'dy'
        else:
            cond = [f for f in f_match if f.find('WX_AT_XY') >= 0]
            anal_pl = 'dx' if cond else 'dy'
        _log.info(
            'There is symmetry y=x, calculation performed in the ' +
            anal_pl[1].upper() + ' plane.')

    if anal_pl not in {'dx', 'dy'}:
        _log.info(
            f'Plane of analysis {anal_pl:s} does not match any of '
            'the possible options')
        raise Exception(
            f'Plane of analysis {anal_pl:s} does not match any of '
            'the possible options')

    elec_symm = False
    if not f_match:
        _log.info('There is no wake files in this folder.')
        elec_fol = _jnpth([path, 'elec'])
        if _os.path.isdir(elec_fol):
            _log.info(
                ' I found a folder named "elec". I will assume the '
                'simulation has this symmetry.')
            f_in_dir = ' '.join(_os.listdir(elec_fol))
            f_match = _re.findall(FNAME_GDFIDL, f_in_dir)
            elec_symm = True

    if not f_match:
        _log.info(' I will assume there is no symmetry.')

        spos, wake, sbun, bun, bunlen, xd, yd = [], [], [], [], [], [], []
        for sub_fol in ['dpl', 'dmi']:
            ext_path = _jnpth([path, anal_pl+sub_fol])
            _log.info('Looking for '+anal_pl+sub_fol+' subfolder:')
            if not _os.path.isdir(ext_path):
                _log.info(
                    'For non-symmetric structures, there must '
                    f'be subfolders {anal_pl:s}dpl {anal_pl:s}dmi '
                    'with the data')
                raise Exception('Files not found')
            # list all the files that match the pattern
            f_in_dir = ' '.join(_os.listdir(ext_path))
            f_match = _re.findall(FNAME_GDFIDL, f_in_dir)
            if not f_match:
                msg = 'No files found for transverse analysis.'
                _log.info(msg)
                raise Exception(msg)

            sp, _, sb, bn, bnln, xdi, ydi = _GdfidL_get_longitudinal_info(
                ext_path, f_match, pl=anal_pl)
            spos.append(sp)
            bun.append(bn)
            sbun.append(sb)
            bunlen.append(bnln)
            xd.append(xdi)
            yd.append(ydi)

            _log.info('Loading {0:s} Dipolar Wake file:'.format(
                'Horizontal' if anal_pl == 'dx' else 'Vertical'))
            # V/C
            wk = _GdfidL_get_transversal_info(
                ext_path, f_match, pl=anal_pl)
            if wk is not None:
                wake.append(wk)
            else:
                _log.info(
                    'Actually there is something wrong, '
                    'these wake files should be here.')
                raise Exception(
                    'Transverse {0:s} dipolar wake files not found'.format(
                        'Horizontal' if anal_pl == 'dx' else 'Vertical'))

        # If the simulation is not ready yet the lenghts may differ.
        # This line is used to truncate the longer wake in the length of
        # the shorter:
        l1 = min(len(spos[0]), len(spos[1]))

        delta = xd if anal_pl == 'dx' else yd
        ndel = yd if anal_pl == 'dx' else xd
        if not (_np.allclose(spos[0][:l1], spos[1][:l1], atol=0) and
                _np.allclose(sbun[0], sbun[1], atol=0) and
                _np.allclose(bun[0], bun[1], atol=0) and
                _np.allclose(ndel[0], ndel[1], atol=0)):
            msg = 'There is a mismatch between the paramenters of the'
            msg += f'simulation in the {anal_pl:s}dpl and {anal_pl:s}dmi '
            msg += 'folders.'
            _log.info(msg)
            raise Exception(msg)
        simul_data.s = spos[0][:l1]
        simul_data.bun = bun[0]
        simul_data.sbun = sbun[0]
        simul_data.bunlen = bunlen[0]
        setattr(
            simul_data, 'W'+anal_pl,
            (wake[0][:l1]-wake[1][:l1])/(delta[0]-delta[1]))  # V/C/m
    else:
        if elec_symm:
            path = elec_fol
        spos, wake, sbun, bun, bunlen, xd, yd = \
            _GdfidL_get_longitudinal_info(
                path, f_match, pl=anal_pl)
        simul_data.s = spos
        simul_data.bun = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        _log.info('Loading {0:s} Dipolar Wake file:'.format(
            'Horizontal' if anal_pl == 'dx' else 'Vertical'))
        wake = _GdfidL_get_transversal_info(
            path, f_match, pl=anal_pl)  # V/C
        if wake is not None:
            delta = xd if anal_pl == 'dx' else yd
            delta *= 2 if elec_symm else 1
            setattr(simul_data, 'W'+anal_pl, wake/delta)  # V/C/m
        else:
            _log.info(
                'Actually there is something wrong, '
                'these wake files should be here.')
            raise Exception(
                'Transverse {0:s} dipolar wake files not found'.format(
                    'Horizontal' if anal_pl == 'dx' else 'Vertical'))
    _log.info('Transverse Data Loaded.')

    if anal_pl_ori:
        anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
        _log.info(
            'There is symmetry. Copying the data from the ' +
            f'{anal_pl[1].upper():s} plane to the ' +
            f'{anal_pl_comp[1].upper():s} plane')
        setattr(
            simul_data, 'W'+anal_pl_comp,
            getattr(simul_data, 'W'+anal_pl).copy())


CODES = {
    'echoz1': _ECHOz1_load_data,
    'echoz2': _ECHOz2_load_data,
    'echo2d': _ECHO2D_load_data,
    'echo3d': _ECHO3D_load_data,
    'echozr': _ECHOzR_load_data,
    'gdfidl': _GdfidL_load_data,
    'ace3p': _ACE3P_load_data,
    'cst': _CST_load_data}


def _GdfidL_load_dados_info(filename):
    dados, info = [], []
    with open(filename) as fh:
        data = fh.read()
    for line in data.splitlines():
        if not line.startswith((' #', ' %', ' $')):
            dados.append(line)
        else:
            info.append(line)
    return dados, info


def _GdfidL_get_charge(info):
    for line in info:
        if line.find('total charge') >= 0:
            lin = line.split(',')[1]
            charge = float(_re.findall(r'[-+]?\d+\.?\d+[eE]?[-+]?\d+', lin)[0])
            break
    return charge


def _GdfidL_get_integration_path(info):
    for line in info:
        if line.find('subtitle=') >= 0:
            x, y = (float(val) for val in _re.findall(
                r'[-+]?\d+\.?\d+[eE]?[-+]?\d+', line))
            break
    return x, y


def _GdfidL_get_longitudinal_info(path, filelist, pl='ll'):
    _log.info('Loading longitunal Wake file:')
    fn = [f for f in filelist if f.find('Wq_AT_XY') >= 0]
    if not fn:
        msg = 'No longitudinal wake file found. It is needed to have one'
        _log.info(msg)
        raise Exception(msg)
    if len(fn) > 1:
        msg = 'More than one longitudinal wake file found. Only 1 is allowed'
        _log.info(msg)
        raise Exception(msg)

    dados, info = _GdfidL_load_dados_info(_jnpth([path, fn[0]]))
    charge = _GdfidL_get_charge(info)
    xd, yd = _GdfidL_get_integration_path(info)
    spos, wake = _np.loadtxt(dados, unpack=True)  # dados is a list of strings
    _log.info(f'Charge of the driving bunch: {charge*1e12:5.3g} pC')
    if pl == 'll' and (abs(xd) > 1e-10 or abs(yd) > 1e-10):
        _log.info(
            'Driving particle not in the origin. '
            'Are you sure this is what you want?')
    elif pl != 'll' and abs(xd) < 1e-10 and abs(yd) < 1e-10:
        _log.info(
            'The driving bunch is too close to origin. '
            'Are you sure this is what you want?')

    a = _np.argmin(_np.diff(spos)) + 1
    sbun = spos[a:]
    bun = wake[a:]*charge/_np.trapz(wake[a:], x=sbun)  # C
    wake = -wake[:a]/charge  # V/C (minus sign because of convention)
    spos = spos[:a]  # m
    bunlen = -sbun[0]/6  # gdfidl uses a bunch with 6-sigma
    _log.info(f'Bunch length of the driving bunch: {bunlen*1e3:7.4g} mm')
    return spos, wake, sbun, bun, bunlen, xd, yd


def _GdfidL_get_transversal_info(path, filelist, pl='qx'):
    stri = 'W{0:s}_AT_XY'.format(pl[1].upper())
    fn = [f for f in filelist if f.find(stri) >= 0]
    if not fn:
        _log.info(f'No W{pl:s} wake file found. Skipping to next')
        return None
    _log.info(f"{len(fn):2d} W{pl:s} wake file found: {', '.join(fn):s}")

    dados, info = _GdfidL_load_dados_info(_jnpth([path, fn[0]]))
    charge = _GdfidL_get_charge(info)
    if pl[1] == 'x':
        delta1, _ = _GdfidL_get_integration_path(info)
    else:
        _, delta1 = _GdfidL_get_integration_path(info)
    _, wake1 = _np.loadtxt(dados, unpack=True)
    _log.info(f'Integration path at {pl[1]:s} = {delta1*1e6:8.4g} um ')

    wake = wake1/delta1/charge  # V/C/m
    if len(fn) > 1:
        dados, info = _GdfidL_load_dados_info(_jnpth([path, fn[1]]))
        if pl[1] == 'x':
            delta2, _ = _GdfidL_get_integration_path(info)
        else:
            _, delta2 = _GdfidL_get_integration_path(info)
        _, wake2 = _np.loadtxt(dados, unpack=True)
        _log.info(f'and {delta2*1e6:8.4g} um')
        if pl[0] == 'd':
            wake = (wake1/delta1 - wake2/delta2)
            wake /= (1/delta1-1/delta2) * charge  # V/C
        else:
            wake = (wake1 - wake2)/(delta1-delta2)/charge  # V/C/m
    else:
        _log.info('')
    return wake


def _ECHO_rect_load_data(simul_data, code, path, anal_pl):
    def _load_dados(fname, mode, bc, code):
        if code == 'echozr':
            len_unit, charge_unit, header = 1e-2, 1e-12, 3  # cm to m, pC to C
            with open(fname) as f:
                f.readline()
                a = f.readline()
            mstep, offset, wid, bunlen = _np.fromstring(a[1:], sep='\t')
            offset = int(offset)
            # I don't know why I have to divide the echozr data by 2;
            arbitrary_factor = 2
            y0 = y = mstep*offset / 100
        elif code == 'echo2d':
            len_unit, charge_unit, header = 1, 1, 6
            with open(fname) as f:
                f.readline()
                mstep, offset = _np.fromstring(f.readline(), sep='\t')
                f.readline()
                wid, bunlen = _np.fromstring(f.readline(), sep='\t')
            offset = int(offset)
            # But I don't have to do this for the echo2d data.
            arbitrary_factor = 1
            y0 = y = mstep*offset
            offset = 0  # has only one column of wake
        spos, Wm = _np.loadtxt(
            fname, skiprows=header, usecols=(0, 1+offset), unpack=True)
        mstep *= len_unit
        wid *= len_unit
        bunlen *= len_unit
        spos *= len_unit
        # minus sign is due to convention:
        Wm *= -len_unit/charge_unit/arbitrary_factor

        Kxm = _np.pi/wid*mode
        if bc == 'elec':
            Wm /= _np.sinh(Kxm*y0)*_np.sinh(Kxm*y)
        else:
            Wm /= _np.cosh(Kxm*y0)*_np.cosh(Kxm*y)
        return spos, Wm, mstep, wid, bunlen

    if anal_pl == 'db':
        msg = 'Problem: All rectangular geometries does not have symmetry.'
        _log.info(msg)
        raise Exception(msg)

    bc = 'magn' if anal_pl == 'll' else 'elec'

    _log.info(f'Looking for data files in subfolder {bc:s}.')
    pname = _jnpth([path, bc])
    if not _os.path.isdir(pname):
        pname = path
        if code == 'echozr':
            _log.info(
                'Subfolder not found. It would be better to ' +
                'create the subfolder and put the files there...')
            _log.info('Looking for files in the current folder:')
        elif code == 'echo2d':
            msg = 'Files not found.'
            _log.info(msg)
            raise Exception(msg)

    f_in_dir = ' '.join(_os.listdir(pname))
    f_match = sorted(_re.findall(FNAME_ECHOZR2D, f_in_dir))
    if not f_match:
        msg = 'Files not found.'
        _log.info(msg)
        raise Exception(msg)

    _log.info(
        'Files found.\n I am assuming the simulation was performed ' +
        'with {0:s} boundary condition.'.format(
            'electric' if bc == 'elec' else 'magnetic'))
    _log.info('Modes found: ' + ', '.join([m for _, m in f_match]))
    _log.info('Loading data from files')

    spos, W, mode, mesh_size, width, bunlen = [], [], [], [], [], []
    for fn, m in f_match:
        if int(m) == 0:
            continue
        s, Wm, ms, wid, bl = _load_dados(_jnpth([pname, fn]), int(m), bc, code)
        mode.append(int(m))
        spos.append(s)
        W.append(Wm)
        mesh_size.append(ms)
        width.append(wid)
        bunlen.append(bl)

    cond = False
    for i in range(1, len(mode)):
        cond |= len(spos[i]) != len(spos[0])
        cond |= not _np.isclose(mesh_size[i], mesh_size[0], rtol=1e-5, atol=0)
        cond |= not _np.isclose(width[i], width[0], rtol=0, atol=1e-7)
        cond |= not _np.isclose(bunlen[i], bunlen[0], rtol=1e-5, atol=0)
        if cond:
            message = 'Parameters of file {0:s} differ from {1:s}.'.format(
                f_match[i][0], f_match[0][0])
            _log.info(message)
            raise Exception(message)

    simul_data.s = spos[0]
    simul_data.bunlen = bunlen[0]
    # I want the bunch to be symmetric:
    a = _np.argmin(_np.abs(spos[0] + spos[0][0])) + 1
    sbun = spos[0][:a]
    simul_data.sbun = sbun
    simul_data.bun = _np.exp(-sbun**2/(2*bunlen[0]**2))
    simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen[0]
    _log.info(f'Length of the driving bunch: {simul_data.bunlen*1e3:7.4g} mm')
    _log.info(f'Width of the simulated geometry: {width[0]*1e3:7.4g} mm')
    _log.info(f'Mesh step used in the simulation: {mesh_size[0]*1e6:7.4g} um')
    _log.info('All Data Loaded.')

    if anal_pl == 'll':
        _log.info('Calculating longitudinal Wake from data:')
        Wll, frac = None, 1
        for i in range(len(mode)):
            if mode[i] == 1:
                Wll = W[i].copy()
            elif mode[i] % 2:
                Wll += W[i]  # only odd terms
                frac = _np.max(_np.abs(W[i]/Wll))
        if Wll is None:
            _log.info('There is none odd mode to calculate Longitudinal Wake.')
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wll *= 2/width[0]
            simul_data.Wll = Wll

        _log.info('Calculating Quadrupolar Wake from data:')
        Wq, frac = None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wq = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wq += W[i] * Kxm**2  # only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wq))
        if Wq is None:
            _log.info('There is none odd mode to calculate Quadrupolar Wake.')
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wq *= 2/width[0]
            # minus sign is due to Panofsky-Wenzel
            Wq = -_scy_int.cumtrapz(Wq, x=spos[0], initial=0)
            simul_data.Wqy = Wq
            simul_data.Wqx = -Wq

        _log.info('Calculating Dipolar Horizontal Wake from data:')
        Wdx, frac = None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 2:
                Wdx = W[i].copy() * Kxm**2
            elif not mode[i] % 2:
                Wdx += W[i] * Kxm**2  # only even terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wdx))
        if Wdx is None:
            _log.info(
                "There's no even mode to calculate Dip. Horizontal Wake.")
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wdx *= 2/width[0]
            # minus sign is due to Panofsky-Wenzel
            Wdx = -_scy_int.cumtrapz(Wdx, x=spos[0], initial=0)
            simul_data.Wdx = Wdx

    elif anal_pl in {'dx', 'dy'}:
        pl = 'Vertical' if anal_pl == 'dy' else 'Horizontal'
        _log.info(f'Calculating Dipolar {pl:s} Wake from data:')
        Wd, frac = None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wd = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wd += W[i] * Kxm**2  # only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wd))
        if Wd is None:
            _log.info(
                f"There's no even mode to calculate Dipolar {pl:s} Wake.")
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wd *= 2/width[0]
            # minus sign is due to Panofsky-Wenzel
            Wd = -_scy_int.cumtrapz(Wd, x=spos[0], initial=0)
            setattr(simul_data, 'W'+anal_pl, Wd)
    else:
        msg = f'Plane of analysis {anal_pl:s} does not match any '
        msg += 'of the possible options'
        _log.info(msg)
        raise Exception(msg)
