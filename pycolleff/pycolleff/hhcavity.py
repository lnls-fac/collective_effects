#!/usr/bin/env python-sirius

import matplotlib.pyplot as plt
import matplotlib.gridspec as mpl_gs
from matplotlib import rcParams, rc

import numpy as np
import mathphys
from numpy.lib.shape_base import expand_dims
from pycolleff.simulate_landau import form_factor
import scipy.integrate as integrate

import pycolleff.impedances as imp
import pycolleff.sirius as si

_c = mathphys.constants.light_speed
_pi = np.pi


class Params:
    """."""

    def __init__(self):
        """."""
        self.U0 = None
        self.Vrf = None
        self.frf = None
        self.h = None
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


class HarmonicCavity:
    """."""

    def __init__(self, params=None):
        """."""
        self.params = params or Params()
        self._form_factor = 1
        self._wr = 0

    @property
    def k_harmonic(self):
        """."""
        nharm = self.params.nharm
        U0 = self.params.U0
        Vrf = self.params.Vrf
        kharm = 1/nharm**2 - ((U0/Vrf)**2)/(nharm**2-1)
        return kharm**(1/2)

    @property
    def phih_harmonic(self):
        """."""
        nharm = self.params.nharm
        U0 = self.params.U0
        Vrf = self.params.Vrf
        q = Vrf/U0
        numer = nharm/q
        denom = (nharm**2-1)**2 - (nharm**2/q)**2
        arg = -numer/(denom)**(1/2)
        return np.arctan(arg)/nharm

    @phih_harmonic.setter
    def phih_harmonic(self, value):
        self.psih_harmonic = _pi/2 - self.params.nharm * value

    @property
    def psih_harmonic(self):
        """."""
        return _pi/2 - self.params.nharm * self.phih_harmonic

    @psih_harmonic.setter
    def psih_harmonic(self, value):
        self.phih_harmonic = (_pi/2 - value)/self.params.nharm

    @property
    def perturbed_sync_phase(self):
        """."""
        nharm = self.params.nharm
        U0 = self.params.U0
        Vrf = self.params.Vrf
        phih = self.phih_harmonic
        kh = self.k_harmonic
        arg = U0/Vrf - kh*np.sin(nharm*phih)
        return _pi - np.arcsin(arg)

    @property
    def shunt_impedance(self):
        """."""
        kharm = self.k_harmonic
        psih = self.psih_harmonic
        Vrf = self.params.Vrf
        I0 = self.params.I0
        ffact = self.form_factor
        return kharm * Vrf / (2 * I0 * abs(ffact) * abs(np.cos(psih)))

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
        return np.arctan(Q * (wr/(nharm*wrf) - nharm*wrf/wr))

    @detune_angle.setter
    def detune_angle(self, value):
        Q = self.params.Q
        nharm = self.params.nharm
        wrf = 2*_pi*self.params.frf
        self.wr = nharm*wrf/(1+np.tan(value)/2/Q)

    # def wr_flat_potential(nharm, wrf, phih, quality):
    #     alpha = nharm * wrf * phih/2/quality
    #     wrp = alpha * (1 + np.sqrt(1+(nharm*wrf/alpha)**2))
    #     return wrp

    def detune_passive_cavity(self, Rs, form_factor=1):
        """."""
        kh = self.k_harmonic()
        I0 = self.params.I0
        Vrf = self.params.Vrf
        arg = kh*Vrf/(2*I0*abs(form_factor)*Rs)
        return -np.arccos(arg)

    def print_flat_potential(self):
        """."""
        st =  'k flat potential              : {:+08.3f} \n'
        st += 'harmonic phase [deg]          : {:+08.3f} \n'
        st += 'harmonic tuning angle [deg]   : {:+08.3f} \n'
        st += 'detuning frequency    [kHz]   : {:+08.3f} \n'
        st += 'shunt impedance flat [M.ohm]  : {:+08.4f} \n'
        st += 'unperturbed sync. phase [deg] : {:+08.3f} \n'
        st += 'perturbed sync. phase [deg]   : {:+08.3f} \n'
        rad2deg = 180/_pi

        kh = self.k_harmonic
        phih = self.phih_harmonic
        psih = self.psih_harmonic
        self.detune_angle = psih
        fr = self.wr/2/_pi
        df = self.params.nharm*self.params.frf - fr
        Rs_fp = self.shunt_impedance
        phis = self.params.sync_phase
        new_phis = self.perturbed_sync_phase
        print(
            st.format(
                kh, phih*rad2deg, psih*rad2deg, df*1e-3, Rs_fp*1e-6,
                phis*rad2deg, new_phis*rad2deg))

    def integrated_potential(self, z, harmonic=True):
        """."""
        wrf = 2*_pi*self.params.frf
        alpha = self.params.alpha
        sigmae = self.params.espread
        phis = self.params.sync_phase
        h = self.params.h
        tunes = self.params.tunes
        new_phis = self.perturbed_sync_phase
        nh = self.params.nharm
        phih = self.phih_harmonic
        kh = 0
        if harmonic:
            kh = self.k_harmonic
        phase = wrf*z/_c
        pot = (alpha**2 * sigmae**2)/(np.cos(phis)*(h*alpha*sigmae/tunes)**2)
        t1 = np.cos(new_phis)
        t2 = np.cos(phase + new_phis)
        t3 = kh/nh * (np.cos(nh*phih) - np.cos(nh*phase + nh*phih))
        t4 = (np.sin(new_phis) + kh*np.sin(nh*phih))*phase
        pot *= (t1 - t2 + t3 - t4)
        return pot

    def calc_distribution(self, z, harmonic=True):
        """."""
        pot = self.integrated_potential(z=z, harmonic=harmonic)
        alpha = self.params.alpha
        espread = self.params.espread
        dist = np.exp(-pot/(alpha**2 * espread**2))
        dist /= np.trapz(dist, z)
        return dist.ravel()

    @staticmethod
    def calc_sync_phase(z, dist):
        """."""
        return np.trapz(z*dist, z)

    @staticmethod
    def calc_bunch_length(z, dist):
        """."""
        zm = HarmonicCavity.calc_sync_phase(z, dist)
        z2 = np.trapz(z**2 * dist, z)
        return np.sqrt(z2 - zm**2)

    @staticmethod
    def complex_form_factor(z, w, rho):
        """."""
        return np.trapz(
            rho*np.exp(1j*w*z/_c), z)/np.trapz(rho, z)

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
        ring.dampte = np.inf
        ring.bunlen = self.params.bunlen
        if not radiation:
            ring.dampe = np.inf
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
        kh = self.k_harmonic
        nh = self.params.nharm
        phih = self.phih_harmonic

        phase = wrf*z/_c
        Vmain0 = Vrf * np.sin(phase + phis0)
        Vmain_pert = Vrf * np.sin(phase + phis_pert)
        Vharm = Vrf*kh*np.sin(nh*phase + nh*phih)
        return Vmain0, Vmain_pert, Vharm

    def calc_passive_voltage(self, z, Rs, detune_phase, form_factor=1):
        """."""
        I0 = self.params.I0
        wrf = 2*_pi*self.params.frf
        phase = wrf*z/_c
        nh = self.params.nharm
        volt = - 2*I0*Rs*form_factor
        volt *= np.cos(detune_phase)*np.cos(nh*phase-detune_phase)
        return volt

    def plot_voltages(self, z, vmain, vmain_pert, vharm):
        """."""
        vtotal = vmain_pert + vharm

        fig = plt.figure(figsize=(6, 4))
        gs = mpl_gs.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0, 0])

        ax1.plot(z*100, vmain * 1e-6, label=r'$V_{rf}$ Main')
        ax1.plot(z*100, vharm * 1e-6, label=r'$V_{3h}$ 3HC')
        ax1.plot(z*100, vtotal * 1e-6, label=r'$V_{rf} + V_{3h}$ Total')
        ax1.axhline(
            y=self.params.U0*1e-6, color='tab:gray', ls='--', label=r'$U_0$')

        ax1.set_xlabel(r'$z$ [cm]')
        ax1.set_ylabel('RF voltage [MV]')
        ax1.legend(loc='upper right')
        ax1.grid(ls='--', alpha=0.5)
        return fig
