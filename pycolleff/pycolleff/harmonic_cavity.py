"""."""
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
        # angle = _np.arctan2(tan, 1)
        angle = _np.arctan(tan)
        return angle

    @detune_angle.setter
    def detune_angle(self, value):
        Q = self.cavity_Q
        nharm = self.cavity_harm_num
        wrf = 2*_PI*self.rf_freq

        delta = _np.tan(value)/2/Q
        self.cavity_ang_freq = nharm*wrf*(delta + (1+delta**2)**(1/2))

        # alpha = _np.tan(value)/Q
        # self.cavity_ang_freq = nharm*wrf/2
        # self.cavity_ang_freq *= (alpha + (alpha**2 + 4)**(1/2))


class LongitudinalEquilibrium:
    """."""

    def __init__(
            self, ring=None,
            resonator=None, fillpattern=_np.array([])):
        """."""
        self.ring = ring or si.create_ring()
        self.harmcav = resonator
        self._zgrid = self.create_zgrid()
        self.fillpattern = fillpattern
        main_volt = self.voltage_main_cavity()
        self.dist = self.calc_distribution_from_voltage(main_volt)

    @property
    def zgrid(self):
        """."""
        return self._zgrid

    @zgrid.setter
    def zgrid(self, value):
        self._zgrid = value
        main_volt = self.voltage_main_cavity()
        self.dist = self.calc_distribution_from_voltage(main_volt)

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

    def harmonic_voltage_for_flat_potential(self):
        """."""
        nharm = self.harmcav.cavity_harm_num
        U0 = self.ring.en_lost_rad
        Vrf = self.ring.gap_voltage
        kharm = 1/nharm**2 - ((U0/Vrf)**2)/(nharm**2-1)
        return kharm**(1/2)

    def detune_angle_for_flat_potential(self):
        """."""
        Rs = self.harmcav.cavity_shunt_impedance
        I0 = _np.sum(self.fillpattern)
        Vrf = self.ring.gap_voltage
        k_fp = self.harmonic_voltage_for_flat_potential()
        ib = -2*I0*_np.sum(_np.abs(self.form_factor))
        arg = k_fp*Vrf/ib/Rs
        return -_np.arccos(arg)

    @property
    def form_factor(self):
        """."""
        wrf = 2*_PI*self.ring.rf_freq
        w = self.harmcav.cavity_harm_num*wrf
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
        voltage -= self.ring.en_lost_rad
        pot = -scy_int.cumtrapz(
            voltage, self.zgrid, initial=0)

        # subtract minimum value for all bunches
        pot -= _np.min(pot, axis=1)[:, None]

        const = self.ring.espread**2
        const *= self.ring.mom_comp
        const *= self.ring.circum
        # normalize by nominal energy
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
        if isinstance(w, (float)):
            # distribution: (h, zgrid)
            # dist*exp: (h, zgrid) * (None, zgrid) = (h, zgrid)
            # integration along zgrid, axis=1
            # fourier transform: (h, )
            arg = dist*_np.exp(1j*w*self.zgrid/_LSPEED)[None, :]
        else:
            # distribution: (h, zgrid)
            # dist*exp: (h, zgrid, None) * (zgrid, modes) = (h, zgrid, modes)
            # integration along zgrid, axis=1
            # fourier transform: (h, modes)
            arg = dist[:, :, None]*_np.exp(
                1j*w[None, :]*self.zgrid[:, None]/_LSPEED)
        return _np.trapz(arg, self.zgrid, axis=1)

    def voltage_main_cavity(self, sync_phase=None):
        """."""
        wrf = 2*_PI*self.ring.rf_freq
        phase = wrf*self.zgrid/_LSPEED
        phase0 = sync_phase or self.ring.sync_phase
        phase += phase0
        voltage = self.ring.gap_voltage*_np.sin(phase)
        return voltage

    def voltage_harmonic_cavity_impedance(self, dist=None, tol=1e-5):
        """."""
        if dist is None:
            dist = self.dist
        ind = _np.arange(self.ring.harm_num)
        ind_diff = ind[:, None] - ind[None, :]
        arg = -2j*_PI/self.ring.harm_num*ind_diff

        fill_fft = _np.fft.rfft(self.fillpattern)
        modes = _np.where(fill_fft > fill_fft[0]*tol)[0]
        harmonics = modes[1:]
        w0 = self.ring.rev_ang_freq
        pmain = int(self.harmcav.cavity_ang_freq/w0)
        # sum over modes (h*n+p) where p in [0, 1, -1, 2, -2, ...]
        ps = [0, ]
        ps += [tp for p in harmonics for tp in (p, -p)]
        ps = pmain + _np.array(ps)
        voltage = _np.zeros(
            (self.ring.harm_num, self.zgrid.size), dtype=_np.complex)
        zl_wps = self.harmcav.longitudinal_impedance(w=ps*w0).conj()
        print('     summing harmonics in voltage')
        for idx, p in enumerate(ps):
            wp = p * w0
            dist_fourier = self.distribution_fourier_transform(w=wp, dist=dist)
            zl_wp = zl_wps[idx] * _np.exp(1j*wp*self.zgrid/_LSPEED)
            bmat_wp = _np.exp(p*arg)
            beam_part = _np.dot(
                bmat_wp, self.fillpattern*dist_fourier.conj())
            # impedance has dim (zgrid, )
            # beam part has dim (h, )
            # voltage must have dim (h, zgrid)
            voltage += zl_wp[None, :] * beam_part[:, None]
            print(f'     {idx:02d}, {abs(p-pmain):03d}')
        # sum over positive frequencies only
        # ...
        voltage *= 2
        return voltage.real

    def calc_longitudinal_equilibrium(self, niter=100, tol=1e-10, beta=0.1):
        """."""
        main_volt = self.voltage_main_cavity()
        dists = [self.dist, ]
        for nit in range(niter):
            harmonic_volt = self.voltage_harmonic_cavity_impedance()
            # loss_cav = _np.mean(_np.max(harmonic_volt, axis=1))
            # ang = (self.ring.en_lost_rad - loss_cav)/self.ring.gap_voltage
            # sync_phase = _PI - _np.arcsin(ang)
            # main_volt = self.voltage_main_cavity(sync_phase)
            total_volt = main_volt[None, :] + harmonic_volt
            dist = self.calc_distribution_from_voltage(total_volt)
            dists.append(dist)
            self.dist = self.anderson_acceleration(dists, beta=beta)
            diff = (dists[-1]-dists[-2])/dists[-2]
            diff = _np.sqrt(_np.trapz(diff*diff, self.zgrid, axis=1))
            diff = _np.max(diff)
            print(nit+1, diff*(1+beta))
            if diff*(1+beta) < tol:
                print('distribution ok!')
                break
        return dists, harmonic_volt

    def anderson_acceleration(self, dists, beta):
        """."""
        return (dists[-1] + beta*dists[-2])/(1+beta)
