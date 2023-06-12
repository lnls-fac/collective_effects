"""."""
import time as _time

import numpy as _np
import numexpr as _ne

from mathphys import constants

_LSPEED = constants.light_speed


# Ring parameters
class Ring:
    """."""

    def __init__(self):
        """."""
        self._use_gaussian_noise = False
        self.energy = 3e9  # [eV]
        self.u0 = 871e3  # [eV]
        self.harm_num = 864
        self.rf_freq = 499663824  # [Hz]
        self.espread = 8.43589e-4
        self.mom_comp = 1.6446e-4
        self.damping_number = 1.6944268554007123

        self._cav_vgap = 0.0  # [V]
        self.cav_pos = _np.linspace(-1/2, 1/2, 1000)
        self.cav_pos *= self.rf_lamb
        self.cav_volt_norm = _np.zeros(self.cav_pos.shape)

        self.cav_vgap = 3e6  # [V]

    def to_dict(self):
        """."""
        return dict(
            use_gaussian_noise=self.use_gaussian_noise,
            energy=self.energy,
            u0=self.u0,
            harm_num=self.harm_num,
            rf_freq=self.rf_freq,
            espread=self.espread,
            mom_comp=self.mom_comp,
            damping_number=self.damping_number,
            cav_pos=self.cav_pos,
            cav_vgap=self._cav_vgap,
            cav_volt_norm=self.cav_volt_norm,
            )

    def from_dict(self, dic):
        """."""
        self.use_gaussian_noise = dic['use_gaussian_noise']
        self.energy = dic['energy']
        self.u0 = dic['u0']
        self.harm_num = dic['harm_num']
        self.rf_freq = dic['rf_freq']
        self.espread = dic['espread']
        self.mom_comp = dic['mom_comp']
        self.damping_number = dic['damping_number']
        self.cav_pos = dic['cav_pos']
        self._cav_vgap = dic['cav_vgap']
        self.cav_volt_norm = dic['cav_volt_norm']

    @property
    def use_gaussian_noise(self):
        """."""
        return self._use_gaussian_noise

    @use_gaussian_noise.setter
    def use_gaussian_noise(self, value):
        self._use_gaussian_noise = bool(value)

    @property
    def u0_norm(self):
        """."""
        return self.u0/self.energy

    @property
    def cav_vgap(self):
        """."""
        return self._cav_vgap

    @cav_vgap.setter
    def cav_vgap(self, value):
        self._cav_vgap = value

        rel_pos = self.cav_pos / self.rf_lamb
        self.cav_volt_norm = _np.sin(2*_np.pi * rel_pos + self.sync_phase)
        self.cav_volt_norm *= self.cav_vgap_norm
        self.cav_volt_norm -= self.u0_norm

    @property
    def sync_phase(self):
        """."""
        return _np.math.pi - _np.math.asin(self.u0/self.cav_vgap)

    @property
    def cav_vgap_norm(self):
        """."""
        return self.cav_vgap/self.energy

    @property
    def rf_lamb(self):
        """."""
        return _LSPEED/self.rf_freq

    @property
    def rev_freq(self):
        """."""
        return self.rf_freq/self.harm_num

    @property
    def rev_time(self):
        """."""
        return 1/self.rev_freq

    @property
    def circum(self):
        """."""
        return _LSPEED * self.rev_time

    @property
    def damping_time(self):
        """."""
        return 2*self.rev_time / self.damping_number / self.u0_norm

    @damping_time.setter
    def damping_time(self, value):
        """."""
        self.damping_number = 2*self.rev_time / value / self.u0_norm

    def track_one_turn(self, beam, excitation=True, damping=True):
        """."""
        damp = (1 - self.damping_number*self.u0_norm)
        if excitation:
            exc_quant = _np.sqrt(1 - damp*damp)*self.espread
            if not self._use_gaussian_noise:
                exc_quant *= _np.math.sqrt(12)
                noise = _np.random.rand(*beam.de.shape)
                noise -= 0.5
            else:
                noise = _np.random.randn(*beam.de.shape)  # 2/7 of exec. time
            noise *= exc_quant

        beam.ss += (self.circum * self.mom_comp) * beam.de
        if damping:
            beam.de *= damp
        if excitation:
            beam.de += noise
        beam.de += _np.interp(
            beam.ss, self.cav_pos, self.cav_volt_norm)  # 2/7 of exec. time


# Wake parameters:
class Wake:
    """."""

    def __init__(self, Q, Rs, wr):
        """."""
        self.Q = Q
        self.Rs = Rs  # [Ohm]
        self.wr = wr  # [rad/s]
        self.pot_phasor = 0j

    @property
    def Ql(self):
        """."""
        return _np.math.sqrt(self.Q*self.Q - 0.25)

    @property
    def kr(self):
        """."""
        return self.wr / _LSPEED

    @property
    def krl(self):
        """."""
        return self.kr * self.Ql / self.Q

    @property
    def cpl_kr(self):
        """."""
        return self.kr/(2*self.Q) + 1j*self.krl

    def cmd_reset_phasor(self):
        """."""
        self.pot_phasor = 0j

    def to_dict(self):
        """."""
        return dict(
            Q=self.Q,
            Rs=self.Rs,
            wr=self.wr,
            pot_phasor=self.pot_phasor,
            )

    def from_dict(self, dic):
        """."""
        self.Q = dic['Q']
        self.Rs = dic['Rs']
        self.wr = dic['wr']
        self.pot_phasor = dic['pot_phasor']

    def track_one_turn(self, beam, ring):
        """."""
        stren = ring.circum / _LSPEED / ring.energy
        stren *= self.wr * self.Rs / self.Q
        curr_p_part = beam.curr_p_bun / beam.num_part
        stren = stren*curr_p_part[:, None]

        rf_lamb = ring.rf_lamb
        cpl_kr = self.cpl_kr
        # First calculate the potential created by every particle
        # in every bunch:
        # expo = _np.exp(cpl_kr*beam.ss)  # 1/4 of the exec time
        ss = beam.ss
        expo = _ne.evaluate('exp(cpl_kr*ss)')  # much faster!
        pot = stren * expo

        # Now make a cumsum to get the total potential each particle
        # will be subjected to.
        pot_sum = _np.zeros((beam.num_buns, beam.num_part+1), dtype=complex)
        pot_sum[:, 1:] += _np.cumsum(pot, axis=1)

        # Now take the last term of pot_sum for each bunch, which is the total
        # potential created by that bunch, and cumsum it properly propagating
        # to the next bunch and also taking into account the potential of
        # previous turns:
        pot_summ = _np.zeros(beam.num_buns+1, dtype=complex)
        # propagate w_pot from the start of the ring to the first bunch:
        delta_bun = 0 - beam.bun_indices[0]
        pot_summ[0] = self.pot_phasor*_np.exp(cpl_kr*rf_lamb*delta_bun)

        for i, ind in enumerate(beam.bun_indices):
            # summ potential of this bunch with the one it generated:
            pot_summ[i+1] = pot_summ[i] + pot_sum[i, -1]
            # propagate the potential from this bunch to the next:
            if i == beam.bun_indices.size-1:
                ind_ip1 = ring.harm_num
            else:
                ind_ip1 = beam.bun_indices[i+1]
            dbun = ind - ind_ip1
            pot_summ[i+1] *= _np.exp(cpl_kr*rf_lamb*dbun)

        # Sum previous potentials to the one the particles will feel:
        pot_sumi = pot_sum[:, :-1]
        pot_sumi += pot_summ[:-1, None]

        # Apply the kick to the particles:
        kik = pot_sumi / expo
        beam.de -= 0.5*stren
        beam.de -= kik.real
        beam.de -= kik.imag/(2.0*self.Ql)

        # update the potential with the contributions of this turn:
        self.pot_phasor = pot_summ[-1]


class Beam():
    """."""

    def __init__(self, num_part, num_buns, current=0.35):
        """."""
        self.curr_p_bun = _np.ones(num_buns, dtype=float)
        self.bun_indices = _np.arange(num_buns)
        self.current = current
        self.de = _np.zeros((num_buns, num_part), dtype=float)
        self.ss = _np.zeros((num_buns, num_part), dtype=float)

    @property
    def current(self):
        """."""
        return float(self.curr_p_bun.sum())

    @current.setter
    def current(self, current):
        """."""
        self.curr_p_bun *= current/self.current

    @property
    def num_buns(self):
        """."""
        return int(self.bun_indices.size)

    @property
    def num_part(self):
        """."""
        return int(self.de.shape[1])

    def oversample_number_of_particles(
            self, mult_factor: int, noise_frac=0.0):
        """Increase number of particles by repeating the existing ones."""
        mult_factor = int(mult_factor)
        if mult_factor <= 1:
            return
        npart = self.num_part
        self.de = _np.tile(self.de, mult_factor)
        self.ss = _np.tile(self.ss, mult_factor)
        if not _np.math.isclose(noise_frac, 0):
            de_noise = self.de[:, :npart].std(axis=1) * noise_frac
            ss_noise = self.ss[:, :npart].std(axis=1) * noise_frac
            self.de[:, npart:] += de_noise[:, None] * _np.random.randn(
                self.num_buns, npart*(mult_factor-1))
            self.ss[:, npart:] += ss_noise[:, None] * _np.random.randn(
                self.num_buns, npart*(mult_factor-1))

    def to_dict(self):
        """."""
        return dict(
            ss=self.ss.copy(), de=self.de.copy(),
            curr_p_bun=self.curr_p_bun.copy(),
            bun_indices=self.bun_indices.copy())

    def from_dict(self, dic):
        """."""
        self.de = dic['de'].copy()
        self.ss = dic['ss'].copy()
        self.curr_p_bun = dic['curr_p_bun'].copy()
        self.bun_indices = dic['bun_indices'].copy()

    def sort(self):
        """."""
        idx = _np.argsort(self.ss, axis=1)  # 1/7 of the exec.time
        self.ss = _np.take_along_axis(self.ss, idx, axis=1)
        self.de = _np.take_along_axis(self.de, idx, axis=1)

    @staticmethod
    def calc_histogram(array, spos=None, nbins=50):
        """."""
        if spos is None:
            bins = _np.linspace(array.min(), array.max(), nbins+1)
        else:
            bins = _np.linspace(spos[0], spos[-1], nbins+1)
        bin_size = bins[1]-bins[0]

        indices = (array - bins[0])/bin_size
        indices = indices.astype(_np.intp)
        indices += _np.arange(array.shape[0])[:, None]*nbins

        hist = _np.bincount(indices.ravel(), minlength=nbins*array.shape[0])
        hist = hist.reshape(-1, nbins)

        if spos is not None:
            hist = _np.interp(spos, bins, hist)
        else:
            spos = (bins[1:] + bins[:-1])/2
        return spos, hist

    def generate_bunches(self, ring, ss_avg=0.0, de_avg=0.0):
        """."""
        # Energy distribution is Gaussian
        self.de = _np.random.randn(self.num_buns, self.num_part)
        self.de *= ring.espread
        self.de -= self.de.mean(axis=1)[:, None]
        self.de += de_avg

        # Draw longitudinal positions from equilibrium potential well:
        pot = -_np.cumsum(ring.cav_volt_norm)
        pot *= ring.cav_pos[1] - ring.cav_pos[0]
        pot -= pot.min()
        pot /= ring.mom_comp * ring.circum * ring.espread*ring.espread
        distr = _np.exp(-pot)
        int_distr = _np.cumsum(distr)
        int_distr *= ring.cav_pos[1] - ring.cav_pos[0]
        int_distr /= int_distr[-1]

        self.ss = _np.random.rand(self.num_buns, self.num_part)
        self.ss = _np.interp(self.ss, int_distr, ring.cav_pos)
        self.ss -= self.ss.mean(axis=1)[:, None]
        self.ss += ss_avg


def track_particles(
        ring, beam, wakes, num_turns=10000, stats_ev_nt=10,
        dist_ev_nt=1000, print_progress=True, save_dist=False,
        excitation=True, damping=True):
    """."""
    avg_ss, avg_de = [], []
    std_ss, std_de = [], []
    pot_wakes, tim = [], []

    t0 = _time.time()
    for turn in range(num_turns):
        ring.track_one_turn(beam, excitation=excitation, damping=damping)
        beam.sort()
        for wake in wakes:
            wake.track_one_turn(beam, ring)

        # save some statistics:
        if not (turn % stats_ev_nt):
            avg_ss.append(beam.ss.mean(axis=1))
            avg_de.append(beam.de.mean(axis=1))
            std_ss.append(beam.ss.std(axis=1))
            std_de.append(beam.de.std(axis=1))
            pot_wakes.append([wake.pot_phasor for wake in wakes])
            tim.append(ring.rev_time * turn)

        if save_dist:
            if not (turn % dist_ev_nt):
                _np.savez_compressed(
                    f'tracking_phase_space_turn_{turn:d}',
                    ss=beam.ss, de=beam.de)

        if print_progress and not (turn % 100):
            pot = 0.0
            if wakes:
                pot = _np.abs(wakes[0].pot_phasor)/ring.cav_vgap_norm
            print(
                f'{turn:06d}/{num_turns:06d}  ->  '
                f'<ss>(0) = {avg_ss[-1][0]*1e3:>8.4f} mm,   '
                f'STD(ss)(0) = {std_ss[-1][0]*1e3:>8.4f} mm,   '
                f'Pot = {pot:>7.4f}   '
                f'ET: {_time.time()-t0:.2f} s')
            t0 = _time.time()

    avg_ss.append(beam.ss.mean(axis=1))
    avg_de.append(beam.de.mean(axis=1))
    std_ss.append(beam.ss.std(axis=1))
    std_de.append(beam.de.std(axis=1))
    pot_wakes.append([wake.pot_phasor for wake in wakes])
    tim.append(ring.rev_time * turn)

    stats = {
        'avg_ss': _np.array(avg_ss), 'avg_de': _np.array(avg_de),
        'std_ss': _np.array(std_ss), 'std_de': _np.array(std_de),
        'pot_wakes': _np.array(pot_wakes), 'time': _np.array(tim)}
    return stats


def merge_stats(stats_list):
    """."""
    stats = dict()
    stats['avg_ss'] = _np.vstack([st['avg_ss'] for st in stats_list])
    stats['std_ss'] = _np.vstack([st['std_ss'] for st in stats_list])
    stats['avg_de'] = _np.vstack([st['avg_de'] for st in stats_list])
    stats['std_de'] = _np.vstack([st['std_de'] for st in stats_list])
    stats['pot_wakes'] = _np.vstack([st['pot_wakes'] for st in stats_list])
    off = 0
    tim = []
    for st in stats_list:
        tim.append(st['time'] + off)
        off += st['time'][-1]
    stats['time'] = _np.hstack(tim)
    return stats
