"""."""

import numpy as _np
import pandas as pd
import mathphys as _mp

_LSPEED = _mp.constants.light_speed
_factorial = _np.math.factorial
_sqrt = _np.math.sqrt


class Ring:
    """."""

    def __init__(self):
        """."""
        self.version = 'version'
        self.energy = 0.0  # energy [eV]
        self.rf_freq = 0.0  # RF frequency in [Hz]
        self.harm_num = 1    # harmonic Number
        self.mom_comp = 0.0  # momentum compaction factor
        self.tunex = 0.0  # horizontal tune
        self.tuney = 0.0  # vertical tune
        self.chromx = 0.0  # horizontal chromaticity
        self.chromy = 0.0  # vertical chromaticity
        self.num_bun = 1    # number of bunches filled

        self.total_current = 0.0  # total current [A]
        self.sync_tune = 0.0  # synchrotron tune
        self.espread = 0.0
        self.bunlen = 0.0  # bunch length in [m]
        self.damptx = 0.0  # horizontal damping time in [s]
        self.dampty = 0.0  # vertical damping time in [s]
        self.dampte = 0.0  # longitudinal damping time in [s]
        self.en_lost_rad = 0.0  # Energy lost per turn in [eV]
        self.gap_voltage = 0.0  # Gap Voltage in [V]

    @property
    def rf_lamb(self):
        """."""
        return _LSPEED/self.rf_freq

    @property
    def rev_freq(self):
        """."""
        return self.rf_freq/self.harm_num

    @property
    def rev_ang_freq(self):
        """."""
        return 2*_np.math.pi*self.rev_freq

    @property
    def rf_ang_freq(self):
        """."""
        return 2*_np.math.pi*self.rf_freq

    @property
    def rev_time(self):
        """."""
        return 1/self.rev_freq

    @property
    def circum(self):
        """."""
        return _LSPEED * self.rev_time

    @property
    def sync_phase(self):
        """Return the natural synchronous phase in rad.

        Returns:
            float: Natural synchronous phase [rad].

        """
        return _np.math.pi - _np.math.asin(self.en_lost_rad/self.gap_voltage)

    def __str__(self):
        """."""
        tmpl_st = '{0:28s}: {1:^20s}\n'
        tmpl_f = '{0:28s}: {1:^20.3f}\n'
        string = ''
        string += tmpl_st.format('Lattice Version', self.version)
        string += tmpl_f.format('Circumference [m]', self.circum)
        string += tmpl_f.format('Revolution Period [us]', self.rev_time*1e6)
        string += tmpl_f.format(
            'Revolution Frequency [kHz]', self.rev_freq/1e3)
        string += tmpl_f.format('Energy [GeV]', self.energy/1e9)
        string += tmpl_f.format('U0 [keV]', self.en_lost_rad/1e3)
        string += tmpl_f.format('Vgap [MV]', self.gap_voltage/1e6)
        string += '{0:28s}: {1:^20.2e}\n'.format(
            'Momentum Compaction', self.mom_comp)
        string += '{0:28s}: {1:^20d}\n'.format(
            'Harmonic Number', self.harm_num)
        string += tmpl_f.format('Current [mA]', self.total_current*1e3)
        string += tmpl_f.format(
            'Current per Bunch [mA]', self.total_current/self.num_bun*1e3)
        string += '{0:28s}: {1:^20.5f}\n'.format(
            'Synchrotron Tune', self.sync_tune)
        string += '{0:28s}: {1:>9.3f}/{2:<10.3f}\n'.format(
            'Tunes x/y', self.tunex, self.tuney)
        string += '{0:28s}: {1:>9.3f}/{2:<10.3f}\n'.format(
            'Chromaticities x/y', self.chromx, self.chromy)
        string += '{0:28s}: {1:>6.1f}/{2:^6.1f}/{3:<6.1f}\n'.format(
            'Damping Times x/y/e [ms]', self.damptx*1e3, self.dampty*1e3,
            self.dampte*1e3)
        string += '{0:28s}: {1:^20.4f}\n'.format(
            'Energy Spread [%]', self.espread*100)
        string += '{0:28s}: {1:^20.3f}\n'.format(
            'Bunch Length [mm]', self.bunlen*1e3)
        return string

    def to_dict(self):
        """Save state to dictionary."""
        return dict(
            version=self.version,
            energy=self.energy,
            rf_freq=self.rf_freq,
            harm_num=self.harm_num,
            mom_comp=self.mom_comp,
            tunex=self.tunex,
            tuney=self.tuney,
            chromx=self.chromx,
            chromy=self.chromy,
            num_bun=self.num_bun,
            total_current=self.total_current,
            sync_tune=self.sync_tune,
            espread=self.espread,
            bunlen=self.bunlen,
            damptx=self.damptx,
            dampty=self.dampty,
            dampte=self.dampte,
            en_lost_rad=self.en_lost_rad,
            gap_voltage=self.gap_voltage)

    def from_dict(self, dic):
        """Load state from dictionary."""
        self.version = dic.get('version', self.version)
        self.energy = dic.get('energy', self.energy)
        self.rf_freq = dic.get('rf_freq', self.rf_freq)
        self.harm_num = dic.get('harm_num', self.harm_num)
        self.mom_comp = dic.get('mom_comp', self.mom_comp)
        self.tunex = dic.get('tunex', self.tunex)
        self.tuney = dic.get('tuney', self.tuney)
        self.chromx = dic.get('chromx', self.chromx)
        self.chromy = dic.get('chromy', self.chromy)
        self.num_bun = dic.get('num_bun', self.num_bun)
        self.total_current = dic.get('total_current', self.total_current)
        self.sync_tune = dic.get('sync_tune', self.sync_tune)
        self.espread = dic.get('espread', self.espread)
        self.bunlen = dic.get('bunlen', self.bunlen)
        self.damptx = dic.get('damptx', self.damptx)
        self.dampty = dic.get('dampty', self.dampty)
        self.dampte = dic.get('dampte', self.dampte)
        self.en_lost_rad = dic.get('en_lost_rad', self.en_lost_rad)
        self.gap_voltage = dic.get('gap_voltage', self.gap_voltage)

    def get_voltage_waveform(self, zgrid, sync_phase=None):
        """."""
        wrf = 2*_np.pi*self.rf_freq
        phase0 = sync_phase or self.sync_phase
        phase = wrf * zgrid / _LSPEED
        phase += phase0
        voltage = self.gap_voltage*_np.sin(phase)
        return voltage

    def budget_summary(self, budget, fillpattern=None):
        """."""
        props = ['lsf', 'zln', 'pls', 'kdx', 'kdy', 'kqx', 'kqy',
                 'ktx', 'kty', 'ndx', 'ndy', 'nqx', 'nqy', 'ntx', 'nty']

        bud_res = pd.DataFrame(dict(
            name=pd.Series([
                'KLoss', 'Zl/n', 'PLoss',
                'Kdx', 'Kdy', 'Kqx', 'Kqx', 'Kx', 'Ky',
                'TuShdx', 'TuShdy', 'TuShqx', 'TuShqy', 'TuShx', 'TuShy'],
                index=props),
            unit=pd.Series([
                '[V/pC]', '[mOhm]', '[W]',
                '[kV/pC]', '[kV/pC]', '[kV/pC]',
                '[kV/pC]', '[kV/pC]', '[kV/pC]',
                '1/10^3', '1/10^3', '1/10^3',
                '1/10^3', '1/10^3', '1/10^3'],
                index=props),
            latex_unit=pd.Series([
                r'[V/pC]', r'[$m\Omega$]', r'[W]', r'[kV/pC]', r'[kV/pC]',
                r'[kV/pC]', r'[kV/pC]', r'[kV/pC]', r'[kV/pC]',
                r'$\times10^{-3}$',  r'$\times10^{-3}$', r'$\times10^{-3}$',
                r'$\times10^{-3}$', r'$\times10^{-3}$', r'$\times10^{-3}$'],
                index=props),
            latex_name=pd.Series([
                r'$\kappa_{Loss}$', r'$Z_L/n|_{eff}$', r'$P_{Loss}$',
                r'$\beta_x\kappa_x^D$', r'$\beta_y\kappa_y^D$',
                r'$\beta_x\kappa_x^Q$', r'$\beta_y\kappa_y^Q$',
                r'$\beta_x\kappa_x$', r'$\beta_y\kappa_y$',
                r'$\Delta\nu_x^D$', r'$\Delta\nu_y^D$', r'$\Delta\nu_x^Q$',
                r'$\Delta\nu_y^Q$', r'$\Delta\nu_x$', r'$\Delta\nu_y$'],
                index=props)
            ))

        convert = dict(
            lsf=1e-12, zln=1e3, pls=1, kdx=1e-15, kdy=1e-15, kqx=1e-15,
            kqy=1e-15, ktx=1e-15, kty=1e-15, ndx=1e3, ndy=1e3, nqx=1e3,
            nqy=1e3, ntx=1e3, nty=1e3)

        for el in budget:
            values = dict()
            w = el.ang_freq

            Zl = el.Zll * el.quantity
            if len(Zl) != 0:
                lossf, pl, zn, *_ = self.loss_factor(
                    w=w, Zl=Zl, fillpattern=fillpattern)
                values['lsf'] = lossf
                values['pls'] = pl
                values['zln'] = zn

            Zd = el.Zdy * el.quantity * el.betay
            if len(Zd) != 0:
                kd, tus = self.kick_factor(
                    w=w, Z=Zd, imp_type='Zdy', fillpattern=fillpattern)
                values['kdy'] = kd
                values['ndy'] = tus

            Zq = el.Zqy * el.quantity * el.betay
            if len(Zq) != 0:
                kd, tus = self.kick_factor(
                    w=w, Z=Zq, imp_type='Zqy', fillpattern=fillpattern)
                values['kqy'] = kd
                values['nqy'] = tus

            if len(Zd) != 0 or len(Zq) != 0:
                values['kty'] = values.get('kdy', 0) + values.get('kqy', 0)
                values['nty'] = values.get('ndy', 0) + values.get('nqy', 0)

            Zd = el.Zdx * el.quantity * el.betax
            if len(Zd) != 0:
                kd, tus = self.kick_factor(
                    w=w, Z=Zd, imp_type='Zdx', fillpattern=fillpattern)
                values['kdx'] = kd
                values['ndx'] = tus

            Zq = el.Zqx * el.quantity * el.betax
            if len(Zq) != 0:
                kd, tus = self.kick_factor(
                    w=w, Z=Zq, imp_type='Zqx', fillpattern=fillpattern)
                values['kqx'] = kd
                values['nqx'] = tus

            if len(Zd) != 0 or len(Zq) != 0:
                values['ktx'] = values.get('kdx', 0) + values.get('kqx', 0)
                values['ntx'] = values.get('ndx', 0) + values.get('nqx', 0)

            for prop in values.keys():
                values[prop] *= convert[prop]

            bud_res[el.name] = pd.Series(values)
        return bud_res

    def loss_factor(
            self, budget=None, element=None, w=None, Zl=None,
            fillpattern=None):
        """Calculate the loss factor and effective impedance.

        The Effective impedance returned is the one of the azimuthal mode m=1
        and radial mode k=0. The implementation follows closely eqs. 6.140
        and 6.143 of ref. [1] and eq. 9.94 of ref. [2]. However, we do keep
        here the more general scenario of ref. [2] of a uniform filling and
        calculate the effective impedance for coupled bunch mode cbmode=0.

        References:
            [1] Chao, A. W. (1993). Physics of Collective Beam Instabilities
                in High Energy Accelerators (1st ed.). New York: John Wiley &
                Sons.
            [2] Ng, K.-Y. (2006). Physics of Intensity Dependent Beam
                Instabilities (1st ed.). https://doi.org/10.1142/9789812703392

        Attributes and properties used by this method:
            total_current
            rev_ang_freq
            rev_time
            bunlen
            num_bun

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            Zl (numpy.ndarray, (M, ), optional): longitudinal impedance in
                units of [Ohm] evaluated at angular frequencies in `w`.
                Only considered if `budget` and `element` are `None`.
                Defaults to None.
            fillpattern (numpy.ndarray, (h, ), optional): vector with the size
                of the harmonic number describing the current profile. All
                entries must be non-negative and sum to one. If not None it
                will overwrite the parameter self.num_bun in calculations.

        Returns:
            lossf (float): Loss factor in [V/C].
            Pl (float): Power loss in [W].
            Zovn_eff (float): Effective Impedance (Z/n) in [Ohm].
            wp (numpy.ndarray, (N, )): vector of angular frequencies where the
                impedance was sampled [rad/s].
            lossfp (numpy.ndarray, (N, )): vector of loss factor with
                contribution of each sampled frequency [V/C].

        """
        w0 = self.rev_ang_freq
        nb = self.num_bun
        bunlen = self.bunlen
        tot_curr = self.total_current
        rev_time = self.rev_time

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        limw = 5 * _LSPEED / bunlen
        wmin = max(w[0], -limw)
        wmax = min(w[-1], limw)
        wp_args = wmin, wmax, w0
        wp_kws = dict()
        nb, wp, dft_sqr = self._process_fillpattern(
            fillpattern, nb, wp_args, wp_kws)

        gmk = self.calc_spectrum(wp, bunlen, max_rad=0, max_azi=1)
        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)

        # Multibunch loss factor and power loss:
        h00 = gmk[(0, 0)]**2
        lossfp = Zl_interp.real * h00
        if fillpattern is not None:
            lossfp *= dft_sqr
        lossf = nb * w0/(2*_np.pi) * lossfp.sum()
        Pl = lossf * rev_time * tot_curr**2 / nb
        if fillpattern is not None:
            lossf /= _np.sum(fillpattern*fillpattern)

        # eq. 9.94 of ref. [2] and eqs. 6.143 + 6.140 of ref. [1]
        # please, notice the factor of 2 here to convert from gmk to hmk:
        # Also note that we are interested in Z/n_eff, not Z_eff, that's the
        # reason for us to multiply by w0.
        h10 = gmk[(1, 0)]**2 * 2
        Zovn_eff = w0*_np.sum(Zl_interp.imag*h10/(wp+1e-4)) / h10.sum()

        return lossf, Pl, Zovn_eff, wp, lossfp

    def kick_factor(
            self, budget=None, element=None, w=None, Z=None, imp_type='Zdy',
            fillpattern=None):
        """Calculate the kick factor, tune shift and effective impedance.

        The kick factor and tune-shifts will be calculated for uniform
        filling, zero chromaticity and the azimuthal mode m=0.

        Attributes and properties used by this method:
            total_current
            rev_ang_freq
            rev_time
            bunlen
            num_bun
            tunex
            tuney

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            Z (numpy.ndarray (M, ), optional): longitudinal impedance in
                units of [Ohm] evaluated at angular frequencies in `w`.
                Only considered if `budget` and `element` are `None`.
                Defaults to None.
            imp_type (str, optional): string with type of impedance to
                consider. Options: Zdy, Zqy, Zdx, Zqx. Defaults to Zdy.
            fillpattern (numpy.ndarray, (h, ), optional): vector with the size
                of the harmonic number describing the current profile. All
                entries must be non-negative and sum to one. If not None it
                will overwrite the parameter self.num_bun in calculations.

        Returns:
            kick_factor (float): Kick Factor factor in [V/C]
            tune_shift (float): Tune shift of azimuthal mode m=0.
            Zt_eff (complex): Effective Impedance in [Ohm]

        """
        if imp_type.lower().startswith(('zqx', 'zqy')):
            nut = 0
        elif imp_type.lower().startswith('zdx'):
            nut = self.tunex % 1
        else:
            nut = self.tuney % 1

        energy = self.energy
        w0 = self.rev_ang_freq
        nb = self.num_bun
        bunlen = self.bunlen
        tot_curr = self.total_current

        w, Z = self._prepare_input_impedance(budget, element, w, Z, imp_type)

        limw = 5 * _LSPEED / bunlen
        wmin = max(w[0], -limw)
        wmax = min(w[-1], limw)
        wp_args = wmin, wmax, w0
        wp_kws = dict(nut=nut)
        nb, wp, dft_sqr = self._process_fillpattern(
            fillpattern, nb, wp_args, wp_kws)

        gmk = self.calc_spectrum(wp, bunlen, max_rad=0, max_azi=1)
        Zt_interp = self._get_interpolated_impedance(wp, w, Z)

        h00 = gmk[(0, 0)]**2
        kick_factorp = Zt_interp.imag*h00
        if fillpattern is not None:
            kick_factorp *= dft_sqr
        kick_factor = nb * w0/(2*_np.pi) * kick_factorp.sum()
        Tush = (tot_curr/nb) / (2*energy) * kick_factor / w0
        if fillpattern is not None:
            kick_factor /= _np.sum(fillpattern*fillpattern)

        return kick_factor, Tush

    def _process_fillpattern(self, fillp, nb, wp_args, wp_kws):
        hnum = self.harm_num
        dft_sqr = None
        if fillp is not None:
            if fillp.ndim != 1:
                raise ValueError('fillpattern must be one-dimensional.')
            elif fillp.size != hnum:
                raise ValueError('Length of fillpattern must be harm_num.')
            elif not _np.allclose(fillp.sum(), 1):
                raise ValueError('fillpattern must sum to unit.')
            elif _np.any(fillp < 0):
                raise ValueError('All entries of fillpattern must be >= 0.')
            nb = 1
            wp, p = self._get_sampling_ang_freq(
                *wp_args, nb, return_p=True, **wp_kws)
            pmin, pmax = abs(int(p[0])), abs(int(p[-1]))+1
            tile_neg, tile_pos = int(pmin/hnum)+1, int(pmax/hnum)+1

            dft_sqr = _np.abs(_np.fft.fft(fillp))**2
            dft_sqr_neg = _np.tile(dft_sqr, tile_neg)[::-1][:pmin][::-1]
            dft_sqr_pos = _np.tile(dft_sqr, tile_pos)[:pmax]
            dft_sqr = _np.r_[dft_sqr_neg, dft_sqr_pos]
        else:
            wp = self._get_sampling_ang_freq(*wp_args, nb, **wp_kws)
        return nb, wp, dft_sqr

    def longitudinal_cbi(
            self, budget=None, element=None, w=None, Zl=None, m=1,
            inverse=False, full=False):
        """Calculate coupled-bunch relative (and complex) tune-shifts.

        The coupled-bunch formulas implemented here are a special case of the
        ones implemented in method `longitudinal_mode_coupling`. They are
        valid in the limit of small currents per bunch, where the mode-coupling
        is negligible, and each mode (m, k) evolves independently of other
        modes in time. The formulas implemented here are valid for radial mode
        k=0. They can be found in eq. 6.129 of ref. [1], for example, when we
        take k = k' = 0.

        To get the real tune-shifts you must take the real part of the
        eigen-values and multiply them by the synchrotron tune, `sync_tune`.
        To get the growth rate in [1/s] you must multiply the imaginary part
        by the angular synchrotron frequency `sync_tune*w0`. This
        implementation also takes into account the radiation damping, so that
        the instability occurs when any of the growth rates are larger than 0.

        References:
            [1] Chao, A. W. (1993). Physics of Collective Beam Instabilities
                in High Energy Accelerators (1st ed.). New York: John Wiley &
                Sons.
            [2] Ng, K.-Y. (2006). Physics of Intensity Dependent Beam
                Instabilities (1st ed.). https://doi.org/10.1142/9789812703392

        Attributes and properties used by this method:
            sync_tune
            rev_ang_freq
            mom_comp
            bunlen
            energy
            num_bun
            total_current

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            Zl (numpy.ndarray (M, ), optional): longitudinal impedance in
                units of [Ohm] evaluated at angular frequencies in `w`.
                Only considered if `budget` and `element` are `None`.
                Defaults to None.
            m (int, optional): azimuthal mode to consider. Defaults to 1.
            inverse (bool, optional): If this flag is True, the miminum
                impedance as function of angular frequency w needed to create
                instability is returned, instead of the tune-shifts. If True,
                budget, element and Zl arguments are ignored.
                Defaults to False.
            full (bool, optional): Whether or not to return wp and Zl_wp (see
                Returns section below.) togheter with the tune-shifts.
                Defaults to False.

        Returns:
            tune_shift (numpy.ndarray, (nbun, )): complex tune-shift for each
                coupled bunch mode.
            wp (numpy.ndarray, (nbun, N)): angular frequencies used in the
                summation for each coupled-bunch mode.
            Zl_wp (numpy.ndarray, (nbun, N)): impedance sampled at wp.

        Raises:
            ValueError: when azimuthal mode m <= 0.

        """
        if abs(m) <= 0:
            raise ValueError('azimuthal mode m must be greater than zero.')
        sync_tune = self.sync_tune
        w0 = self.rev_ang_freq
        ws = sync_tune * w0
        eta = self.mom_comp
        bunlen = self.bunlen
        energy = self.energy
        nb = self.num_bun
        tot_curr = self.total_current
        curr_p_bun = tot_curr / nb
        alpe = 1/self.dampte

        # Look at eq. 6.157 of ref. [1] and the discussion that follows to
        # understand the motivation behind the definition of xi:
        xi = eta * curr_p_bun / sync_tune**2 / energy
        sig_theta = w0*bunlen/_LSPEED
        K = nb * xi * w0/(2*_np.pi) / (sig_theta)**2

        if inverse:
            hm0 = self.calc_spectrum(
                w, bunlen, max_rad=0, max_azi=m, only=True)
            hm0 *= hm0
            Z_min = alpe/ws / K * w/hm0
            return Z_min

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        # Calculate Effective Impedance
        cbmodes = _np.arange(nb)
        wp = self._get_sampling_ang_freq(
            w[0], w[-1], w0, nb, m, sync_tune, cbmodes)
        gm0 = self.calc_spectrum(wp, bunlen, max_rad=0, max_azi=m, only=True)
        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)
        Zl_sum = (Zl_interp/wp * gm0 * gm0).sum(axis=1)

        # Returns the relative Tune Shift:
        tune_shift = 1j*(K*m*Zl_sum - abs(m)*alpe/ws)

        if full:
            return tune_shift, wp, Zl_interp, gm0
        return tune_shift

    def transverse_cbi(
            self, budget=None, element=None, w=None, beta_Zt=None, m=0,
            plane='y', inverse=False, full=False):
        """Calculate coupled-bunch relative (and complex) tune-shifts.

        The coupled-bunch formulas implemented here are a special case of the
        ones implemented in method `transverse_mode_coupling`. They are valid
        in the limit of small currents per bunch, where the mode-coupling is
        negligible, and each mode (m, k) evolves independently of other modes
        in time. The formulas implemented here are valid for radial mode k=0.
        They can be found in eq. 6.198 of ref. [1] and in eq. 9.38 of ref. [2],
        for example, when we take k = k' = 0.

        To get the real tune-shifts you must take the real part of the
        eigen-values and multiply by the synchrotron tune, `sync_tune`. To get
        the growth rate  in [1/s] you must multiply the imaginary part by the
        angular synchrotron frequency `sync_tune*w0`. This implementation also
        takes into account the radiation damping, so that the instability
        occurs when any of the growth rates are larger than 0.

        References:
            [1] Chao, A. W. (1993). Physics of Collective Beam Instabilities
                in High Energy Accelerators (1st ed.). New York: John Wiley &
                Sons.
            [2] Ng, K.-Y. (2006). Physics of Intensity Dependent Beam
                Instabilities (1st ed.). https://doi.org/10.1142/9789812703392

        Attributes and properties used by this method:
            sync_tune
            rev_ang_freq
            mom_comp
            energy
            num_bun
            total_current
            dampte
            bunlen
            damptx, tunex, chromx
            dampty, tuney, chromy

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            beta_Zt (numpy.ndarray (M, ), optional): beta weighted transverse
                dipolar impedance in units of [Ohm] evaluated at angular
                frequencies in `w`. Only considered if `budget` and `element`
                are `None`. Defaults to None.
            m (int, optional): azimuthal mode to consider. Defaults to 1.
            plane (str, optional): transverse plane where to make the analysis.
                Can assume {'x', 'y'}. Defaults to 'y'.
            inverse (bool, optional): If this flag is True, the miminum
                impedance as function of angular frequency w needed to create
                instability is returned, instead of the tune-shifts. If True,
                budget, element and Zl arguments are ignored.
                Defaults to False.
            full (bool, optional): Whether or not to return wp and Zl_wp (see
                Returns section below.) togheter with the tune-shifts.
                Defaults to False.

        Returns:
            tune_shift (numpy.ndarray, (nbun, )): complex normalized
                tune-shift for each coupled bunch mode.
            wp (numpy.ndarray, (nbun, N)): angular frequencies used in the
                summation for each coupled-bunch mode.
            Zl_wp (numpy.ndarray, (nbun, N)): impedance sampled at wp.

        """
        sync_tune = self.sync_tune
        w0 = self.rev_ang_freq
        ws = sync_tune * w0
        eta = self.mom_comp
        energy = self.energy
        nb = self.num_bun
        tot_curr = self.total_current
        taue = self.dampte
        bunlen = self.bunlen
        if plane.lower().startswith(('x', 'h')):
            taut, nut, chrom, imp = self.damptx, self.tunex, self.chromx, 'Zdx'
        else:
            taut, nut, chrom, imp = self.dampty, self.tuney, self.chromy, 'Zdy'

        alpe = 1/taue
        alpt = 1/taut
        w_crom = chrom/eta*w0  # chrom. frequency shift of the bunch spectrum

        # Notice this definition has an extra 1/2/pi in relation to eq. 9.38
        # of ref. [2]. However, notice that equation 9.64 defines labmda_mk
        # with an extra 1/sqrt(2*pi) in relation to our g_mk.
        K = tot_curr*w0/(4*_np.pi)/ws/energy

        if inverse:
            gm0 = self.calc_spectrum(
                w - w_crom, bunlen, max_rad=0, max_azi=m, only=True)
            gm0 *= gm0
            Z = (alpt + abs(m)*alpe)/ws / K / gm0
            return Z

        w, beta_Zt = self._prepare_input_impedance(
            budget, element, w, beta_Zt, imp)

        # Calculate Effective Impedance
        cbmodes = _np.arange(nb)
        wp = self._get_sampling_ang_freq(
            w[0], w[-1], w0, nb, m, sync_tune, cbmodes, nut)
        wpcro = wp - w_crom

        gm0 = self.calc_spectrum(
            wpcro, bunlen, max_rad=0, max_azi=m, only=True)
        gm0 *= gm0
        Zt_interp = self._get_interpolated_impedance(wp, w, beta_Zt)
        Zt_sum = (Zt_interp * gm0).sum(axis=1)

        # Calculate Coupled_bunch Instability
        # Returns the relative Tune Shift:
        tune_shift = -1j*(alpt/ws + abs(m)*alpe/ws + K * Zt_sum)

        if full:
            return tune_shift, wp, Zt_interp, gm0
        return tune_shift

    def longitudinal_mode_coupling(
        self,
        budget=None,
        element=None,
        w=None,
        Zl=None,
        max_azi=10,
        max_rad=12,
        cbmode=0,
        use_fokker=True,
        modecoup_matrix=None,
        fokker_matrix=None,
        delete_m0=True,
        delete_m0k0=True
    ):
        """Calculate the longitudinal mode-coupling eigen-values.

        The eigen-values returned here are normalized by the synchrotron
        frequency, which means they are adimensional. Besides, the
        implementation does not guarantee any specific ordering.

        To get the tune-shifts you must take the real part of the eigen-values.
        To get the growth rate in [1/s] you must multiply the imaginary part
        by the angular synchrotron frequency `sync_tune*w0`.
        Instability occurs when any of the growth rates are larger than zero.

        We follow ref.[1] closely, such that the comments in the
        implementation refer to equation of this reference. Ref. [2] gives
        identical results of ref. [1] and ref.[3] doesn't give explicit
        formulas for gaussian beams, so we leave it here only for the
        interested user.

        We make a small generalization of ref. [1], by considering that the
        ring is filled with several equally spaced identical bunches and that
        the eigen-values are being calculated for a given coupled-bunch mode
        of the beam.

        References:
            [1] Suzuki, T., Chin, Y.-H., & Satoh, K. (1983). Mode-coupling
                theory and bunch lengthening in spear. Particle Accelerators,
                13, 179-198.
            [2] Suzuki, T. (1983). Theory of longitudinal bunched-beam
                instabilities based on the fokker-planck equation. Particle
                Accelerators, 14, 91-108.
            [3] Cai, Y. (2011). A Linear Theory of Microwave Instability in
                Electron Storage Rings. Physical Review Special Topics -
                Accelerators and Beams, 14(6), 061002.
                https://doi.org/10.1103/PhysRevSTAB.14.061002

        Attributes and properties used by this method:
            total_current
            energy
            rev_ang_freq
            sync_tune
            bunlen
            mom_comp
            num_bun
            dampte

        The size N of each dimension of the mode-coupling and fokker plank
        matrices is given by:
            N = (2*max_azi+1)*(max_rad+1))

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            Zl (numpy.ndarray (M, ), optional): longitudinal impedance in
                units of [Ohm] evaluated at angular frequencies in `w`.
                Only considered if `budget` and `element` are `None`.
                Defaults to None.
            max_azi (int, optional): Maximum azimuthal mode to consider in the
                truncation. The solution will consider all `2*max_azi + 1`
                azimuthal modes whose absolute value are <= than `max_azi`.
                Defaults to 3.
            max_rad (int, optional): Maximum radial mode to consider in the
                truncation. The solution will consider all `max_rad + 1`
                radial modes <= than `max_rad`. Defaults to 4.
            cbmode (int, optional): Coupled-bunch mode for which the
                eigen-values will be found. Must be lower than `num_bun`.
                Defaults to 0.
            use_fokker (bool, optional): Whether or not to include the
                mode-coupling fokker-planck terms described in ref. [2].
                If False, only the diagonal damping terms will be included.
                Defaults to True.
            modecoup_matrix (numpy.ndarray, (N, N), optional): the
                mode-coupling matrix to be used in calculations. If None, then
                it will be calculated internally. Defaults to None.
            fokker_matrix (numpy.ndarray, (N, N), optional): the fokker planck
                matrix to be used in calculations. If None, then it will be
                calculated internally. Defaults to None.
            delete_m0 (bool, optional): Whether or not to remove mode m=0 from
                calculations. Defaults to True.
            delete_m0k0 (bool, optional): Whether or not to remove mode m=0
                and k=0 from calculations. Only relevant when delete_m0 is
                False. Defaults to True.

        Returns:
            eigenvals (numpy.ndarray, (N, ): normalized eigen-modes of the
                mode-coupling problem.
            modecoup_matrix (numpy.ndarray, (N, N)): the mode-coupling matrix
                used in calculations.
            fokker_matrix (numpy.ndarray, (N, N)): the fokker planck matrix
                used in calculations.

        """
        energy = self.energy
        w0 = self.rev_ang_freq
        bunlen = self.bunlen
        sync_tune = self.sync_tune
        eta = self.mom_comp
        nb = self.num_bun
        I_b = self.total_current / nb
        alpe = 1/self.dampte
        if cbmode >= nb:
            cbmode = 0
            print(
                'Coupled bunch mode greater than number of bunchs.\n',
                'Reseting cbmode to 0.')

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        # There is an approximation here which is only valid for broad-band
        # impedances. I sample the impedance only at the azimuthal mode m=1.
        wp = self._get_sampling_ang_freq(
            w[0], w[-1], w0, nb, 1, sync_tune, [cbmode])

        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)
        Zl_wp = Zl_interp / wp

        # Calculate the fokker-planck matrix:
        if fokker_matrix is None:
            fokker_matrix = self._calc_fokker_planck(
                max_azi, max_rad, alpe, use_fokker)
        # Calculate the mode coupling matrix:
        if modecoup_matrix is None:
            modecoup_matrix = self._calc_vlasov(
                Zl_wp, wp, bunlen, max_azi, max_rad)
        # Calculate the current independent diagonal matrix:
        ms = _np.arange(-max_azi, max_azi+1)
        D = _np.einsum('mn,kl->mknl', _np.diag(ms), _np.eye(max_rad+1))
        D = self._reshape_coupling_matrix(D)

        # We separated the K value from the defition of M so that M
        # could be current independent in case of a constant bunch
        # length. To understand the equation for K implemented here,
        # please look at eqs. 41 and 43 of ref. [2]. and notice that
        # the definition of M_mlnh in 2.26 of ref. [1] has Z/p instead
        # of Z/wp as we have here:
        xi = eta * I_b/sync_tune**2/energy
        sig_theta = w0*bunlen/_LSPEED
        K = xi * w0/(2*_np.pi)/(sig_theta)**2
        # Lastly we need to multiply by `nb` because of our
        # generalization to include n equally spaced bunches.
        K *= nb

        # Please, check eq. 43 of ref. [2]:
        A = D + 1j*K*modecoup_matrix + 1j*fokker_matrix / (sync_tune*w0)

        if delete_m0:
            nr_azi = 2*max_azi + 1
            idx_m0 = max_azi+0 + nr_azi*_np.arange(max_rad+1)
            A = _np.delete(_np.delete(A, idx_m0, axis=0), idx_m0, axis=1)
        elif delete_m0k0:
            nr_azi = 2*max_azi + 1
            idx_m0k0 = max_azi+0 + nr_azi*0
            A = _np.delete(_np.delete(A, idx_m0k0, axis=0), idx_m0k0, axis=1)

        return _np.linalg.eigvals(A), modecoup_matrix, fokker_matrix

    @classmethod
    def _calc_vlasov(cls, Z_wp, wp, bunlen, max_azi, max_rad):
        """Calculate the longitudinal mode-coupling matrix.

        We follow ref.[1] closely, such that the comments in the
        implementation refer to equation of this reference. Ref. [2] gives
        identical results of ref. [1] and ref.[3] doesn't give explicit
        formulas for gaussian beams, so we leave it here only for the
        interested user.

        The mode-coupling matrix is actually an object with 4 indices: Mmknl,
        because it relates the expansion coefficients of the oscillation modes,
        which are 2-indices objects, a_mk, being m for the azimuthal modes
        and k for the radial modes.

        However, after filling this 4-indices "matrix" we cast it in a
        regular 2D matrix by mapping the m and k coefficients into a single
        dimension, so we can solve the eigen-value problem with standard
        techniques.

        References:
            [1] Suzuki, T., Chin, Y.-H., & Satoh, K. (1983). Mode-coupling
                theory and bunch lengthening in spear. Particle Accelerators,
                13, 179-198.
            [2] Suzuki, T. (1983). Theory of longitudinal bunched-beam
                instabilities based on the fokker-planck equation. Particle
                Accelerators, 14, 91-108.
            [3] Cai, Y. (2011). A Linear Theory of Microwave Instability in
                Electron Storage Rings. Physical Review Special Topics -
                Accelerators and Beams, 14(6), 061002.
                https://doi.org/10.1103/PhysRevSTAB.14.061002

        Args:
            Z_wp (numpy.ndarray, (N, )): Impedance divided sampled at wp
                angular frequency and divided by wp. Unit in [Ohm/(rad/s)]
            wp (numpy.ndarray, (N, )): angular frequency in units of [rad/s]
            bunlen (float): bunch length in units of [m]
            max_azi (int): maximum azimuthal mode to consider.
            max_rad (int): maximum radial mode to consider.

        """
        def fill_symmetric_terms(m, k, n, l, Mmknl):
            # Please, note the symmetry relations in eq. 2.46 and also notice
            # by 2.26 that Mmknl = Mnlmk * (-1)**(m-n):
            M[max_azi+m, k, max_azi+n, l] = m*Mmknl
            M[max_azi-m, k, max_azi+n, l] = -m*Mmknl
            M[max_azi+m, k, max_azi-n, l] = m*Mmknl
            M[max_azi-m, k, max_azi-n, l] = -m*Mmknl
            Mnlmk = Mmknl * (-1)**(m-n)
            M[max_azi+n, l, max_azi+m, k] = n*Mnlmk
            M[max_azi+n, l, max_azi-m, k] = n*Mnlmk
            M[max_azi-n, l, max_azi+m, k] = -n*Mnlmk
            M[max_azi-n, l, max_azi-m, k] = -n*Mnlmk

        M = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        spectrum = cls.calc_spectrum(wp, bunlen, max_rad, max_azi)
        for m in range(max_azi+1):
            for k in range(max_rad+1):
                Imk = spectrum[(abs(m), k)]
                # For l=k it is only necessary to cover n=m..max_azi
                # because of the symmetry relations of eq. 2.46:
                for n in range(m, max_azi+1):
                    Ink = spectrum[(abs(n), k)]
                    Mmknk = (1j)**(m-n)*_np.dot(Z_wp, Imk*Ink)
                    fill_symmetric_terms(m, k, n, k, Mmknk)

                # For l!=k it we need to cover n=0..max_azi because of eq. 2.38
                # but we just need to cover l=k+1..max_rad, because of eq. 2.39
                for n in range(max_azi+1):
                    for l in range(k+1, max_rad+1):
                        Inl = spectrum[(abs(n), l)]
                        Mmknl = (1j)**(m-n)*_np.dot(Z_wp, Imk*Inl)
                        fill_symmetric_terms(m, k, n, l, Mmknl)

        return cls._reshape_coupling_matrix(M)

    def reduced_longitudinal_mode_coupling(
            self, budget=None,  element=None, w=None, Zl=None,
            max_azi=10, max_rad=12, cbmode=0, modecoup_matrix=None):
        """Calculate the longitudinal mode-coupling eigen-values.

        This implementation uses a symmetry of the mode-coupling matrix to
        make a dimension reduction of the problem. It differs from the method
        `longitudinal_mode_coupling` because the azimuthal expansion only
        consider m>0 modes. However, this consideration excludes the
        possibility of including the fokker planck terms, since they don't
        have this symmetry.

        The eigen-values returned here are normalized by the synchrotron
        frequency, which means they are adimensional. Besides, the
        implementation does not guarantee any specific ordering.

        According to ref. [1], the stability will occur if all eigen-value are
        real and positive.

        We follow ref.[1] closely, such that the comments in the
        implementation refer to equation of this reference.

        We make a small generalization of ref. [1], by considering that the
        ring is filled with several equally spaced identical bunches and that
        the eigen-values are being calculated for a given coupled-bunch mode
        of the beam.

        References:
            [1] Suzuki, T., Chin, Y.-H., & Satoh, K. (1983). Mode-coupling
                theory and bunch lengthening in spear. Particle Accelerators,
                13, 179-198.

        Attributes and properties used by this method:
            total_current
            energy
            rev_ang_freq
            sync_tune
            bunlen
            mom_comp
            num_bun

        The size N of each dimension of the mode-coupling and fokker plank
        matrices is given by:
            N = max_azi*(max_rad+1))

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            Zl (numpy.ndarray (M, ), optional): longitudinal impedance in
                units of [Ohm] evaluated at angular frequencies in `w`.
                Only considered if `budget` and `element` are `None`.
                Defaults to None.
            max_azi (int, optional): Maximum azimuthal mode to consider in the
                truncation. The solution will consider all `max_azi`
                azimuthal modes whose absolute value are <= than `max_azi`.
                Defaults to 3.
            max_rad (int, optional): Maximum radial mode to consider in the
                truncation. The solution will consider all `max_rad + 1`
                radial modes <= than `max_rad`. Defaults to 4.
            cbmode (int, optional): Coupled-bunch mode for which the
                eigen-values will be found. Must be lower than `num_bun`.
                Defaults to 0.
            modecoup_matrix (numpy.ndarray, (N, N), optional): the
                mode-coupling matrix to be used in calculations. If None, then
                it will be calculated internally. Defaults to None.

        Returns:
            eigenvals (numpy.ndarray, (N, )): squared normalized eigen-modes
                of the mode-coupling problem.
            modecoup_matrix (numpy.ndarray, (N, N)): the mode-coupling matrix
                used in calculations.

        """
        E = self.energy
        w0 = self.rev_ang_freq
        bunlen = self.bunlen
        sync_tune = self.sync_tune
        eta = self.mom_comp
        nb = self.num_bun
        I_b = self.total_current / nb
        if cbmode >= nb:
            cbmode = 0
            print(
                'Coupled bunch mode greater than number of bunchs.\n',
                'Reseting cbmode to 0.')

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        # There is an approximation here which is only valid for broad-band
        # impedances. I sample the impedance only at the azimuthal mode m=1.
        wp = self._get_sampling_ang_freq(
            w[0], w[-1], w0, nb, 1, sync_tune[0], [cbmode])

        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)
        Zl_wp = Zl_interp / wp

        # Calculate the mode coupling matrix:
        if modecoup_matrix is None:
            modecoup_matrix = self._calc_vlasov_reduced(
                Zl_wp, wp, bunlen, max_azi, max_rad)
        # Calculate the current independent diagonal matrix:
        ms = _np.arange(1, max_rad+1)
        D = _np.einsum('mn,kl->mknl', _np.diag(ms*ms), _np.eye(max_rad+1))
        D = self._reshape_coupling_matrix(D)

        # We separated the K value from the defition of M so that M
        # could be current independent in case of a constant bunch
        # length. To understand the equation for K implemented here,
        # please look at eqs. 2.49 and 2.52 of ref. [2]. and notice
        # that the definition of M_mlnh in 2.26 of ref. [1] has Z/p
        # instead of Z/wp as we have here:
        xi = eta * I_b/sync_tune**2/E
        sig_theta = w0*bunlen/_LSPEED
        K = xi * w0/_np.pi/(sig_theta)**2
        # Lastly we need to multiply by `nb` because of our
        # generalization to include n equally spaced bunches.
        K *= nb

        # Please, check eq. 2.52 of ref. [1]:
        A = D + 1j*K*modecoup_matrix

        return _np.linalg.eigvals(A), modecoup_matrix

    @classmethod
    def _calc_vlasov_reduced(cls, Z_wp, wp, bunlen, max_azi, max_rad):
        """Calculate the longitudinal mode-coupling matrix.

        We follow ref.[1] closely, such that the comments in the
        implementation refer to equation of this reference.

        The mode-coupling matrix is actually an object with 4 indices: Mmknl,
        because it relates the expansion coefficients of the oscillation modes,
        which are 2-indices objects, a_mk, being m for the azimuthal modes
        and k for the radial modes.

        However, after filling this 4-indices "matrix" we cast it in a
        regular 2D matrix by mapping the m and k coefficients into a single
        dimension, so we can solve the eigen-value problem with standard
        techniques.

        References:
            [1] Suzuki, T., Chin, Y.-H., & Satoh, K. (1983). Mode-coupling
                theory and bunch lengthening in spear. Particle Accelerators,
                13, 179-198.

        Args:
            Z_wp (numpy.ndarray, (N, )): Impedance divided sampled at wp
                angular frequency and divided by wp. Unit in [Ohm/(rad/s)]
            wp (numpy.ndarray, (N, )): angular frequency in units of [rad/s]
            bunlen (float): bunch length in units of [m]
            max_azi (int): maximum azimuthal mode to consider.
            max_rad (int): maximum radial mode to consider.

        """
        def fill_symmetric_terms(m, k, n, l, Mmknl):
            # by 2.26 that Mnlmk = Mmknl * (-1)**(m-n):
            M[m-1, k, n-1, l] = m*m*Mmknl
            M[n-1, l, m-1, k] = n*n*Mmknl * (-1)**(m-n)

        M = _np.zeros([max_azi, max_rad+1, max_azi, max_rad+1], dtype=complex)
        spectrum = cls.calc_spectrum(wp, bunlen, max_rad, max_azi)
        for m in range(1, max_azi+1):
            for k in range(max_rad+1):
                Imk = spectrum[(abs(m), k)]
                # For l=k it is only necessary to cover n=m..max_azi
                # because of the symmetry relations of eq. 2.46:
                for n in range(m, max_azi+1):
                    Ink = spectrum[(abs(n), k)]
                    Mmknk = (1j)**(m-n)*_np.dot(Z_wp, Imk*Ink)
                    fill_symmetric_terms(m, k, n, k, Mmknk)

                # For l!=k it we need to cover n=0..max_azi because of eq. 2.38
                # but we just need to cover l=k+1..max_rad, because of eq. 2.39
                for n in range(1, max_azi+1):
                    for l in range(k+1, max_rad+1):
                        Inl = spectrum[(abs(n), l)]
                        Mmknl = (1j)**(m-n)*_np.dot(Z_wp, Imk*Inl)
                        fill_symmetric_terms(m, k, n, l, Mmknl)

        return cls._reshape_coupling_matrix(M)

    def transverse_mode_coupling(
            self, budget=None, element=None, w=None, Zt=None, plane='y',
            max_azi=3, max_rad=4, cbmode=0, use_fokker=True,
            modecoup_matrix=None, fokker_matrix=None):
        """Calculate the transverse mode-coupling eigen-values.

        The eigen-values returned here are normalized by the synchrotron
        frequency, which means they are adimensional. Besides, the
        implementation does not guarantee any specific ordering.

        To get the tune-shifts you must take the real part of the eigen-values.
        To get the growth rate in [1/s] you must multiply the imaginary part
        by the angular synchrotron frequency `sync_tune*w0`.
        Instability occurs if any of the growth rates are larger than zero.

        We follow ref.[1] closely, such that the comments in the
        implementation refer to equation of this reference. Ref. [2] gives
        identical results of ref. [1] and ref.[3] makes a different expansion
        of oscillation modes, so we leave it here only for the interested user.

        We make a small generalization of ref. [1], by considering that the
        ring is filled with several equally spaced identical bunches and that
        the eigen-values are being calculated for a given coupled-bunch mode
        of the beam.

        References:
            [1] Chin, Y. H. (1985). Transverse Mode Coupling Instabilities in
                the SPS. In CERN/SPS/85-2 (DI-MST).
            [2] Suzuki, T. (1986). Fokker-planck theory of transverse
                mode-coupling instability. Particle Accelerators, 20, 79-96.
            [3] Lindberg, R. R. (2016). Fokker-Planck analysis of transverse
                collective instabilities in electron storage rings. Phys. Rev.
                Accel. Beams, 19(12), 124402.
                https://doi.org/10.1103/PhysRevAccelBeams.19.124402

        Attributes and properties used by this method:
            total_current
            energy
            rev_ang_freq
            sync_tune
            bunlen
            mom_comp
            num_bun
            dampte
            damptx
            dampty
            tunex
            tuney
            chromx
            chromy

        The size N of each dimension of the mode-coupling and fokker plank
        matrices is given by:
            N = (2*max_azi+1)*(max_rad+1))

        Args:
            budget (pycolleff.impedances.Budget, optional): Impdedance budget
                of the machine. If not None, `element`, `w` and `Zt` inputs
                are ignored. Defaults to None.
            element (pycolleff.impedances.Element, optional): Impedance
                element. Only considered if `budget is `None`. If not `None`,
                `w` and `Zt` are ignored. Defaults to None.
            w (numpy.ndarray, (M, ), optional): Angular frequency in [rad/s].
                Only considered if `budget` and `element` are `None`. Defaults
                to `None`.
            Zt (numpy.ndarray (M, ), optional): Beta weighted transverse
                impedance in [Ohm] evaluated at angular frequencies in `w`.
                Only considered if `budget` and `element` are `None`. Defaults
                to None.
            plane (str, optional): Plane where to evaluate the analysis. May
                be in {'x', 'y'}. Defaults to 'y'.
            max_azi (int, optional): Maximum azimuthal mode to consider in the
                truncation. The solution will consider all `2*max_azi + 1`
                azimuthal modes whose absolute value are <= than `max_azi`.
                Defaults to 3.
            max_rad (int, optional): Maximum radial mode to consider in the
                truncation. The solution will consider all `max_rad + 1`
                radial modes <= than `max_rad`. Defaults to 4.
            cbmode (int, optional): Coupled-bunch mode for which the
                eigen-values will be found. Must be lower than `num_bun`.
                Defaults to 0.
            use_fokker (bool, optional): Whether or not to include the
                mode-coupling fokker-planck terms described in ref. [2].
                If False, only the diagonal damping terms will be included.
                Defaults to True.

            modecoup_matrix (numpy.ndarray, (N, N), optional): the
                mode-coupling matrix to be used in calculations. If None, then
                it will be calculated internally. Defaults to None.
            fokker_matrix (numpy.ndarray, (N, N), optional): the fokker planck
                matrix to be used in calculations. If None, then it will be
                calculated internally. Defaults to None.

        Returns:
            eigenvals (numpy.ndarray, (N, ): normalized eigen-modes of the
                mode-coupling problem.
            modecoup_matrix (numpy.ndarray, (N, N)): the mode-coupling matrix
                used in calculations.
            fokker_matrix (numpy.ndarray, (N, N)): the fokker planck matrix
                used in calculations.

        """
        energy = self.energy
        w0 = self.rev_ang_freq
        bunlen = self.bunlen
        sync_tune = self.sync_tune
        eta = self.mom_comp
        nb = self.num_bun
        I_b = self.total_current / nb
        alpe = 1/self.dampte
        if cbmode >= nb:
            cbmode = 0
            print(
                'Coupled bunch mode greater than number of bunchs.\n',
                'Reseting cbmode to 0.')

        if plane.lower().startswith(('x', 'h')):
            alpt, nut = 1/self.damptx, self.tunex
            chrom, imp = self.chromx, 'Zdx'
        else:
            alpt, nut = 1/self.dampty, self.tuney
            chrom, imp = self.chromy, 'Zdy'

        w, Zt = self._prepare_input_impedance(budget, element, w, Zt, imp)

        # There is an approximation here which is only valid for broad-band
        # impedances. I sample the impedance only at the azimuthal mode m=1.
        wp = self._get_sampling_ang_freq(
            w[0], w[-1], w0, nb, 1, sync_tune, [cbmode], nut)
        w_crom = chrom/eta * w0
        wpcro = wp - w_crom

        Zt_wp = self._get_interpolated_impedance(wp, w, Zt)

        # Calculate the fokker-planck matrix:
        if fokker_matrix is None:
            fokker_matrix = self._calc_fokker_planck(
                max_azi, max_rad, alpe, use_fokker=use_fokker, alpt=alpt)
        # Calculate the mode coupling matrix:
        if modecoup_matrix is None:
            modecoup_matrix = self._calc_vlasov_transverse(
                Zt_wp, wpcro, bunlen, max_azi, max_rad)
        # D is the current independent diagonal matrix of eq. 2.37:
        ms = _np.arange(-max_azi, max_azi+1)
        D = _np.einsum('mn,kl->mknl', _np.diag(ms), _np.eye(max_rad+1))
        D = self._reshape_coupling_matrix(D)

        # We separated the K value from the defition of M so that M
        # could be current independent in case of a constant bunch
        # length. To understand the equation for K implemented here,
        # please look at eqs. 2.35 and 2.41 of ref. [1]., or to eq. 37
        # of ref. [2].
        K = I_b * w0/(4*_np.pi)/(sync_tune * w0)/energy
        # Lastly we need to multiply by `nb` because of our
        # generalization to include n equally spaced bunches.
        K *= nb

        # This is the complete matrix described in eq. 2.37 of ref.
        # [1] with the aditional fokker-plank terms of eq. 37 of
        # ref. [2] normalized by ws:
        A = D + K*modecoup_matrix + 1j*fokker_matrix/(sync_tune * w0)
        return _np.linalg.eigvals(A), modecoup_matrix, fokker_matrix

    @classmethod
    def _calc_vlasov_transverse(cls, Z_wp, wpcro, bunlen, max_azi, max_rad):
        """Calculate the transverse mode-coupling matrix.

        We follow ref.[1] closely, such that the comments in the
        implementation refer to equation of this reference. Ref. [2] gives
        identical results of ref. [1] and ref.[3] makes a different expansion
        of oscillation modes, so we leave it here only for the interested user.

        The mode-coupling matrix is actually an object with 4 indices: Mmknl,
        because it relates the expansion coefficients of the oscillation modes,
        which are 2-indices objects, a_mk, being m for the azimuthal modes
        and k for the radial modes.

        However, after filling this 4-indices "matrix" we cast it in a
        regular 2D matrix by mapping the m and k coefficients into a single
        dimension, so we can solve the eigen-value problem with standard
        techniques.

        References:
            [1] Chin, Y. H. (1985). Transverse Mode Coupling Instabilities in
                the SPS. In CERN/SPS/85-2 (DI-MST).
            [2] Suzuki, T. (1986). Fokker-planck theory of transverse
                mode-coupling instability. Particle Accelerators, 20, 79-96.
            [3] Lindberg, R. R. (2016). Fokker-Planck analysis of transverse
                collective instabilities in electron storage rings. Phys. Rev.
                Accel. Beams, 19(12), 124402.
                https://doi.org/10.1103/PhysRevAccelBeams.19.124402

        Args:
            Z_wp (numpy.ndarray, (N, )): Beta weigthed impedance sampled at wp
                angular frequency [Ohm]
            wpcro (numpy.ndarray, (N, )): chromatic shifted angular frequency
                in units of [rad/s]
            bunlen (float): bunch length in units of [m]
            max_azi (int): maximum azimuthal mode to consider.
            max_rad (int): maximum radial mode to consider.

        """
        def fill_symmetric_terms(m, k, n, l, Mmknl):
            # Please, note the symmetry relations in eq. 2.38 e 2.39:
            M[max_azi+m, k, max_azi+n, l] = Mmknl
            M[max_azi-m, k, max_azi+n, l] = Mmknl
            M[max_azi+m, k, max_azi-n, l] = Mmknl
            M[max_azi-m, k, max_azi-n, l] = Mmknl
            Mnlmk = Mmknl * (-1)**(m-n)
            M[max_azi+n, l, max_azi+m, k] = Mnlmk
            M[max_azi+n, l, max_azi-m, k] = Mnlmk
            M[max_azi-n, l, max_azi+m, k] = Mnlmk
            M[max_azi-n, l, max_azi-m, k] = Mnlmk

        M = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        spectrum = cls.calc_spectrum(wpcro, bunlen, max_rad, max_azi)
        for m in range(max_azi+1):
            for k in range(max_rad+1):
                Imk = spectrum[(abs(m), k)]
                # For l=k it is only necessary to cover n=m..max_azi
                # because of the symmetry relations of eq. 2.38 and 2.39:
                for n in range(m, max_azi+1):
                    Ink = spectrum[(abs(n), k)]
                    Mmknk = -1j*(1j)**(m-n)*_np.dot(Z_wp, Imk*Ink)
                    fill_symmetric_terms(m, k, n, k, Mmknk)

                # For l!=k it we need to cover n=0..max_azi because of eq. 2.38
                # but we just need to cover l=k+1..max_rad, because of eq. 2.39
                for n in range(max_azi+1):
                    for l in range(k+1, max_rad+1):
                        Inl = spectrum[(abs(n), l)]
                        Mmknl = -1j*(1j)**(m-n)*_np.dot(Z_wp, Imk*Inl)
                        fill_symmetric_terms(m, k, n, l, Mmknl)

        return cls._reshape_coupling_matrix(M)

    @classmethod
    def _calc_fokker_planck(
            cls, max_azi, max_rad, alpe, use_fokker=True, alpt=0):
        """Calculate Fokker-Planck Matrix.

        This implementation is based on the theory developed in refs [1,2].
        Note that the expansion made in ref.[3] is different than the one
        performed at refs.[1,2]. Besides the matrix terms are simpler and the
        author claim they are equal if we expand the ones of the other refs.,
        but here we decided to follow refs.[1,2].

        To understand the implementation, please follow eq. 37 of ref.[1].
        Note that this equation is a matrix equation for the coefficients
        a_ml of the expansion of the perturbed distribution.

        Since these coefficients have to symbols, m and l, being natural to
        think of them as 2D matrices, it is also natural to think of a four
        indices "matrix" to compose the eigen-value problem, relating the
        a_ml to a_kn.
        However, after filling this four indices "matrix" we cast it in a
        regular 2D matrix by vectorizing the m and l coefficients so we can
        solve the eigen-value problem with standard techniques.

        Refs:
            [1] Suzuki, T. (1986). Fokker-planck theory of transverse
                mode-coupling instability. Particle Accelerators, 20, 79-96.
            [2] Suzuki, T. (1983). Theory of longitudinal bunched-beam
                instabilities based on the fokker-planck equation. Particle
                Accelerators, 14, 91-108.
            [3] Lindberg, R. R. (2016). Fokker-Planck analysis of transverse
                collective instabilities in electron storage rings. Phys. Rev.
                Accel. Beams, 19(12), 124402.
                https://doi.org/10.1103/PhysRevAccelBeams.19.124402

        Args:
            max_azi (int): Maximum azimuthal number before truncation.
            max_rad (int): Maximum radial number before truncation.
            alpe (float): longitudinal damping rate [1/s]
            alpt (int, optional): transverse damping rate [1/s]. Defaults to 0.
            use_fokker (bool, optional): Whether or not to calculate the
                mode-coupling fokker planck terms. If False returns only the
                diagonal damping terms. Defaults to True.

        Returns:
            (numpy.ndarray, (
                (2*max_azi+1)x(max_rad+1), (2*max_azi+1)x(max_rad+1))):
                    Square Fokker-Planck matrix that relates a_ml to a_kn.

        """
        F = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        for m in range(-max_azi, max_azi+1):
            amm2 = abs(m-2)
            amp2 = abs(m+2)
            am = abs(m)
            for l in range(max_rad+1):
                F[max_azi+m, l, max_azi+m, l] = - alpt - alpe*(abs(m) + 2*l)
                if not use_fokker:
                    continue

                Fml = F[max_azi+m, l]
                for k in range(max_rad+1):
                    if m-2 >= -max_azi:
                        tm1 = cls._kappa(amm2, am, k, l, 0)
                        tm2 = cls._kappa(amm2, am, k, l, -1)
                        tm3 = cls._kappa(amm2, am, k-1, l, -1)
                        tm1 *= (k + amm2/2 - (m-2)/2) / 2
                        tm2 *= -(m-1)/4 * (m-2 - amm2 - 2*k)
                        tm3 *= -(m-1)/2 * _sqrt(k*(amm2 + k))
                        Fml[max_azi+m-2, k] = (tm1 + tm2 + tm3) * 2 * alpe
                    if m+2 <= max_azi:
                        tm1 = cls._kappa(amp2, am, k, l, 0)
                        tm2 = cls._kappa(amp2, am, k, l, -1)
                        tm3 = cls._kappa(amp2, am, k-1, l, -1)
                        tm1 *= (k + amp2/2 + (m+2)/2)/2
                        tm2 *= -(m+1)/4 * (m+2 + amp2 + 2*k)
                        tm3 *= +(m+1)/2 * _sqrt(k*(amp2 + k))
                        Fml[max_azi+m+2, k] = (tm1 + tm2 + tm3) * 2 * alpe
        return cls._reshape_coupling_matrix(F)
        # return F

    @classmethod
    def _kappa(cls, alpha, beta, k, l, N):
        """Calculate Kappa factor for Fokker-Planck Matrix.

        [1] Suzuki, T. (1986). Fokker-planck theory of transverse
            mode-coupling instability. Particle Accelerators, 20, 79-96.
        [2] Suzuki, T. (1983). Theory of longitudinal bunched-beam
            instabilities based on the fokker-planck equation. Particle
            Accelerators, 14, 91-108.
        [3] Wikipedia contributors. Laguerre polynomials. Wikipedia, The Free
            Encyclopedia. April 5, 2022, 20:42 UTC. Available at:
            https://en.wikipedia.org/w/index.php?title=Laguerre_polynomials&oldid=1081186293.
            Accessed April 14, 2022.

        This implementation follow eqs. 42-45 of ref.[1] plus the
        symmetry relation that can be infered from eq. 39.

        Args:
            alpha (int): alpha from ref.[1]
            beta (int): beta from ref.[1]. Must be in {alpha, alpha+2}.
            k (int): k from ref.[1]
            l (int): l from ref.[1]
            N (int): N from ref.[1]. Must be in {0, -1}

        Raises:
            ValueError: when beta not in {alpha-2, alpha, alpha+2} or
                N not in {0, -1}.

        Returns:
            float: kappa value

        """
        kappa = None
        # Refs [1-2] doesn't mention anything about this, but the Laguerre
        # polynomials are not defined for negative indices. Looking at
        # ref. [3] I noticed there is no non-trivial solution for n<0, so I
        # considered the integral that originates kappa is also zero in this
        # case:
        if k < 0:
            kappa = 0
        elif l < 0:
            kappa = 0
        elif alpha - 2 == beta:
            # note the symmetry relation in eq. 39 or ref.[1]
            kappa = cls._kappa(beta, alpha, l, k, N)
        elif alpha == beta:
            if not N:
                kappa = 1 if k == l else 0
            elif N == -1:
                kappa = _factorial(alpha+l) / _factorial(alpha+k)
                kappa *= _factorial(k) / _factorial(l)
                kappa = _sqrt(kappa)
                if k <= l:
                    kappa = 1/kappa
                kappa /= alpha
        elif alpha + 2 == beta:
            if not N:
                kappa = 0
                if k <= l:
                    kappa = _factorial(alpha+k) / _factorial(alpha+l+2)
                    kappa *= _factorial(l) / _factorial(k)
                    kappa = _sqrt(kappa)
                    kappa *= alpha + 1
                elif k == l + 1:
                    kappa = -_sqrt((l+1)/(alpha+l+2))
            elif N == -1:
                kappa = 0
                if k <= l:
                    kappa = _factorial(alpha+k) / _factorial(alpha+l+2)
                    kappa *= _factorial(l) / _factorial(k)
                    kappa = _sqrt(kappa)
                    kappa *= l + 1 - k
        if kappa is None:
            raise ValueError(
                'beta must be in {alpha-2, alpha, alpha+2} and N in {0, -1}.')
        return kappa

    @staticmethod
    def calc_spectrum(wp, bunlen, max_rad=4, max_azi=3, only=False):
        """Calculate the bunch spectrum based on eq. 2.43 of ref [1].

        It calculates all bunch spectra up the radial mode defined by `max_rad`
        and the azimuthal mode defined by `max_azi`.

        This implementation does not add the sign defined by the epsilon_m term
        of eq. 2.43. This is done to speedup the calculation of the terms of
        the mode-coupling matrix.

        References:
            [1] Chin, Y. H. (1985). Transverse Mode Coupling Instabilities in
                the SPS. In CERN/SPS/85-2 (DI-MST).

        Args:
            wp (numpy.ndarray, (N, )): frequencies where to calculate the
                spectrum [rad/s]
            bunlen (float): bunch length in [m].
            max_rad (int, optional): Number of radial modes. Defaults to 4.
            max_azi (int, optional): Number of azimuthal modes. Defaults to 3.
            only (bool, optional): This flag indicates whether only the mode
                `m=max_azi` and `k=max_rad` will the returned.
                Defaults to False.

        Returns:
            dict: Each key, defined by `(abs(m) + k)`, has as value the bunch
                spectrum for oscillation mode with azimuthal mode `m` and
                radial mode `k`. Each bunch spectrum ins a `numpy.ndarray` of
                shape (N, ).

        """
        def my_pow(vetor, n):
            pows = dict()
            res = _np.ones(vetor.shape, dtype=float)
            for i in range(n):
                pows[i] = res.copy()
                res *= vetor
            pows[n] = res
            return pows

        max_azi0, max_rad0 = 0, 0
        if only:
            max_azi0, max_rad0 = max_azi, max_rad
        sigW = wp*bunlen/_LSPEED/_sqrt(2)
        spectrum = dict()
        powers = my_pow(sigW, max(max_azi + 2*max_rad, 2))
        expo = _np.exp(-powers[2])
        for azi in range(max_azi0, max_azi+1):
            for rad in range(max_rad0, max_rad+1):
                chave2 = abs(azi) + 2*rad
                fact = _factorial(abs(azi) + rad) * _factorial(rad)
                fact = 1/_sqrt(fact)
                chave = (abs(azi), rad)
                spectrum[chave] = fact * powers[chave2] * expo
        if only:
            spectrum = spectrum[chave]
        return spectrum

    def kicker_power(
            self, gr, Rshunt=15e3, betak=5, Ab=1e-3, betab=5, cbmode=None,
            plane='long'):
        '''Calculate the minimum transverse bunch by bunch feedback power
            necessary to damp coupled-bunch oscillations, assuming perfect
            tunning of the system.
            Formula taken from CAS_Digital_Signal_Processing, pag. 314.

        INPUTS:
          gr     : growth rate you want to damp [Hz].
          Rshunt : shunt resistence of the kicker [Ohm].
          betak  : betatron function at kicker position [m].
          betab  : betatron function at bpm postiion [m].
          Ab     : Initial amplitude of the oscilation at the BPM position [m].
          cbmode : coupled bunch mode of the instability. Not used yet.

        OUTPUT: RF power needed to damp instability [W].
        '''
        if plane.lower().startswith('l'):
            P = 2*_np.pi*self.sync_tune * self.energy / self.mom_comp
            P *= gr * Ab / _LSPEED
            P *= P
        else:
            P = 1/betak * self.energy**2 * (self.rev_time*gr)**2
            P *= (Ab**2/betab)
        return P*2/Rshunt

    def calc_energy_spread_increase(
            self, N=10000000, bunlen=3e-3, wake=1e17, larg=50e-6, Ib=1e-3,
            zz=None):
        """."""
        zz = zz or bunlen
        po = _np.exp(-zz*zz/2/bunlen/bunlen)/_np.sqrt(2*_np.pi)/bunlen*larg

        var1 = po*(1-po)*N

        # radiation kick
        fde = (1 - self.rev_time/self.dampte/2)
        srde = _np.sqrt(1-fde*fde)*self.espread(Ib)

        # wake kick
        factor = _np.sqrt(var1)/N
        kick = self.rev_time/self.energy * Ib * wake * factor
        # kick = np.sqrt(var1)*1.6e-19/3e9 * Wake

        return _np.sqrt(kick**2/(1-fde*fde))

    def _prepare_input_impedance(self, budget, element, w, Z, imp='Zll'):
        if budget is not None:
            return budget.ang_freq, getattr(budget, imp)
        elif element is not None:
            return element.ang_freq, getattr(element, imp)
        elif w is not None and Z is not None:
            return w, Z
        else:
            raise Exception('Incorrect impedance input.')

    @staticmethod
    def _get_sampling_ang_freq(
            wmin, wmax, w0, nb, m=0, sync_tune=0, cbmodes=(0, ), nut=0,
            return_p=False):
        wp = []
        cbmodes = _np.array(cbmodes)
        shift = cbmodes + m*sync_tune + nut
        pmin = _np.ceil((wmin/w0 - shift.min())/nb)
        pmax = _np.floor((wmax/w0 - shift.max())/nb)
        p = _np.arange(pmin, pmax+1)
        wp = nb*p[None, :] + shift[:, None]
        wp = wp.squeeze()
        wp *= w0
        if return_p:
            return wp, p
        return wp

    @staticmethod
    def _get_interpolated_impedance(wp, w, Z):
        if callable(Z):
            Z_interp = Z(wp)
        else:
            # Complex interpolation is ill-defined because it gives different
            # results depending whether you do it in polar or cartesian
            # coordinates, so we must do it by parts.
            # Imaginary must come first, so that interpolated variable is
            # recognized as of complex type:
            Z_interp = _np.interp(wp, w, Z.imag) * 1j
            Z_interp += _np.interp(wp, w, Z.real)
        return Z_interp

    @staticmethod
    def _reshape_coupling_matrix(M):
        shap = M.shape
        siz = shap[0] * shap[1]
        return M.swapaxes(0, 3).swapaxes(1, 2).reshape([siz, siz]).transpose()
