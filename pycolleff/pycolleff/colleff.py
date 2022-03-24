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
        self.circ = 0.0  # ring circumference in [m]
        self.T0 = 0.0  # revolution period in [s]
        self.f0 = 0.0  # revolution frequency [Hz]
        self.w0 = 0.0  # revolution angular frequency [rad/s]
        self.mom_cmpct = 0.0  # momentum compaction factor
        self.E = 0.0  # energy [eV]
        self.nuy = 0.0  # vertical tune
        self.nux = 0.0  # horizontal tune
        self.chromx = 0.0  # horizontal chromaticity
        self.chromy = 0.0  # vertical chromaticity
        self.harm_num = 1    # harmonic Number
        self.nbun = 1    # number of bunches filled

        # bunch current vector in [A]
        self.cur_bun = _np.linspace(0, 4, num=40)*0.0
        self.nom_cur = 0.00  # total current [A]
        self.nus = 0.00  # synchrotron tune
        self.espread = 0.0
        self.bunlen = 0.0  # bunch length in [m]
        self.emitx = 0.0
        self.emity = 0.0
        self.damptx = 0.0  # horizontal damping time in [s]
        self.dampty = 0.0  # vertical damping time in [s]
        self.dampte = 0.0  # longitudinal damping time in [s]
        self.en_lost_rad = 0.0  # Energy lost per turn in [eV]

    def __str__(self):
        """."""
        tmpl_st = '{0:28s}: {1:^20s}\n'
        tmpl_f = '{0:28s}: {1:^20.3f}\n'
        string = ''
        string += tmpl_st.format('Lattice Version', self.version)
        string += tmpl_f.format('Circumference [m]', self.circ)
        string += tmpl_f.format('Revolution Period [us]', self.T0*1e6)
        string += tmpl_f.format('Revolution Frequency [kHz]', self.f0/1e3)
        string += tmpl_f.format('Energy [GeV]', self.E/1e9)
        string += '{0:28s}: {1:^20.2e}\n'.format(
            'Momentum Compaction', self.mom_cmpct)
        string += '{0:28s}: {1:^20d}\n'.format(
            'Harmonic Number', self.harm_num)
        string += tmpl_f.format('Current [mA]', self.nom_cur*1e3)
        string += tmpl_f.format(
            'Current per Bunch [mA]', self.nom_cur/self.nbun*1e3)
        string += '{0:28s}: {1:^20.3e}\n'.format('Synchrotron Tune', self.nus)
        string += '{0:28s}: {1:>9.3f}/{2:<10.3f}\n'.format(
            'Tunes x/y', self.nux, self.nuy)
        string += '{0:28s}: {1:>6.1f}/{2:^6.1f}/{3:<6.1f}\n'.format(
            'Damping Times x/y/e [ms]', self.damptx*1e3, self.dampty*1e3,
            self.dampte*1e3)
        string += '{0:28s}: {1:^20.2e}\n'.format('Energy Spread', self.espread)
        string += '{0:28s}: {1:^20.2e}\n'.format(
            'Bunch Length [mm]', self.bunlen*1e3)
        return string

    def loss_factor(
            self, budget=None, element=None, w=None, Zl=None, bunlen=None):
        """Calculate the loss factor and effective impedance.

        Inputs:
          budget  = instance of Budget class
                        or
          element = instance of Element class
                        or
          w  = angular frequency [rad/s]
          Zl = Longitudinal Impedance [Ohm]

          (optional) bunlen = Longitudinal beamsize [m]

        Outputs:
          lossf  = Loss factor in V/C
          Pl     = Power loss in W
          Zl_eff = Effective Impedance in Ohm
          wp     = vector of angular frequencies where the impedance
                   was sampled
          lossfp = vector of loss factor with contribution of each sampled
                   frequency
        """
        w0 = self.w0
        nb = self.nbun
        bunlen = bunlen or self.bunlen(self.nom_cur/nb)

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        limw = 5 * _LSPEED / bunlen
        maxw = min(w[-1], limw)
        minw = max(w[0], -limw)

        pmin = _np.ceil(minw / (w0*nb))   # arredonda em direcao a +infinito
        pmax = _np.floor(maxw / (w0*nb))  # arredonda em direcao a -infinito
        p = _np.arange(pmin, pmax+1)
        wp = w0*p*nb

        # h = _np.exp(-(wp*bunlen/_LSPEED)**2)
        specs = self.calc_spectrum(wp, bunlen, max_rad=0, max_azi=1)
        h = specs[(0, 0)]**2
        interpol_Z = _np.interp(wp, w, Zl.real)

        # Loss factor and Power loss:
        lossfp = nb*(w0/2/_np.pi)*interpol_Z*h
        lossf = sum(lossfp)
        Pl = lossf * (self.T0*self.nom_cur**2/nb)

        # Effective Impedance: pag 223 K.Y.Ng book
        interpol_Z = _np.interp(wp, w, Zl.imag)
        h = specs[(0, 0)]**2
        Zl_eff = sum(interpol_Z*h/((wp+1e-4)/w0)) / sum(h)

        return lossf, Pl, Zl_eff, wp, lossfp

    def kick_factor(
            self, budget=None, element=None, w=None, Z=None, bunlen=None,
            Imp='Zdy'):
        """Calculate the kick factor, tune shift and effective impedance.

        Inputs:
          budget  = instance of Budget class
                        or
          element = instance of Element class
                        or
          w  = angular frequency [rad/s]
          Zl = Longitudinal Impedance [Ohm]

          (optional) bunlen = Longitudinal beamsize [m]
          (optional) Imp   = string with type of impedance to consider
                             from element or budget
                             options: Zll(default), Zdy, Zqy, Zdx, Zqx

        Outputs:
          Kick_f = Kick Factor factor in V/C
          Tush   = Tune shift
          Zt_eff = Effective Impedance in Ohm
        """
        if Imp.lower().startswith(('zqx', 'zqy')):
            nut = 0
        elif Imp.lower().startswith('zdx'):
            nut = self.nux % 1
        else:
            nut = self.nuy % 1

        w0 = self.w0
        nb = self.nbun
        bunlen = bunlen or self.bunlen(self.nom_cur/nb)

        w, Z = self._prepare_input_impedance(budget, element, w, Z, Imp)

        limw = 5 * _LSPEED / bunlen
        maxw = min(w[-1], limw)
        minw = max(w[0], -limw)

        pmin = _np.ceil(minw / (w0*nb))   # arredonda em direcao a +infinito
        pmax = _np.floor(maxw / (w0*nb))  # arredonda em direcao a -infinito
        p = _np.arange(pmin, pmax+1)
        wp = w0*(p*nb + nut)

        specs = self.calc_spectrum(wp, bunlen, max_rad=0, max_azi=1)
        h = specs[(0, 0)]**2
        interpol_Z = _np.interp(wp, w, Z.imag)

        Kick_f = nb*(w0/2/_np.pi) * sum(interpol_Z*h)
        Tush = (self.nom_cur/nb)/(2*self.E) * Kick_f / w0

        h = specs[(1, 0)]**2
        Zt_eff = sum(interpol_Z*h) / sum(h)  # pag 383 K.Y.Ng Book

        return Kick_f, Tush, Zt_eff

    def _prepare_input_impedance(self, budget, element, w, Z, imp='Zll'):
        if budget is not None:
            return budget.w, getattr(budget, imp)
        elif element is not None:
            return element.w, getattr(element, imp)
        elif w is not None and Z is not None:
            return w, Z
        # elif self.budget is not None:
        #     return self.budget.w, getattr(self.budget,imp)
        else:
            raise Exception('Incorrect impedance input.')

    def longitudinal_cbi(
            self, budget=None, element=None, w=None,
            Zl=None, bunlen=None, m=1, inverse=False, full=False):
        """Calculate the complex coeherent frequency shifts of the beam for all
        the oscilation modes, considering a Gaussian beam and only azimuthal
        mode k=0;
        """

        if abs(m) <= 0:
            raise ValueError('azimuthal mode m must be greater than zero.')
        nus = self.nus
        w0 = self.w0
        eta = self.mom_cmpct
        E = self.E
        nb = self.nbun
        I_tot = self.nom_cur
        alpe = 1/self.dampte/nus/w0
        if bunlen is None:
            bunlen = self.bunlen

        if inverse:
            h = self.calc_spectrum(w, bunlen, max_rad=0, max_azi=m, only=True)**2
            Z = w*alpe
            Z /= (I_tot*w0*eta/(2*_np.pi)/(nus*w0)**2/E/(bunlen/_LSPEED)**2)*h
            return Z

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        # Calculate Effective Impedance
        cbmodes = _np.arange(nb)
        wp = self._get_sampling_ang_freq(w, w0, nb, m, nus, cbmodes)
        h = self.calc_spectrum(wp, bunlen, max_rad=0, max_azi=m, only=True)**2
        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)
        Zl_eff = (Zl_interp/wp * h).sum(1)

        # Returns the relative Tune Shift:
        deltaw = 1j*(
            -abs(m)*alpe +
            m*I_tot*w0*eta/(2*_np.pi)/E/(nus*w0*bunlen/_LSPEED)**2*Zl_eff)

        if full:
            return deltaw, wp, Zl_interp, Zl_eff
        return deltaw

    def transverse_cbi(
            self, budget=None, element=None, w=None, Zt=None,
            bunlen=None, m=0,  plane='y', inverse=False, full=False):
        """Calcula a impedância transversal efetiva dos nb modos de oscilacao,
        considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        E calcula as instabilidades de Coupled_bunch a partir dela.

        deltaw = transverse_cbi(
            w, Z, bunlen, nb, w0, nus, nut, chrom, eta, m, E, I_tot)
        """
        nus = self.nus
        w0 = self.w0
        eta = self.mom_cmpct
        E = self.E
        nb = self.nbun
        I_tot = self.nom_cur
        taue = self.dampte
        if plane.lower().startswith(('x', 'h')):
            taut, nut, chrom, imp = self.damptx, self.nux, self.chromx, 'Zdx'
        else:
            taut, nut, chrom, imp = self.dampty, self.nuy, self.chromy, 'Zdy'

        if bunlen is None:
            bunlen = self.bunlen

        alpe = 1/taue/w0/nus
        alpt = 1/taut/w0/nus
        w_crom = chrom/eta * w0

        if inverse:
            h = self.calc_spectrum(
                w - w_crom, bunlen, max_rad=0, max_azi=m, only=True)**2
            Z = (alpt + abs(m)*alpe)/(I_tot*w0/(4*_np.pi)/(nus*w0)/E)/h
            return Z

        w, Zt = self._prepare_input_impedance(budget, element, w, Zt, imp)

        # Calculate Effective Impedance
        cbmodes = _np.arange(nb)
        wp = self._get_sampling_ang_freq(w, w0, nb, m, nus, cbmodes, nut)
        wpcro = wp - w_crom

        h = self.calc_spectrum(wpcro, bunlen, max_rad=0, max_azi=m, only=True)
        h *= h
        Zt_interp = self._get_interpolated_impedance(wp, w, Zt)
        Zt_eff = (Zt_interp * h).sum(1)

        # Calculate Coupled_bunch Instability
        # Returns the relative Tune Shift:
        deltaw = -1j*(
            alpt + abs(m)*alpe + I_tot*w0/(4*_np.pi)/(nus*w0)/E * Zt_eff)

        if full:
            return deltaw, wp, Zt_interp, Zt_eff
        return deltaw

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
            res = _np.ones(vetor.shape)
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
        powers = my_pow(sigW, max_azi + 2*max_rad)
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

    def longitudinal_mode_coupling(
            self, budget=None, element=None, w=None, Zl=None,
            max_azi=10, max_rad=12, cbmode=0, use_fokker=True):
        """Calculate the longitudinal mode-coupling eigen-values.

        The eigen-values returned here are normalized by the synchrotron
        frequency, which means they are adimensional. Besides, the
        implementation does not guarantee any specific ordering.

        To get the tune-shifts you must take the real part of the eigen-values.
        To get the growth rate in [1/s] you must multiply the imaginary part
        by the angular synchrotron frequency `nus*w0`.
        If any of the growth rates are larger than 0, It means instability if
        `use_fokker` is `True`, otherwise, a comparison with the damping times
        of the machine must be made after the calculation.

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
            cur_bun
            E
            w0
            nus
            bunlen
            mom_cmpct
            nbun
            dampte

        The problem will be solved for all N values of currents in the vector
        `cur_bun`. I `nus` or `bunlen` are arrays, they must be of the same
        size as `cur_bun`. If the bunch length is an array, the execution time
        will be larger because of the need of recalculating the mode-coupling
        matrix for each current.

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
                eigen-values will be found. Must be lower than `nbun`.
                Defaults to 0.
            use_fokker (bool, optional): Whether or not to include the
                fokker-planck terms described in ref. [2]. If False, the
                Vlasov equation will be considered, like in ref. [1].
                Defaults to True.

        Returns:
            (numpy.ndarray, (N, (2*max_azi+1)*(max_rad+1))): normalized
                eigen-modes of the mode-coupling problem.

        """
        I_b = self.cur_bun
        E = self.E
        w0 = self.w0
        bunlen = self.bunlen
        nus = self.nus
        eta = self.mom_cmpct
        nb = self.nbun
        alpe = 1/self.dampte
        if cbmode >= nb:
            cbmode = 0
            print(
                'Coupled bunch mode greater than number of bunchs.\n',
                'Reseting cbmode to 0.')

        # Calculate the fokker-planck matrix:
        F = self._calc_fokker_planck(max_azi, max_rad, alpe, use_fokker)

        if isinstance(bunlen, (float, _np.float_)):
            bunlen = [bunlen]
        if isinstance(nus, (float, _np.float_)):
            nus = _np.ones(I_b.shape) * nus
        ws = nus * w0

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        # There is an approximation here which is only valid for broad-band
        # impedances. I sample the impedance only at the azimuthal mode m=1.
        wp = self._get_sampling_ang_freq(w, w0, nb, 1, nus[0], [cbmode])

        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)
        Zl_wp = Zl_interp / wp

        delta = _np.zeros([len(I_b), (max_rad+1)*(2*max_azi+1)], dtype=complex)
        # if the bunch length is a vector the mode-coupling matrix must be
        # calculated for each frequency:
        if not len(bunlen) == 1:
            for ii in range(len(I_b)):
                D, M = self._calc_vlasov(
                    Zl_wp, wp, bunlen[ii], max_azi, max_rad)

                # We separated the K value from the defition of M so that M
                # could be current independent in case of a constant bunch
                # length. To understand the equation for K implemented here,
                # please look at eqs. 41 and 43 of ref. [2]. and notice that
                # the definition of M_mlnh in 2.26 of ref. [1] has Z/p instead
                # of Z/wp as we have here:
                xi = eta * I_b[ii]/nus[ii]**2/E
                sig_theta = w0*bunlen[ii]/_LSPEED
                K = xi * w0/(2*_np.pi)/(sig_theta)**2
                # Lastly we need to multiply by `nb` because of our
                # generalization to include n equally spaced bunches.
                K *= nb

                # Please, check eq. 43 of ref. [2]:
                A = D + 1j*K*M + 1j*F/ws[ii]
                delta[ii, :] = _np.linalg.eigvals(A)
        else:
            D, M = self._calc_vlasov(Zl_wp, wp, bunlen[0], max_azi, max_rad)
            for ii in range(len(I_b)):
                xi = eta * I_b[ii]/nus[ii]**2/E
                sig_theta = w0*bunlen[0]/_LSPEED
                K = xi * w0/(2*_np.pi)/(sig_theta)**2
                K *= nb

                A = D + 1j*K*M + 1j*F/ws[ii]
                delta[ii, :] = _np.linalg.eigvals(A)
        return delta

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

        D = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        M = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        spectrum = cls.calc_spectrum(wp, bunlen, max_rad, max_azi)
        for m in range(max_azi+1):
            for k in range(max_rad+1):
                # D is the current independent diagonal matrix:
                D[max_azi+m, k, max_azi+m, k] = m
                D[max_azi-m, k, max_azi-m, k] = -m

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
        D = cls._reshape_coupling_matrix(D)
        M = cls._reshape_coupling_matrix(M)
        return D, M

    def reduced_longitudinal_mode_coupling(
            self, budget=None,  element=None, w=None, Zl=None,
            max_azi=10, max_rad=12, cbmode=0, use_fokker=True):
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
            cur_bun
            E
            w0
            nus
            bunlen
            mom_cmpct
            nbun

        The problem will be solved for all N values of currents in the vector
        `cur_bun`. I `nus` or `bunlen` are arrays, they must be of the same
        size as `cur_bun`. If the bunch length is an array, the execution time
        will be larger because of the need of recalculating the mode-coupling
        matrix for each current.

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
                eigen-values will be found. Must be lower than `nbun`.
                Defaults to 0.

        Returns:
            (numpy.ndarray, (N, max_azi*(max_rad+1))): squared normalized
                eigen-modes of the mode-coupling problem.

        """
        I_b = self.cur_bun
        E = self.E
        w0 = self.w0
        bunlen = self.bunlen
        nus = self.nus
        eta = self.mom_cmpct
        nb = self.nbun
        if cbmode >= nb:
            cbmode = 0
            print(
                'Coupled bunch mode greater than number of bunchs.\n',
                'Reseting cbmode to 0.')

        if isinstance(bunlen, (float, _np.float_)):
            bunlen = [bunlen]
        if isinstance(nus, (float, _np.float_)):
            nus = _np.ones(I_b.shape) * nus

        w, Zl = self._prepare_input_impedance(budget, element, w, Zl, 'Zll')

        # There is an approximation here which is only valid for broad-band
        # impedances. I sample the impedance only at the azimuthal mode m=1.
        wp = self._get_sampling_ang_freq(w, w0, nb, 1, nus[0], [cbmode])

        Zl_interp = self._get_interpolated_impedance(wp, w, Zl)
        Zl_wp = Zl_interp / wp

        delta_sqr = _np.zeros([len(I_b), (max_rad+1)*max_azi], dtype=complex)
        # if the bunch length is a vector the mode-coupling matrix must be
        # calculated for each frequency:
        if not len(bunlen) == 1:
            for ii in range(len(I_b)):
                D, M = self._calc_vlasov_reduced(
                    Zl_wp, wp, bunlen[ii], max_azi, max_rad)

                # We separated the K value from the defition of M so that M
                # could be current independent in case of a constant bunch
                # length. To understand the equation for K implemented here,
                # please look at eqs. 2.49 and 2.52 of ref. [2]. and notice
                # that the definition of M_mlnh in 2.26 of ref. [1] has Z/p
                # instead of Z/wp as we have here:
                xi = eta * I_b[ii]/nus[ii]**2/E
                sig_theta = w0*bunlen[ii]/_LSPEED
                K = xi * w0/_np.pi/(sig_theta)**2
                # Lastly we need to multiply by `nb` because of our
                # generalization to include n equally spaced bunches.
                K *= nb

                # Please, check eq. 2.52 of ref. [1]:
                A = D + 1j*K*M
                delta_sqr[ii, :] = _np.linalg.eigvals(A)
        else:
            D, M = self._calc_vlasov_reduced(
                Zl_wp, wp, bunlen[0], max_azi, max_rad)
            for ii in range(len(I_b)):
                xi = eta * I_b[ii]/nus[ii]**2/E
                sig_theta = w0*bunlen[0]/_LSPEED
                K = xi * w0/_np.pi/(sig_theta)**2
                K *= nb

                A = D + 1j*K*M
                delta_sqr[ii, :] = _np.linalg.eigvals(A)
        return delta_sqr

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

        D = _np.zeros([max_azi, max_rad+1, max_azi, max_rad+1], dtype=complex)
        M = _np.zeros([max_azi, max_rad+1, max_azi, max_rad+1], dtype=complex)
        spectrum = cls.calc_spectrum(wp, bunlen, max_rad, max_azi)
        for m in range(1, max_azi+1):
            for k in range(max_rad+1):
                # D is the current independent diagonal matrix:
                D[m-1, k, m-1, k] = m*m

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
        D = cls._reshape_coupling_matrix(D)
        M = cls._reshape_coupling_matrix(M)
        return D, M

    def transverse_mode_coupling(
            self, budget=None, element=None, w=None, Zt=None, plane='y',
            max_azi=3, max_rad=4, cbmode=0, use_fokker=True):
        """Calculate the transverse mode-coupling eigen-values.

        The eigen-values returned here are normalized by the synchrotron
        frequency, which means they are adimensional. Besides, the
        implementation does not guarantee any specific ordering.

        To get the tune-shifts you must take the real part of the eigen-values.
        To get the growth rate in [1/s] you must multiply the imaginary part
        by the angular synchrotron frequency `nus*w0`.
        If any of the growth rates are larger than 0, It means instability if
        `use_fokker` is `True`, otherwise, a comparison with the damping times
        of the machine must be made after the calculation.

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
            cur_bun
            E
            w0
            nus
            bunlen
            mom_cmpct
            nbun
            dampte
            damptx
            dampty
            nux
            nuy
            chromx
            chromy

        The problem will be solved for all N values of currents in the vector
        `cur_bun`. I `nus` or `bunlen` are arrays, they must be of the same
        size as `cur_bun`. If the bunch length is an array, the execution time
        will be larger because of the need of recalculating the mode-coupling
        matrix for each current.

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
                eigen-values will be found. Must be lower than `nbun`.
                Defaults to 0.
            use_fokker (bool, optional): Whether or not to include the
                fokker-planck terms described in ref. [2]. If False, the
                Vlasov equation will be considered, like in ref. [1].
                Defaults to True.

        Returns:
            (numpy.ndarray, (N, (2*max_azi+1)*(max_rad+1))): normalized
                eigen-modes of the mode-coupling problem.

        """
        I_b = self.cur_bun
        E = self.E
        w0 = self.w0
        bunlen = self.bunlen
        nus = self.nus
        eta = self.mom_cmpct
        nb = self.nbun
        alpe = 1/self.dampte
        if cbmode >= nb:
            cbmode = 0
            print(
                'Coupled bunch mode greater than number of bunchs.\n',
                'Reseting cbmode to 0.')

        if plane.lower().startswith(('x', 'h')):
            alpt, nut, chrom, imp = 1/self.damptx, self.nux, self.chromx, 'Zdx'
        else:
            alpt, nut, chrom, imp = 1/self.dampty, self.nuy, self.chromy, 'Zdy'

        # Calculate the fokker-planck matrix:
        F = self._calc_fokker_planck(max_azi, max_rad, alpe, alpt, use_fokker)

        if isinstance(bunlen, (float, _np.float_)):
            bunlen = [bunlen]
        if isinstance(nus, (float, _np.float_)):
            nus = _np.ones(I_b.shape) * nus
        ws = nus * w0

        w, Zt = self._prepare_input_impedance(budget, element, w, Zt, imp)

        # There is an approximation here which is only valid for broad-band
        # impedances. I sample the impedance only at the azimuthal mode m=1.
        wp = self._get_sampling_ang_freq(w, w0, nb, 1, nus[0], [cbmode], nut)
        w_crom = chrom/eta * w0
        wpcro = wp - w_crom

        Zt_interp = self._get_interpolated_impedance(wp, w, Zt)

        delta = _np.zeros([len(I_b), (max_rad+1)*(2*max_azi+1)], dtype=complex)
        # if the bunch length is a vector the mode-coupling matrix must be
        # calculated for each frequency:
        if not len(bunlen) == 1:
            for ii in range(len(I_b)):
                D, M = self._calc_vlasov_transverse(
                    Zt_interp, wpcro, bunlen[ii], max_azi, max_rad)

                # We separated the K value from the defition of M so that M
                # could be current independent in case of a constant bunch
                # length. To understand the equation for K implemented here,
                # please look at eqs. 2.35 and 2.41 of ref. [1]., or to eq. 37
                # of ref. [2].
                K = I_b[ii] * w0/(4*_np.pi)/ws[ii]/E
                # Lastly we need to multiply by `nb` because of our
                # generalization to include n equally spaced bunches.
                K *= nb

                # This is the complete matrix described in eq. 2.37 of ref.
                # [1] with the aditional fokker-plank terms of eq. 37 of
                # ref. [2] normalized by ws:
                A = D + K*M + 1j*F/ws[ii]
                delta[ii, :] = _np.linalg.eigvals(A)
        else:
            D, M = self._calc_vlasov_transverse(
                Zt_interp, wpcro, bunlen[0], max_azi, max_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0/(4*_np.pi)/ws[ii]/E
                A = D + K*M + 1j*F/ws[ii]
                delta[ii, :] = _np.linalg.eigvals(A)
        return delta

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

        D = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        M = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        spectrum = cls.calc_spectrum(wpcro, bunlen, max_rad, max_azi)
        for m in range(max_azi+1):
            for k in range(max_rad+1):
                # D is the current independent diagonal matrix of eq. 2.37:
                D[max_azi+m, k, max_azi+m, k] = m
                D[max_azi-m, k, max_azi-m, k] = -m

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
        D = cls._reshape_coupling_matrix(D)
        M = cls._reshape_coupling_matrix(M)
        return D, M

    @classmethod
    def _calc_fokker_planck(
            cls, max_azi, max_rad, alpe, alpt=0, use_fokker=True):
        """Calculate Fokker-Planck Matrix.

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

        Args:
            max_azi (int): Maximum azimuthal number before truncation.
            max_rad (int): Maximum radial number before truncation.
            alpe (float): longitudinal damping rate [1/s]
            alpt (int, optional): transverse damping rate [1/s]. Defaults to 0.
            use_fokker (bool, optional): Whether or not to calculate the
                fokker planck terms. If False returns a zero matrix.
                Defaults to True.

        Returns:
            (numpy.ndarray, (
                (2*max_azi+1)x(max_rad+1), (2*max_azi+1)x(max_rad+1))):
                    Square Fokker-Planck matrix that relates a_ml to a_kn.

        """
        F = _np.zeros(
            [1+2*max_azi, max_rad+1, 1+2*max_azi, max_rad+1], dtype=complex)
        if not use_fokker:
            return cls._reshape_coupling_matrix(F)

        for m in range(-max_azi, max_azi+1):
            amm2 = abs(m-2)
            amp2 = abs(m+2)
            am = abs(m)
            for l in range(max_rad+1):
                F[max_azi+m, l, max_azi+m, l] = - alpt - alpe*(abs(m) + 2*l)
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
        if alpha - 2 == beta:
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
    def _get_sampling_ang_freq(w, w0, nb, m, nus, cbmodes, nut=0):
        wp = []
        for cbmode in cbmodes:
            # round towards +infinity
            pmin = _np.ceil((w[0]/w0 - (cbmode + m*nus + nut))/nb)
            # round towards -infinity
            pmax = _np.floor((w[-1]/w0 - (cbmode + m*nus + nut))/nb)

            p = _np.arange(pmin, pmax+1)
            wp.append(nb*p + cbmode + nut + m*nus)
        wp = _np.array(wp).squeeze()
        wp *= w0
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

    def budget_summary(self, budget):
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
                r'[kV/pC]', r'[kV/pC]', r'[kV/pC]',
                r'[kV/pC]', r'$10^{-3}$',  r'$10^{-3}$',  r'$10^{-3}$',
                r'$10^{-3}$', r'$10^{-3}$',  r'$10^{-3}$'],
                index=props),
            latex_name=pd.Series([
                '$\kappa_{Loss}$', '$Z_L/n|_{eff}$', '$P_{Loss}$',
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
                nqy=1e3, ntx=1e3, nty=1e3
                )

        for el in budget:
            values = dict()
            w = el.w

            Zl = el.Zll * el.quantity
            if len(Zl) != 0:
                lossf, pl, zn, *_ = self.loss_factor(w=w, Zl=Zl)
                values['lsf'] = lossf
                values['pls'] = pl
                values['zln'] = zn

            Zd = el.Zdy * el.quantity * el.betay
            if len(Zd) != 0:
                kd, tus, *_ = self.kick_factor(w=w, Z=Zd, Imp='Zdy')
                values['kdy'] = kd
                values['ndy'] = tus

            Zq = el.Zqy * el.quantity * el.betay
            if len(Zq) != 0:
                kd, tus, *_ = self.kick_factor(w=w, Z=Zq, Imp='Zqy')
                values['kqy'] = kd
                values['nqy'] = tus

            if len(Zd) != 0 or len(Zq) != 0:
                values['kty'] = values.get('kdy', 0) + values.get('kqy', 0)
                values['nty'] = values.get('ndy', 0) + values.get('nqy', 0)

            Zd = el.Zdx * el.quantity * el.betax
            if len(Zd) != 0:
                kd, tus, *_ = self.kick_factor(w=w, Z=Zd, Imp='Zdx')
                values['kdx'] = kd
                values['ndx'] = tus

            Zq = el.Zqx * el.quantity * el.betax
            if len(Zq) != 0:
                kd, tus, *_ = self.kick_factor(w=w, Z=Zq, Imp='Zqx')
                values['kqx'] = kd
                values['nqx'] = tus

            if len(Zd) != 0 or len(Zq) != 0:
                values['ktx'] = values.get('kdx', 0) + values.get('kqx', 0)
                values['ntx'] = values.get('ndx', 0) + values.get('nqx', 0)

            for prop in values.keys():
                values[prop] *= convert[prop]

            bud_res[el.name] = pd.Series(values)
        return bud_res

    def kicker_power(
            self, gr, Rshunt=15e3, betak=5, Ab=1e-3, betab=5,
            coupled_mode=None, plane='long'):
        '''Calculate the minimum transverse bunch by bunch feedback power
            necessary to damp coupled-bunch oscillations, assuming perfect
            tunning of the system.
            Formula taken from CAS_Digital_Signal_Processing, pag. 314.

        INPUTS:
          gr     : growth rate you want to damp [Hz].
          Rshunt : shunt resistence of the kicker [Ohm].
          betak  : betatron function at kicker position [m].
          betab  : betatron function at bpm postiion [m].
          Ab     : Initial amplitude of the the oscilation at the BPM position [m].
          coupled_mode : coupled bunch mode of the instability. Not used yet.

        OUTPUT: RF power needed to damp instability [W].
        '''
        if plane.lower().startswith('l'):
            P = (2*_np.pi*self.nus(0)*self.E/self.mom_cmpct*gr*Ab/_LSPEED)**2
        else:
            P = 1/betak * self.E**2 * (self.T0*gr)**2 * (Ab**2/betab)
        return P*2/Rshunt

    def calc_energy_spread_increase(
            self, N=10000000, bunlen=3e-3, wake=1e17, larg=50e-6, Ib=1e-3,
            zz=None):
        """."""
        zz = zz or bunlen
        po = _np.exp(-zz*zz/2/bunlen/bunlen)/_np.sqrt(2*_np.pi)/bunlen*larg

        var1 = po*(1-po)*N

        # radiation kick
        fde = (1 - self.T0/self.dampte/2)
        srde = _np.sqrt(1-fde*fde)*self.espread(Ib)

        # wake kick
        factor = _np.sqrt(var1)/N
        kick = self.T0/self.E * Ib * wake * factor
        # kick = np.sqrt(var1)*1.6e-19/3e9 * Wake

        return _np.sqrt(kick**2/(1-fde*fde))
