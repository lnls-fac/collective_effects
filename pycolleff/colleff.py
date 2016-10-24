#!/usr/bin/env python3
import numpy as _np
import pandas as pd
import mathphys as _mp
# import impedances as _imp

_c   = _mp.constants.light_speed
_factorial = _np.math.factorial

_TYPES = dict(version = str,
          circ        = (float,_np.float_),
          T0          = (float,_np.float_),
          f0          = (float,_np.float_),
          w0          = (float,_np.float_),
          mom_cmpct   = (float,_np.float_),
          E           = (float,_np.float_),
          nuy         = (float,_np.float_),
          nux         = (float,_np.float_),
          chromx      = (float,_np.float_),
          chromy      = (float,_np.float_),
          harm_num    = (int,_np.int_),
          nbun        = (int,_np.int_),
        #   budget      = _imp.Budget,
          cur_bun     = _np.ndarray,
          nom_cur     = (float,_np.float_),
          nus         = (float,_np.float_,type(lambda x:x)),
          espread     = type(lambda x:x),
          sigma       = type(lambda x:x),
          emitx       = type(lambda x:x),
          emity       = type(lambda x:x),
          damptx      = (float,_np.float_),
          dampty      = (float,_np.float_),
          dampte      = (float,_np.float_),
          en_lost_rad = (float,_np.float_))


class Ring:
    def __init__(self):
        self.version   = 'version'
        self.circ      = 0.0  # ring circumference in [m]
        self.T0        = 0.0  # revolution period in [s]
        self.f0        = 0.0  # revolution frequency [Hz]
        self.w0        = 0.0  # revolution angular frequency [rad/s]
        self.mom_cmpct = 0.0  # momentum compaction factor
        self.E         = 0.0  # energy [eV]
        self.nuy       = 0.0  # vertical tune
        self.nux       = 0.0  # horizontal tune
        self.chromx    = 0.0  # horizontal chromaticity
        self.chromy    = 0.0  # vertical chromaticity
        self.harm_num  = 1    # harmonic Number
        self.nbun      = 1    # number of bunches filled
        # self.budget    = _imp.Budget() # impedance Budget

        I = _np.linspace(0,4,num=40)
        self.cur_bun     = I*0.0 # bunch current vector in [A]
        self.nom_cur     = 0.00 # total current [A]
        self._nus         = 0.00 # synchrotron tune
        self._espread     = lambda x:0.0+0*x
        self._sigma       = lambda x:0.0+0*x # bunch length in [m]
        self._emitx       = lambda x:0.0+0*x
        self._emity       = lambda x:0.0+0*x
        self.damptx      = 0.0 # horizontal damping time in [s]
        self.dampty      = 0.0 # vertical damping time in [s]
        self.dampte      = 0.0 # longitudinal damping time in [s]
        self.en_lost_rad = 0.0 # Energy lost per turn in [eV]

    def nus(self,I=None):
        if isinstance(self._nus,(int,_np.int,float,_np.float)):
            return self._nus
        elif isinstance(self._nus,type(lambda x:x)):
            if I is None:
                return self._nus(self.nom_cur/self.nbun)
            else:
                return self._nus(I)
        elif isinstance(self._nus,(_np.ndarray,list,tuple)):
            if I is None:
                return _np.interp(self.nom_cur/self.nbun,fp=self._nus,xp=self.cur_bun)
            else:
                return _np.interp(I,fp=self._nus,xp=self.cur_bun)

    def espread(self,I=None):
        if isinstance(self._espread,(int,_np.int,float,_np.float)):
            return self._espread
        elif isinstance(self._espread,type(lambda x:x)):
            if I is None:
                return self._espread(self.nom_cur/self.nbun)
            else:
                return self._espread(I)
        elif isinstance(self._espread,(_np.ndarray,list,tuple)):
            if I is None:
                return _np.interp(self.nom_cur/self.nbun,fp=self._espread,xp=self.cur_bun)
            else:
                return _np.interp(I,fp=self._esprea,xp=self.cur_bun)

    def sigma(self,I=None):
        if isinstance(self._sigma,(int,_np.int,float,_np.float)):
            return self._sigma
        elif isinstance(self._sigma,type(lambda x:x)):
            if I is None:
                return self._sigma(self.nom_cur/self.nbun)
            else:
                return self._sigma(I)
        elif isinstance(self._sigma,(_np.ndarray,list,tuple)):
            if I is None:
                return _np.interp(self.nom_cur/self.nbun,fp=self._sigma,xp=self.cur_bun)
            else:
                return _np.interp(I,fp=self._sigma,xp=self.cur_bun)

    def emitx(self,I=None):
        if isinstance(self._emitx,(int,_np.int,float,_np.float)):
            return self._emitx
        elif isinstance(self._emitx,type(lambda x:x)):
            if I is None:
                return self._emitx(self.nom_cur/self.nbun)
            else:
                return self._emitx(I)
        elif isinstance(self._emitx,(_np.ndarray,list,tuple)):
            if I is None:
                return _np.interp(self.nom_cur/self.nbun,fp=self._emitx,xp=self.cur_bun)
            else:
                return _np.interp(I,fp=self._emitx,xp=self.cur_bun)

    def emity(self,I=None):
        if isinstance(self._emity,(int,_np.int,float,_np.float)):
            return self._emity
        elif isinstance(self._emity,type(lambda x:x)):
            if I is None:
                return self._emity(self.nom_cur/self.nbun)
            else:
                return self._emity(I)
        elif isinstance(self._emity,(_np.ndarray,list,tuple)):
            if I is None:
                return _np.interp(self.nom_cur/self.nbun,fp=self._emity,xp=self.cur_bun)
            else:
                return _np.interp(I,fp=self._emity,xp=self.cur_bun)

    def __str__(self):
        string = ''
        string += '{0:28s}: {1:^20s}\n'.format('Lattice Version',self.version)
        string += '{0:28s}: {1:^20.3f}\n'.format('Circumference [m]',self.circ)
        string += '{0:28s}: {1:^20.3f}\n'.format('Revolution Period [us]',self.T0*1e6)
        string += '{0:28s}: {1:^20.4f}\n'.format('Revolution Frequency [kHz]',self.f0/1e3)
        string += '{0:28s}: {1:^20.1f}\n'.format('Energy [GeV]',self.E/1e9)
        string += '{0:28s}: {1:^20.2e}\n'.format('Momentum Compaction',self.mom_cmpct)
        string += '{0:28s}: {1:^20d}\n'.format('Harmonic Number',self.harm_num)
        string += '{0:28s}: {1:^20.1f}\n'.format('Current [mA]',self.nom_cur*1e3)
        string += '{0:28s}: {1:^20.3f}\n'.format('Current per Bunch [mA]',self.nom_cur/self.nbun*1e3)
        string += '{0:28s}: {1:^20.3e}\n'.format('Synchrotron Tune',self.nus())
        string += '{0:28s}: {1:>9.3f}/{2:<10.3f}\n'.format('Tunes x/y',self.nux,self.nuy)
        string += '{0:28s}: {1:>6.1f}/{2:^6.1f}/{3:<6.1f}\n'.format('Damping Times x/y/e [ms]',self.damptx*1e3,self.dampty*1e3,self.dampte*1e3)
        string += '{0:28s}: {1:^20.2e}\n'.format('Energy Spread',self.espread())
        string += '{0:28s}: {1:^20.2e}\n'.format('Bunch Length [mm]',self.sigma(self.nom_cur/self.nbun)*1e3)
        return string

    def __getattr__(self,attr):
        if attr in {'nus','sigma','emitx','emity','espread'}:
            attr = '_' + attr
            if isinstance(self.__dict__[attr],(float,_np.float)):
                return self.__dict__[attr]
            elif isinstance(self.__dict__[attr],type(lambda x:x)):
                return self.__dict__[attr](self.nom_cur/self.nbun)
            elif isinstance(self.__dict__[attr],(_np.ndarray,list,tuple)):
                return _np.interp(self.nom_cur/self.nbun,xp=self.__dict__[attr],fp=self.cur_bun)
        else:
            return self.__dict__[attr]

    def __setattr__(self,attr,value):
        if attr in {'nus','sigma','emitx','emity','espread',}:
            if isinstance(value,(_np.ndarray,list,tuple)) and len(value) != len(self.cur_bun):
                raise Exception('Length of input must match length of self.cur_bun.')
            self.__dict__['_'+attr] = value
        else:
            self.__dict__[attr] = value

    def loss_factor(self,budget=None, element=None, w=None, Zl=None, sigma=None):
        """ Calculate the loss factor and effective impedance.

        Inputs:
          budget  = instance of Budget class
                        or
          element = instance of Element class
                        or
          w  = angular frequency [rad/s]
          Zl = Longitudinal Impedance [Ohm]

          (optional) sigma = Longitudinal beamsize [m]

        Outputs:
          lossf  = Loss factor in V/C
          Pl     = Power loss in W
          Zl_eff = Effective Impedance in Ohm
          wp     = vector of angular frequencies where the impedance was sampled
          lossfp = vector of loss factor with contribution of each sampled frequency
        """

        w0 = self.w0
        nb = self.nbun
        sigma = sigma or self.sigma(self.nom_cur/nb)

        w, Zl = self._prepare_input_impedance(budget,element,w,Zl,'Zll')

        pmin = _np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = _np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = _np.arange(pmin, pmax+1)
        wp = w0*p*nb

        # h = _np.exp(-(wp*sigma/_c)**2)
        specs = self.calc_spectrum(wp,sigma,n_rad=0,n_azi=1)
        h = specs[(0,0)]**2
        interpol_Z = _np.interp(wp,w,Zl.real)

        # Loss factor and Power loss:
        lossfp = nb*(w0/2/_np.pi)*interpol_Z*h
        lossf  = sum(lossfp)
        Pl     = lossf * (self.T0*self.nom_cur**2/nb)

        #Effective Impedance:
        interpol_Z = _np.interp(wp,w,Zl.imag)
        h = specs[(0,0)]**2
        Zl_eff =  sum(interpol_Z*h/( (wp+1e-4)/w0 )) / sum(h) #pag 223 K.Y.Ng book

        return lossf, Pl, Zl_eff, wp, lossfp

    def kick_factor(self,budget=None, element=None, w=None, Z=None, sigma=None, Imp='Zdy'):
        """ Calculate the kick factor, tune shift and effective impedance.

        Inputs:
          budget  = instance of Budget class
                        or
          element = instance of Element class
                        or
          w  = angular frequency [rad/s]
          Zl = Longitudinal Impedance [Ohm]

          (optional) sigma = Longitudinal beamsize [m]
          (optional) Imp   = string with type of impedance to consider from element or budget
                             options: Zll(default), Zdy, Zqy, Zdx, Zqx

        Outputs:
          Kick_f = Kick Factor factor in V/C
          Tush   = Tune shift
          Zt_eff = Effective Impedance in Ohm
        """

        w0 = self.w0
        nb = self.nbun
        sigma = sigma or self.sigma(self.nom_cur/nb)

        w, Z = self._prepare_input_impedance(budget,element,w,Z,Imp)

        pmin = _np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = _np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = _np.arange(pmin, pmax+1)
        wp = w0*p*nb

        specs = self.calc_spectrum(wp,sigma,n_rad=0,n_azi=1)
        h = specs[(0,0)]**2
        interpol_Z = _np.interp(wp,w,Z.imag)

        Kick_f = nb*(w0/2/_np.pi)* sum(interpol_Z*h)
        Tush = (self.nom_cur/nb)/(2*self.E) * Kick_f / w0

        h = specs[(1,0)]**2
        Zt_eff = sum(interpol_Z*h) / sum(h) #pag 383 K.Y.Ng Book

        return Kick_f, Tush, Zt_eff

    def _prepare_input_impedance(self, budget,element,w,Z,imp='Zll'):
        if budget is not None:
            return budget.w, getattr(budget,imp)
        elif element is not None:
            return element.w, getattr(element,imp)
        elif w is not None and Z is not None:
            return w, Z
        # elif self.budget is not None:
        #     return self.budget.w, getattr(self.budget,imp)
        else:
            raise Exception('Incorrect impedance input.')

    def longitudinal_cbi(self, budget=None, element=None, w=None, Zl=None, sigma=None, m=0):
        """Calculate the complex coeherent frequency shifts of the beam for all the oscilation modes,
        considering a Gaussian beam and only azimuthal mode k=0;
        """

        assert m > 0, 'azimuthal mode m must be greater than zero.'
        nus  = self.nus()
        w0   = self.w0
        eta  = self.mom_cmpct
        E    = self.E
        nb   = self.nbun
        I_tot= self.nom_cur

        if sigma is None: sigma = self.sigma(I_tot/nb)
        w, Zl = self._prepare_input_impedance(budget,element,w,Zl,'Zl')

        # Calculate Effective Impedance
        pmin = _np.ceil( (w[ 0] -    m*nus*w0    )/(w0*nb)) #arredonda em direcao a +infinito
        pmax = _np.floor((w[-1] - (nb-1+m*nus)*w0)/(w0*nb)) #arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + _np.arange(0,nb)[:,None] + m*nus)

        h = self.calc_spectrum(wp,sigma,n_rad=0,n_azi=m,only=True)**2
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Zl.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Zl.real)
        Zl_eff = (interpol_Z/wp * h).sum(1)

        #Returns the relative Tune Shift:
        deltaw = 1j*abs(m)* I_tot*w0*eta/(2*_np.pi)/(nus*w0)**2/E/(sigma/_c)**2 * Zl_eff
        return deltaw

    def transverse_cbi(self, budget=None, element=None, w=None, Zt=None, sigma=None, m=0,  plane='y'):
        """Calcula a impedância transversal efetiva dos nb modos de oscilacao,
        considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        E calcula as instabilidades de Coupled_bunch a partir dela.

        deltaw = transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot)
        """
        nus  = self.nus()
        w0   = self.w0
        eta  = self.mom_cmpct
        E    = self.E
        nb   = self.nbun
        I_tot= self.nom_cur
        if plane.lower().startswith(('x','h')):
            nut, chrom, imp   = self.nux, self.chromx, 'Zdx'
        else:
            nut, chrom, imp   = self.nuy, self.chromy, 'Zdy'

        if sigma is None: sigma = self.sigma(I_tot/nb)
        w, Zt = self._prepare_input_impedance(budget,element,w,Zt,imp)

        ## Calculate Effective Impedance
        nucro = chrom/eta
        pmin = _np.ceil( (w[0]  -        (m*nus + nut)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = _np.floor((w[-1] - (nb-1 + nut + m*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + _np.arange(0,nb)[:,None] + nut + m*nus)
        wpcro = wp - nucro*w0

        h = self.calc_spectrum(wpcro,sigma,n_rad=0,n_azi=m,only=True)**2
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Zt.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Zt.real)
        Zt_eff = (interpol_Z * h).sum(1)

        ## Calculate Coupled_bunch Instability
        #Returns the relative Tune Shift:
        deltaw = -1j*I_tot*w0/(4*_np.pi)/(nus*w0)/E * Zt_eff
        return deltaw

    def calc_spectrum(self,wp,sigma,n_rad=4,n_azi=3,only=False):
        def my_pow(vetor,n):
            res = _np.ones(vetor.shape)
            for _ in range(n): res *= vetor
            return res

        n_azi0, n_rad0 = 0, 0
        if only: n_azi0, n_rad0 = n_azi, n_rad
        sigW = wp*sigma/_c/_np.sqrt(2)
        spectrum = dict()
        spect    = dict()
        expo     = _np.exp(-my_pow(sigW,2))
        for azi in range(n_azi0, n_azi+1):
            for rad in range(n_rad0, n_rad+1):
                chave = (abs(azi),rad)
                chave2 = abs(azi)+2*rad
                if chave2 not in spect.keys():
                    spect[chave2] = my_pow(sigW,chave2)
                spectrum[chave] = (1/_np.math.sqrt(float(_factorial(abs(azi)+rad) *
                                           _factorial(rad))) * spect[chave2] * expo)
        if only: spectrum = spectrum[chave]
        return spectrum

    def longitudinal_mode_coupling(self, budget=None, element=None, w=None, Zl=None, sigma=None, n_azi=10, n_rad=12,mu=0):

        def calc_M(interpol_Z, wp, sigma, n_azi=7, n_rad=6):
            def fill_M(m,k,n,l,Mmknl):
                M[n_azi+m, k, n_azi+n, l] =              m*Mmknl
                M[n_azi-m, k, n_azi+n, l] =      -       m*Mmknl
                M[n_azi+m, k, n_azi-n, l] =              m*Mmknl
                M[n_azi-m, k, n_azi-n, l] =      -       m*Mmknl
                M[n_azi+n, l, n_azi+m, k] =  (-1)**(m-n)*n*Mmknl
                M[n_azi+n, l, n_azi-m, k] =  (-1)**(m-n)*n*Mmknl
                M[n_azi-n, l, n_azi+m, k] = -(-1)**(m-n)*n*Mmknl
                M[n_azi-n, l, n_azi-m, k] = -(-1)**(m-n)*n*Mmknl

            A = _np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            M = _np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            spectrum = self.calc_spectrum(wp,sigma,n_rad,n_azi)
            for k in range(n_rad+1):
                for m in range(n_azi+1):
                    Imk = spectrum[(abs(m),k)]
                    A[n_azi+m, k, n_azi+m, k] =  m
                    A[n_azi-m, k, n_azi-m, k] = -m
                    for n in range(m,n_azi+1):
                        Inl =  spectrum[(abs(n),k)]
                        Mmknl = 1j*(1j)**(m-n)*_np.dot(interpol_Z/wp,Imk*Inl)
                        fill_M(m,k,n,k,Mmknl)
                for l in range(k+1,n_rad+1):
                    for m in range(n_azi+1):
                        Imk =  spectrum[(abs(m),k)]
                        for n in range(n_azi+1):
                            Inl =  spectrum[(abs(n),l)]
                            Mmknl = 1j*(1j)**(m-n)*_np.dot(interpol_Z/wp,Imk*Inl)
                            fill_M(m,k,n,l,Mmknl)
            return (A.swapaxes(0,3).swapaxes(1,2).reshape([(1+2*n_azi)*(1+n_rad),(1+2*n_azi)*(1+n_rad)]).transpose(),
                    M.swapaxes(0,3).swapaxes(1,2).reshape([(1+2*n_azi)*(1+n_rad),(1+2*n_azi)*(1+n_rad)]).transpose())

        I_b   = self.cur_bun
        E     = self.E
        w0    = self.w0
        nus   = self.nus()
        eta   = self.mom_cmpct
        nb    = self.nbun

        if sigma is None: sigma = self.sigma(I_b)
        if isinstance(sigma,(float,_np.float_)): sigma = [sigma]

        w, Zl = self._prepare_input_impedance(budget,element,w,Zl,'Zl')

        pmin = _np.ceil( (w[0] -(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = _np.floor((w[-1]-(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu + 1*nus)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Zl.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Zl.real)

        delta = _np.zeros([len(I_b),(n_rad+1)*(2*n_azi+1)],dtype=complex)
        if not len(sigma)==1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wp, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0*eta/(2*_np.pi)/(nus*w0)**2/E/(sigma[ii]/_c)**2
                B    = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        else:
            A, M = calc_M(interpol_Z, wp, sigma[0], n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0*eta/(2*_np.pi)/(nus*w0)**2/E/(sigma[0]/_c)**2
                B = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        return delta

    def transverse_mode_coupling(self, budget=None, element=None, w=None, Zt=None, sigma=None, plane='y', n_azi=3, n_rad=4, mu=0):
        def calc_M(interpol_Z, wpcro, sigma, n_azi, n_rad):
            def fill_M(m,k,n,l,Mmknl):
                M[n_azi+m, k, n_azi+n, l] =             Mmknl
                M[n_azi-m, k, n_azi+n, l] =             Mmknl
                M[n_azi+m, k, n_azi-n, l] =             Mmknl
                M[n_azi-m, k, n_azi-n, l] =             Mmknl
                M[n_azi+n, l, n_azi+m, k] = (-1)**(m-n)*Mmknl
                M[n_azi+n, l, n_azi-m, k] = (-1)**(m-n)*Mmknl
                M[n_azi-n, l, n_azi+m, k] = (-1)**(m-n)*Mmknl
                M[n_azi-n, l, n_azi-m, k] = (-1)**(m-n)*Mmknl
            A = _np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            M = _np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            spectrum = self.calc_spectrum(wpcro,sigma,n_rad,n_azi)
            for k in range(n_rad+1):
                for m in range(n_azi+1):
                    Imk = spectrum[(abs(m),k)]
                    A[n_azi+m, k, n_azi+m, k] =  m
                    A[n_azi-m, k, n_azi-m, k] = -m
                    for n in range(m,n_azi+1):
                        Inl =  spectrum[(abs(n),k)]
                        Mmknl = -1j*(1j)**(m-n)*_np.dot(interpol_Z,Imk*Inl)
                        fill_M(m,k,n,k,Mmknl)
                for l in range(k+1,n_rad+1):
                    for m in range(n_azi+1):
                        Imk =  spectrum[(abs(m),k)]
                        for n in range(n_azi+1):
                            Inl =  spectrum[(abs(n),l)]
                            Mmknl = -1j*(1j)**(m-n)*_np.dot(interpol_Z,Imk*Inl)
                            fill_M(m,k,n,l,Mmknl)
            return (A.swapaxes(0,3).swapaxes(1,2).reshape([(1+2*n_azi)*(1+n_rad),(1+2*n_azi)*(1+n_rad)]).transpose(),
                    M.swapaxes(0,3).swapaxes(1,2).reshape([(1+2*n_azi)*(1+n_rad),(1+2*n_azi)*(1+n_rad)]).transpose())

        I_b   = self.cur_bun
        E     = self.E
        w0    = self.w0
        nus   = self.nus()
        eta   = self.mom_cmpct
        nb    = self.nbun
        if mu >= nb:
            mu = 0
            print('Coupled Bunch Mode greater than Number of Bunchs.\n',
                  'Reseting mu to 0.')
        if plane.lower().startswith(('x','h')):
            nut, chrom, imp = self.nux, self.chromx, 'Zdx'
        else:
            nut, chrom, imp = self.nuy, self.chromy, 'Zdy'

        if sigma is None: sigma = self.sigma(I_b)
        if isinstance(sigma,(float,_np.float_)): sigma = [sigma]

        w, Zt = self._prepare_input_impedance(budget,element,w,Zt,imp)

        nucro = chrom/eta
        pmin = _np.ceil( (w[0] -(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = _np.floor((w[-1]-(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu +nut + 1*nus)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Zt.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Zt.real)
        wpcro = wp - nucro*w0

        delta = _np.zeros([len(I_b),(n_rad+1)*(2*n_azi+1)],dtype=complex)
        if not len(sigma)==1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wpcro, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0/(4*_np.pi)/(nus*w0)/E
                B    = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        else:
            A, M = calc_M(interpol_Z, wpcro, sigma[0], n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0/(4*_np.pi)/(nus*w0)/E
                B = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        return delta

    def budget_summary(self,budget=None ):
        props = ['lsf', 'zln', 'pls', 'kdx', 'kdy', 'kqx', 'kqy',
                 'ktx', 'kty', 'ndx', 'ndy', 'nqx', 'nqy', 'ntx', 'nty']

        bud_res = pd.DataFrame(dict(
        name=pd.Series([
            'KLoss', 'Zl/n', 'PLoss', 'Kdx', 'Kdy', 'Kqx', 'Kqx','Kx', 'Ky',
            'TuShdx', 'TuShdy', 'TuShqx', 'TuShqy', 'TuShx', 'TuShy'],
            index=props),
        unit=pd.Series([
            '[mV/pC]', '[mOhm]', '[W]', '[kV/pC]', '[kV/pC]', '[kV/pC]', '[kV/pC]',
            '[kV/pC]', '[kV/pC]', '1/10^3',  '1/10^3',  '1/10^3',  '1/10^3',  '1/10^3',  '1/10^3'],
            index=props),
        latex_unit=pd.Series([
            r'[mV/pC]', r'[m$\Omega$]', r'[W]', r'[kV/pC]', r'[kV/pC]', r'[kV/pC]', r'[kV/pC]', r'[kV/pC]',
            r'[kV/pC]', r'$10^{-3}$',  r'$10^{-3}$',  r'$10^{-3}$',  r'$10^{-3}$',  r'$10^{-3}$',  r'$10^{-3}$'],
            index=props),
        latex_name=pd.Series([
            '$\kappa_{Loss}$', '$Z_L/n|_{eff}$', '$P_{Loss}$',
            r'$\kappa_x^D$', r'$\kappa_y^D$', r'$\kappa_x^Q$', r'$\kappa_y^Q$',r'$\kappa_x$', r'$\kappa_y$',
            r'$\Delta\nu_x^D$', r'$\Delta\nu_y^D$', r'$\Delta\nu_x^Q$', r'$\Delta\nu_y^Q$', r'$\Delta\nu_x$', r'$\Delta\nu_y$'],
            index=props)
        ))

        convert = dict(
                lsf=1e-9, zln=1e3, pls=1, kdx=1e-15, kdy=1e-15, kqx=1e-15, kqy=1e-15,
                ktx=1e-15, kty=1e-15, ndx=1e3, ndy=1e3, nqx=1e3, nqy=1e3, ntx=1e3, nty=1e3
                )

        for el in budget:
            values = dict()
            w = el.w

            Zl = el.Zll * el.quantity
            if len(Zl) != 0:
                values['lsf'], values['pls'],values['zln'],*_ = self.loss_factor(w=w,Zl=Zl)

            Zd = el.Zdy * el.quantity * el.betay
            if len(Zd) != 0:
                values['kdy'], values['ndy'],*_ = self.kick_factor(w=w,Z=Zd,Imp='Zdy')
            Zq = el.Zqy * el.quantity * el.betay
            if len(Zq) != 0:
                values['kqy'], values['nqy'],*_ = self.kick_factor(w=w,Z=Zq,Imp='Zqy')
            if len(Zd) != 0 or len(Zq) != 0:
                values['kty'] = values.get('kdy',0) + values.get('kqy',0)
                values['nty'] = values.get('ndy',0) + values.get('nqy',0)

            Zd = el.Zdx * el.quantity * el.betax
            if len(Zd) != 0:
                values['kdx'], values['ndx'],*_ = self.kick_factor(w=w,Z=Zd,Imp='Zdx')
            Zq = el.Zqx * el.quantity * el.betax
            if len(Zq) != 0:
                values['kqx'], values['nqx'],*_ = self.kick_factor(w=w,Z=Zq,Imp='Zqx')
            if len(Zd) != 0 or len(Zq) != 0:
                values['ktx'] = values.get('kdx',0) + values.get('kqx',0)
                values['ntx'] = values.get('ndx',0) + values.get('nqx',0)

            for prop in values.keys():   values[prop] *= convert[prop]

            bud_res[el.name] = pd.Series(values)
        return bud_res

    def kicker_power(self, gr, Rshunt=15e3, betak=5, Ab=1e-3, betab=5, coupled_mode=None):
        ''' Calculate the minimum transverse bunch by bunch feedback power necessary
        to damp coupled-bunch oscillations, assuming perfect tunning of the system.
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

        P = 2/Rshunt/betak * self.E**2 * (self.T0*gr)**2 * (Ab**2/betab)
        return P
