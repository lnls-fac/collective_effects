#!/usr/bin/env python3
import numpy as _np
import impedances as _imp

_factorial = _np.math.factorial

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
        self.budget    = _imp.Budget() # impedance Budget

        I = _np.linspace(0,4,num=40)
        self.cur_bun     = I*0.0 # bunch current vector in [A]
        self.nom_cur     = 0.00 # total current [A]
        self.nus         = 0.00 # synchrotron tune
        self.espread     = lambda x:0.0+0*x
        self.sigma       = lambda x:0.0+0*x # bunch length in [m]
        self.emitx       = lambda x:0.0+0*x
        self.emity       = lambda x:0.0+0*x
        self.damptx      = 0.0 # horizontal damping time in [s]
        self.dampty      = 0.0 # vertical damping time in [s]
        self.dampte      = 0.0 # longitudinal damping time in [s]
        self.en_lost_rad = 0.0 # Energy lost per turn in [eV]

    def __setattr__(self,attr,value):
        if attr not in self.__dict__:
            raise AttributeError('Invalid attribute '+attr+' for object of class '+self.__class__+'.')
        if not isinstance(value,type(getattr(self,attr))):
            raise TypeError(attr+' must be '+type(getattr(self,attr))+' not '+ type(value)+'.')
        self.__dict__[attr] = value

    def loss_factor(self,w,Z,sigma):
        """ Calcula o loss factor and effective impedance para nb pacotes com
        comprimento longitudinal sigma igualmente espacados.

        Chamada:
          lossf = loss_factor(w,Z,sigma,w0,nb)

        Inputs:
          w = frequencia angular [rad/s]
          Z = Impedancia longitudinal [Ohm]
          sigma = tamanho longitudinal do feixe [m]
          w0 = frequencia angular de revolucao no anel [rad/s]
          nb = numero de pacotes preenchidos
        """

        w0 = self.w0
        nb = self.nbun

        pmin = _np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = _np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = _np.arange(pmin, pmax+1)
        wp = w0*p*nb

        h = _np.exp(-(wp*sigma/c)**2)
        interpol_Z = _np.interp(wp,w,Z.real)

        lossf = nb*(w0/2/_np.pi)*sum(interpol_Z*h)
        return wp, h, lossf

    def kick_factor(self, w,Z,sigma):
        """Calcula o kick factor para nb pacotes com comprimento longitudinal sigma
        igualmente espacados.
        """

        w0 = self.w0
        nb = self.nbun

        pmin = _np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = _np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = _np.arange(pmin, pmax+1)

        wp = w0*p*nb

        h = _np.exp(-(wp*sigma/c)**2)
        interpol_Z = _np.interp(wp,w,Z.imag)

        Zt_eff = sum(interpol_Z*h)
        return nb*(w0/2/_np.pi)*Zt_eff

    def _prepare_input_impedance(self, budget,element,w,Z,imp='Zl'):
        if budget is not None:
            return budget.w, getattr(budget,imp)
        elif element is not None:
            return element.w, getattr(element,imp)
        elif w is not None and Z is not None:
            return w, Z
        elif self.budget is not None:
            return self.budget.w, getattr(self.budget,imp)
        else:
            raise Exception('Incorrect impedance input.')

    def longitudinal_cbi(self, budget=None, element=None, w=None, Zl=None, sigma=None, m=0):
        """Calcula a impedancia longitudinal efetiva dos nb modos de oscilacao,
        considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        E calcula as instabilidades de Coupled_bunch a partir dela.
        """

        assert m > 0, 'azimuthal mode m must be greater than zero.'
        nus  = self.nus
        w0   = self.w0
        eta  = self.mom_cmpct
        E    = self.E
        nb   = self.nbun
        I_tot= self.nom_cur

        if sigma is None: sigma = self.sigma(I_tot/nb)
        w, Zl = prepare_input_impedance(budget,element,w,Zl,'Zl')

        # Calculate Effective Impedance
        pmin = _np.ceil( (w[ 0] -    m*nus*w0    )/(w0*nb)) #arredonda em direcao a +infinito
        pmax = _np.floor((w[-1] - (nb-1+m*nus)*w0)/(w0*nb)) #arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + _np.arange(0,nb)[:,None] + m*nus)

        h = (wp*sigma/c)**(2*abs(m))*_np.exp(-(wp*sigma/c)**2)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Zl.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Zl.real)
        Zl_eff = (interpol_Z/wp * h).sum(1)

        deltaw = 1j/(2*_np.pi*2**m*_factorial(m - 1)) * I_tot*eta/(E*nus*(sigma/c)**2) * Zl_eff
        return deltaw

    def transverse_cbi(self, budget=None, element=None, w=None, Zt=None, sigma=None, m,  plane='y'):
        """Calcula a impedância transversal efetiva dos nb modos de oscilacao,
        considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        E calcula as instabilidades de Coupled_bunch a partir dela.

        deltaw = transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot)
        """
        nus  = self.nus
        w0   = self.w0
        eta  = self.mom_cmpct
        E    = self.E
        nb   = self.nbun
        I_tot= self.nom_cur
        if plane.lower().startswith(('x','h')):
            nut, chrom, imp   = self.nux, self.chromx, 'Zdh'
        else:
            nut, chrom, imp   = self.nuy, self.chromy, 'Zdv'

        if sigma is None: sigma = self.sigma(I_tot/nb)
        w, Zt = prepare_input_impedance(budget,element,w,Zt,imp)

        ## Calculate Effective Impedance
        nucro = chrom/eta
        pmin = _np.ceil( (w[0]  -        (m*nus + nut)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = _np.floor((w[-1] - (nb-1 + nut + m*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + _np.arange(0,nb)[:,None] + nut + m*nus)
        wpcro = wp - nucro*w0

        h = (wpcro*sigma/c)**(2*abs(m))*_np.exp(-(wpcro*sigma/c)**2)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Zt.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Zt.real)
        Zt_eff = (interpol_Z * h).sum(1)

        ## Calculate Coupled_bunch Instability

        deltaw = -1j*I_tot/(4*_np.pi*2**abs(m)*_factorial(abs(m))) / E * w0 * Zt_eff
        return deltaw

    def calc_spectrum(self,wp,sigma,n_rad=4,n_azi=3):
        def my_pow(vetor,n):
            res = _np.ones(vetor.shape)
            for _ in range(n): res *= vetor
            return res

        sigW = wp*sigma/c/_np.sqrt(2)
        spectrum = dict()
        spect    = dict()
        expo     = _np.exp(-my_pow(sigW,2))
        for azi in range(n_azi+1):
            for rad in range(n_rad+1):
                chave = (abs(azi),rad)
                chave2 = abs(azi)+2*rad
                if chave2 not in spect.keys():
                    spect[chave2] = my_pow(sigW,chave2)
                spectrum[chave] = (1/_np.math.sqrt(float(_factorial(abs(azi)+rad) *
                                           _factorial(rad))) * spect[chave2] * expo)
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
        nus   = self.nus
        eta   = self.mom_cmpct
        nb    = self.nbun

        if sigma is None: sigma = self.sigma(I_b)
        w, Zl = prepare_input_impedance(budget,element,w,Zl,'Zl')

        pmin = _np.ceil( (w[0] -(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = _np.floor((w[-1]-(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu + 1*nus)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Z.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Z.real)

        delta = _np.zeros([len(I_b),(n_rad+1)*(2*n_azi+1)],dtype=complex)
        if not len(sigma)==1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wp, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0*eta/(2*_np.pi)/(nus*w0)**2/E/(sigma[ii]/c)**2
                B    = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        else:
            A, M = calc_M(interpol_Z, wp, sigma, n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0*eta/(2*_np.pi)/(nus*w0)**2/E/(sigma/c)**2
                B = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        return delta

    def transverse_mode_coupling(self, budget=None, element=None, w=None, Zl=None, sigma=None, plane='y', n_azi=3, n_rad=4, mu=0):
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
        nus   = self.nus
        eta   = self.mom_cmpct
        nb    = self.nbun
        if mu >= nb:
            mu = 0
            print('Coupled Bunch Mode greater than Number of Bunchs.\n',
                  'Reseting mu to 0.')
        if plane.lower().startswith(('x','h')):
            nut, chrom, imp = self.nux, self.chromx, 'Zdh'
        else:
            nut, chrom, imp = self.nuy, self.chromy, 'Zdv'

        if sigma is None: sigma = self.sigma(I_b)
        w, Zt = prepare_input_impedance(budget,element,w,Zt,imp)

        nucro = chrom/eta
        pmin = _np.ceil( (w[0] -(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = _np.floor((w[-1]-(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = _np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu +nut + 1*nus)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*_np.interp(wp, w, Z.imag) # imaginary must come first
        interpol_Z +=    _np.interp(wp, w, Z.real)
        wpcro = wp - nucro*w0

        delta = _np.zeros([len(I_b),(n_rad+1)*(2*n_azi+1)],dtype=complex)
        if not len(sigma)==1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wpcro, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0/(4*_np.pi)/(nus*w0)/E
                B    = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        else:
            A, M = calc_M(interpol_Z, wpcro, sigma, n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0/(4*_np.pi)/(nus*w0)/E
                B = A + K*M
                delta[ii,:] = _np.linalg.eigvals(B)
        return delta
