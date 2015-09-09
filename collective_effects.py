#!/usr/bin/env python3

import numpy as np

c = 299792458  # Velocity of light [m/s]

class Element:
    def __init__(self,name=None, folder=None, betax=None,betay=None, quantity=1):
        self.name     = name or ''
        self.folder   = folder or os.path.abspath(os.path.curdir)
        self.quantity = quantity
        self.betax    = betax or 7.2
        self.betay    = betay or 11.0
        self.w        = np.array([],dtype=float)
        self.Zl       = np.array([],dtype=float)
        self.Zdv      = np.array([],dtype=float)
        self.Zdh      = np.array([],dtype=float)
        self.Zqv      = np.array([],dtype=float)
        self.Zqh      = np.array([],dtype=float)
        self.z        = np.array([],dtype=float)
        self.Wl       = np.array([],dtype=float)
        self.Wdv      = np.array([],dtype=float)
        self.Wdh      = np.array([],dtype=float)
        self.Wqv      = np.array([],dtype=float)
        self.Wqh      = np.array([],dtype=float)

class Budget(list):
    def __init__(self, lista=None):
        lista = lista or []
        super().__init__(lista)
        assert all(isinstance(x,Element) for x in self)

    def __setitem__(self,k,v):
        assert isinstance(v,Element)
        super().__setitem__(k,v)


class Ring:
    def __init__(self,phase=None):
        self.phase     = phase or 'phase_1'
        self.version   = 'SI.v10.c01'
        self.circ      = 518.396
        self.T0        = self.circ/299792458
        self.f0        = 1/self.T0  # revolution frequency [Hz]
        self.w0        = 2*np.pi*self.f0 # revolution angular frequency [rad/s]
        self.mom_cmpct = 1.7e-4    # momentum compaction factor
        self.E         = 3e9       # energy [eV]
        self.nuy       = 13.116    # vertical tune
        self.nux       = 48.131    # horizontal tune
        self.chromx    = 0.0       # horizontal chromaticity
        self.chromy    = 0.0       # vertical chromaticity
        self.harm_num  = 864       # harmonic Number
        self.nbun      = 864       # number of bunches filled

        I = np.linspace(0,4,num=40)
        self.cur_bun   = I*1e-3
        if self.phase.startswith('commissioning'):
            self.phase       = 'Commissioning'
            self.nom_cur     = 0.100       # total current [A]
            self.nus         = 0.00435    # synchrotron tune
            self.espread     = lambda x:7.64e-4 +0*x
            self.sigma       = lambda x:3e-3    +0*x
            self.emitx       = lambda x:271e-12 +0*x
            self.emity       = lambda x:2.71e-12+0*x
            self.damptx      = 17.1e-3
            self.dampty      = 22.7e-3
            self.dampte      = 13.6e-3
            self.en_lost_rad = 456740.6 #eV
        elif self.phase.startswith('phase_1'):
            self.phase       = 'Phase 1 with IDS'
            self.nom_cur     = 0.10         # total current [A]
            self.nus         = 0.00435        # synchrotron tune
            self.espread     = lambda x:1e-2*(9.4e-2+3.80e-2*x-1.83e-2*x**2+4.78e-3*x**3-4.73e-4*x**4)
            self.sigma       = lambda x:3e-3    +0*x
            self.emitx       = lambda x:1e-9*(2.3e-1+1.57e-1*x-6.36e-2*x**2+1.60e-2*x**3-1.55e-3*x**4)
            self.emity       = lambda x:1e-12*(2.15 +1.87   *x-8.49e-1*x**2+2.25e-1*x**3-2.25e-3*x**4)
            self.damptx      = 12.4e-3
            self.dampty      = 15.1e-3
            self.dampte      =  8.5e-3
            self.en_lost_rad = 685374.1 #eV
        elif self.phase.startswith('phase_2'):
            self.phase       = 'Phase 2 with IDS'
            self.nom_cur     = 0.35        # total current [A]
            self.nus         = 0.00435    # synchrotron tune
            self.espread     = lambda x:1e-2*(8.87e-2+1.58e-2*x-5.48e-3*x**2+1.25e-3*x**3-1.14e-4*x**4)
            self.sigma       = lambda x:3e-3    +0*x
            self.emitx       = lambda x:1e-9*(1.89e-1+5.61e-2*x-1.59e-2*x**2+3.44e-3*x**3-3.10e-4*x**4)
            self.emity       = lambda x:1e-12*(1.6497+1.04220*x-5.15e-1*x**2+1.45e-1*x**3-1.51e-2*x**4)
            self.damptx      = 10.6e-3
            self.dampty      = 12.5e-3
            self.dampte      =  6.9e-3
            self.en_lost_rad = 829761.9 #eV
        elif self.phase.startswith('phase_2_HC'):
            self.phase       = 'Phase 2 with IDS High Current'
            self.nom_cur     = 0.5        # total current [A]
            self.nus         = 0.00435    # synchrotron tune
            self.espread     = lambda x:1e-2*(8.87e-2+1.58e-2*x-5.48e-3*x**2+1.25e-3*x**3-1.14e-4*x**4)
            self.sigma       = lambda x:3e-3    +0*x
            self.emitx       = lambda x:1e-9*(1.89e-1+5.68e-2*x-1.59e-2*x**2+3.45e-3*x**3-3.10e-4*x**4)
            self.emity       = lambda x:1e-12*(1.6497+1.04220*x-5.15e-1*x**2+1.45e-1*x**3-1.51e-2*x**4)
            self.damptx      = 10.6e-3
            self.dampty      = 12.5e-3
            self.dampte      =  6.9e-3
            self.en_lost_rad = 829761.9 #eV


    def loss_factor(self,w,Z,sigma):
        # Calcula o loss factor and effective impedance para nb pacotes com
        # comprimento longitudinal sigma igualmente espacados.
        #
        # Chamada:
        #   lossf = loss_factor(w,Z,sigma,w0,nb)
        #
        # Inputs:
        #   w = frequencia angular [rad/s]
        #   Z = Impedancia longitudinal [Ohm]
        #   sigma = tamanho longitudinal do feixe [m]
        #   w0 = frequencia angular de revolucao no anel [rad/s]
        #   nb = numero de pacotes preenchidos

        w0 = self.w0
        nb = self.nbun

        pmin = np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = np.arange(pmin, pmax+1)
        wp = w0*p*nb

        h = np.exp(-(wp*sigma/c)**2)
        interpol_Z = np.interp(wp,w,Z.real)

        lossf = nb*(w0/2/np.pi)*sum(interpol_Z*h)
        return wp, h, lossf

    def kick_factor(self, w,Z,sigma):
        # Calcula o kick factor para nb pacotes com comprimento longitudinal sigma
        # igualmente espacados.

        w0 = self.w0
        nb = self.nbun

        pmin = np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = np.arange(pmin, pmax+1)

        wp = w0*p*nb

        h = np.exp(-(wp*sigma/c)**2)
        interpol_Z = np.interp(wp,w,Z.imag)

        Zt_eff = sum(interpol_Z*h)
        return nb*(w0/2/np.pi)*Zt_eff

    def longitudinal_cbi(self, w, Zl, sigma, m):
        # Calcula a impedancia longitudinal efetiva dos nb modos de oscilacao,
        # considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        # E calcula as instabilidades de Coupled_bunch a partir dela.

        assert m > 0, 'azimuthal mode m must be greater than zero.'
        nus  = self.nus
        w0   = self.w0
        eta  = self.mom_cmpct
        E    = self.E
        nb   = self.nbun
        I_tot= self.nom_cur

        # Calculate Effective Impedance
        pmin = np.ceil( (w[ 0] -    m*nus*w0    )/(w0*nb)) #arredonda em direcao a +infinito
        pmax = np.floor((w[-1] - (nb-1+m*nus)*w0)/(w0*nb)) #arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + np.arange(0,nb)[:,None] + m*nus)

        h = (wp*sigma/c)**(2*abs(m))*np.exp(-(wp*sigma/c)**2)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*np.interp(wp, w, Zl.imag) # imaginary must come first
        interpol_Z +=    np.interp(wp, w, Zl.real)
        Zl_eff = (interpol_Z/wp * h).sum(1)

        deltaw = 1j/(2*np.pi*2**m*np.math.factorial(m - 1)) * I_tot*eta/(E*nus*(sigma/c)**2) * Zl_eff
        return deltaw

    def transverse_cbi(self, w,Z, sigma, m,  plane='y'):
        # Calcula a impedância transversal efetiva dos nb modos de oscilacao,
        # considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        # E calcula as instabilidades de Coupled_bunch a partir dela.
        #
        # deltaw = transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot)

        nus  = self.nus
        w0   = self.w0
        eta  = self.mom_cmpct
        E    = self.E
        nb   = self.nbun
        I_tot= self.nom_cur
        if plane.lower().startswith(('x','h')):
            nut, chrom   = self.nux, self.chromx
        else:
            nut, chrom   = self.nuy, self.chromy

        ## Calculate Effective Impedance
        nucro = nut/eta*chrom
        pmin = np.ceil( (w[0]  -        (m*nus + nut)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = np.floor((w[-1] - (nb-1 + nut + m*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + np.arange(0,nb)[:,None] + nut + m*nus)
        wpcro = wp - nucro*w0

        h = (wpcro*sigma/c)**(2*abs(m))*np.exp(-(wpcro*sigma/c)**2)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*np.interp(wp, w, Z.imag) # imaginary must come first
        interpol_Z +=    np.interp(wp, w, Z.real)
        Zt_eff = (interpol_Z * h).sum(1)

        ## Calculate Coupled_bunch Instability

        deltaw = -1j*I_tot/(4*np.pi*2**abs(m)*np.math.factorial(abs(m))) / E * w0 * Zt_eff
        return deltaw

    def longitudinal_mode_coupling(self,w,Z, n_azi=10, n_rad=12,mu=0):

        def calc_M(interpol_Z, wp, sigma, n_azi=7, n_rad=6):
            spectrum = dict()
            def my_pow(vetor,n):
                res = np.ones(vetor.shape)
                for _ in range(n): res *= vetor
                return res
            def bunch_spectrum(wp,sigma,azi,rad):
                sigW = wp*sigma/c/np.sqrt(2)
                chave = (abs(azi),rad)
                if chave not in spectrum.keys():
                    spectrum[chave] = (1/np.sqrt(float(np.math.factorial(abs(azi)+rad) *
                                                       np.math.factorial(rad))) *
                                       my_pow(sigW,(abs(azi)+2*rad)) * np.exp(-my_pow(sigW,2)))
                return spectrum[chave]

            sigW = wp*sigma/c/np.sqrt(2)
            A = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            M = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            for k in range(n_rad+1):
                for m in range(n_azi+1):
                    Imk = bunch_spectrum(wp,sigma,m,k)
                    A[n_azi+m, k, n_azi+m, k] =  m
                    A[n_azi-m, k, n_azi-m, k] = -m
                    for n in range(m,n_azi+1):
                        Inl =  bunch_spectrum(wp,sigma,n,k)
                        Mmknl = 1j*(1j)**(m-n)*((interpol_Z/wp)*Imk*Inl).sum()
                        M[n_azi+m, k, n_azi+n, k] =              m*Mmknl
                        M[n_azi-m, k, n_azi+n, k] =      -       m*Mmknl
                        M[n_azi+m, k, n_azi-n, k] =              m*Mmknl
                        M[n_azi-m, k, n_azi-n, k] =      -       m*Mmknl
                        M[n_azi+n, k, n_azi+m, k] =  (-1)**(m-n)*n*Mmknl
                        M[n_azi+n, k, n_azi-m, k] =  (-1)**(m-n)*n*Mmknl
                        M[n_azi-n, k, n_azi+m, k] = -(-1)**(m-n)*n*Mmknl
                        M[n_azi-n, k, n_azi-m, k] = -(-1)**(m-n)*n*Mmknl
                for l in range(k+1,n_rad+1):
                    for m in range(n_azi+1):
                        Imk =  bunch_spectrum(wp,sigma,m,k)
                        for n in range(n_azi+1):
                            Inl =  bunch_spectrum(wp,sigma,n,l)
                            Mmknl = 1j*(1j)**(m-n)*((interpol_Z/wp)*Imk*Inl).sum()
                            M[n_azi+m, k, n_azi+n, l] =              m*Mmknl
                            M[n_azi-m, k, n_azi+n, l] =      -       m*Mmknl
                            M[n_azi+m, k, n_azi-n, l] =              m*Mmknl
                            M[n_azi-m, k, n_azi-n, l] =      -       m*Mmknl
                            M[n_azi+n, l, n_azi+m, k] =  (-1)**(m-n)*n*Mmknl
                            M[n_azi+n, l, n_azi-m, k] =  (-1)**(m-n)*n*Mmknl
                            M[n_azi-n, l, n_azi+m, k] = -(-1)**(m-n)*n*Mmknl
                            M[n_azi-n, l, n_azi-m, k] = -(-1)**(m-n)*n*Mmknl
            return (A.reshape([(1 + 2*n_azi)*(1+n_rad),(1 + 2*n_azi)*(1+n_rad)]).transpose(),
                    M.reshape([(1 + 2*n_azi)*(1+n_rad),(1 + 2*n_azi)*(1+n_rad)]).transpose())

        I_b   = self.cur_bun
        sigma = self.sigma(I_b)
        E     = self.E
        w0    = self.w0
        nus   = self.nus
        eta   = self.mom_cmpct
        nb    = self.nbun

        pmin = np.ceil( (w[0] -(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = np.floor((w[-1]-(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu + 1*nus)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*np.interp(wp, w, Z.imag) # imaginary must come first
        interpol_Z +=    np.interp(wp, w, Z.real)

        delta = np.zeros([len(I_b),(n_rad+1)*(2*n_azi+1)],dtype=complex)
        if not len(sigma)==1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wp, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0*eta/(2*np.pi)/(nus*w0)**2/E/(sigma[ii]/c)**2
                B    = A + K*M
                delta[ii,:] = eig(B)
        else:
            A, M = calc_M(interpol_Z, wp, sigma, n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0*eta/(2*np.pi)/(nus*w0)**2/E/(sigma/c)**2
                B = A + K*M
                delta[ii,:] = np.linalg.eigvals(B)
        return delta

    def transverse_mode_coupling(self,w,Z, plane='y', n_azi=3, n_rad=4, mu=0):

        def calc_M(interpol_Z, wpcro, sigma, n_azi, n_rad):
            spectrum = dict()
            def my_pow(vetor,n):
                res = np.ones(vetor.shape)
                for _ in range(n): res *= vetor
                return res
            def bunch_spectrum(wp,sigma,azi,rad):
                sigW = wp*sigma/c/np.sqrt(2)
                chave = (abs(azi),rad)
                if chave not in spectrum.keys():
                    spectrum[chave] = (1/np.sqrt(float(np.math.factorial(abs(azi)+rad) *
                                                       np.math.factorial(rad))) *
                                       my_pow(sigW,(abs(azi)+2*rad)) * np.exp(-my_pow(sigW,2)))
                return spectrum[chave]

            A = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            M = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1],dtype=complex)
            for k in range(n_rad+1):
                for m in range(n_azi+1):
                    Imk = bunch_spectrum(wpcro,sigma,m,k)
                    A[n_azi+m, k, n_azi+m, k] =  m
                    A[n_azi-m, k, n_azi-m, k] = -m
                    for n in range(m,n_azi+1):
                        Inl =  bunch_spectrum(wpcro,sigma,n,k)
                        Mmknl = -1j*(1j)**(m-n)*(interpol_Z*Imk*Inl).sum()
                        M[n_azi+m, k, n_azi+n, k] =             Mmknl
                        M[n_azi-m, k, n_azi+n, k] =             Mmknl
                        M[n_azi+m, k, n_azi-n, k] =             Mmknl
                        M[n_azi-m, k, n_azi-n, k] =             Mmknl
                        M[n_azi+n, k, n_azi+m, k] = (-1)**(m-n)*Mmknl
                        M[n_azi+n, k, n_azi-m, k] = (-1)**(m-n)*Mmknl
                        M[n_azi-n, k, n_azi+m, k] = (-1)**(m-n)*Mmknl
                        M[n_azi-n, k, n_azi-m, k] = (-1)**(m-n)*Mmknl
                for l in range(k+1,n_rad+1):
                    for m in range(n_azi+1):
                        Imk =  bunch_spectrum(wpcro,sigma,m,k)
                        for n in range(n_azi+1):
                            Inl =  bunch_spectrum(wpcro,sigma,n,l)
                            Mmknl = -1j*(1j)**(m-n)*(interpol_Z*Imk*Inl).sum()
                            M[n_azi+m, k, n_azi+n, l] =             Mmknl
                            M[n_azi-m, k, n_azi+n, l] =             Mmknl
                            M[n_azi+m, k, n_azi-n, l] =             Mmknl
                            M[n_azi-m, k, n_azi-n, l] =             Mmknl
                            M[n_azi+n, l, n_azi+m, k] = (-1)**(m-n)*Mmknl
                            M[n_azi+n, l, n_azi-m, k] = (-1)**(m-n)*Mmknl
                            M[n_azi-n, l, n_azi+m, k] = (-1)**(m-n)*Mmknl
                            M[n_azi-n, l, n_azi-m, k] = (-1)**(m-n)*Mmknl
            return (A.reshape([(1 + 2*n_azi)*(1+n_rad),(1 + 2*n_azi)*(1+n_rad)]).transpose(),
                    M.reshape([(1 + 2*n_azi)*(1+n_rad),(1 + 2*n_azi)*(1+n_rad)]).transpose())

        I_b   = self.cur_bun
        sigma = self.sigma(I_b)
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
            nut, chrom   = self.nux, self.chromx
        else:
            nut, chrom   = self.nuy, self.chromy

        nucro = nut/eta*chrom
        pmin = np.ceil( (w[0] -(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = np.floor((w[-1]-(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu +nut + 1*nus)
        # Complex interpolation is ill-defined
        interpol_Z  = 1j*np.interp(wp, w, Z.imag) # imaginary must come first
        interpol_Z +=    np.interp(wp, w, Z.real)
        wpcro = wp - nucro*w0

        delta = np.zeros([len(I_b),(n_rad+1)*(2*n_azi+1)],dtype=complex)
        if not len(sigma)==1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wpcro, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0/(4*np.pi)/(nus*w0)/E
                B    = A + K*M
                delta[ii,:] = np.linalg.eigvals(B)
        else:
            A, M = calc_M(interpol_Z, wpcro, sigma, n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0/(4*np.pi)/(nus*w0)/E
                B = A + K*M
                delta[ii,:] = np.linalg.eigvals(B)
        return delta
