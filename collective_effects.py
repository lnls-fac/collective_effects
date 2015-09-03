#!/usr/bin/env python3

import numpy as np

class Element:
    def __init__(self,name=None, folder=None, betax=None,betay=None):
        self.name     = name
        self.folder   = folder
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
        self.f0        = 1/self.T0
        self.w0        = 2*np.pi*self.f0 # revolution angular frequency [Hz]
        self.mom_cmpct = 1.7e-4       # momentum compaction factor
        self.E         = 3e9          # energy [eV]
        self.nuy       = 13.116
        self.nux       = 48.131
        self.harm_num  = 864           # harmonic Number

        I = np.linspace(0,4,num=40)
        self.cur_bun   = I*1e-3
        if self.phase.startswith('commissioning'):
            self.phase       = 'Commissioning'
            self.nom_cur     = 0.100       # total current [A]
            self.nus         = 0.00435    # synchrotron tune
            self.espread     = 7.64e-4 + 0*I
            self.emitx       = 271e-12 + 0*I
            self.emity       = 2.71e-12 + 0*I
            self.damptx      = 17.1e-3
            self.dampty      = 22.7e-3
            self.dampte      = 13.6e-3
            self.en_lost_rad = 456740.6 #eV
        elif self.phase.startswith('phase_1'):
            self.phase       = 'Phase 1 with IDS'
            self.nom_cur     = 0.10         # total current [A]
            self.nus         = 0.00435        # synchrotron tune
            self.espread     = lambda x:1e-2*(9.4e-2+3.80e-2*x-1.83e-2*x**2+4.78e-3*x**3-4.73e-4*x**4)
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
            self.espread     = 1e-02*(0.088704 + 0.015765*I -0.005477*I**2 + 0.0012452*I**3 -0.00011434*I**4)
            self.emitx       = 1e-09*(0.18859  + 0.056781*I -0.015909*I**2 + 0.003445* I**3 -0.00031039*I**4)
            self.emity       = 1e-12*(1.6497   + 1.0422*  I -0.51454 *I**2 + 0.14498*  I**3 -0.015059 * I**4)
            self.damptx      = 10.6e-3
            self.dampty      = 12.5e-3
            self.dampte      =  6.9e-3
            self.en_lost_rad = 829761.9 #eV
        elif self.phase.startswith('phase_2_HC'):
            self.phase       = 'Phase 2 with IDS High Current'
            self.nom_cur     = 0.5        # total current [A]
            self.nus         = 0.00435    # synchrotron tune
            self.espread     = 1e-02*(0.088704 + 0.015765*I -0.005477*I**2 + 0.0012452*I**3 -0.00011434*I**4)
            self.emitx       = 1e-09*(0.18859  + 0.056781*I -0.015909*I**2 + 0.003445* I**3 -0.00031039*I**4)
            self.emity       = 1e-12*(1.6497   + 1.0422*  I -0.51454 *I**2 + 0.14498*  I**3 -0.015059 * I**4)
            self.damptx      = 10.6e-3
            self.dampty      = 12.5e-3
            self.dampte      =  6.9e-3
            self.en_lost_rad = 829761.9 #eV


    def loss_factor(w,Z,sigma,w0,nb):
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


        c = 299792458

        pmin = np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = np.arange(pmin, pmax+1)
        wp = w0*p*nb

        h = exp(-(wp*sigma/c).^2)
        interpol_Z = np.interp(wp,w,Z.real)

        lossf = nb*(w0/2/np.pi)*sum(interpol_Z*h)
        return lossf

    def kick_factor(w,Z,sigma,w0,nb):
        # Calcula o kick factor para nb pacotes com comprimento longitudinal sigma
        # igualmente espacados.

        c = 299792458

        pmin = np.ceil( w[0] /(w0*nb))  # arredonda em direcao a +infinito
        pmax = np.floor(w[-1]/(w0*nb)) # arredonda em direcao a -infinito
        p = np.arange(pmin, pmax+1)

        wp = w0*p*nb

        h = np.exp(-(wp*sigma/c)**2)
        interpol_Z = np.interp(wp,w,Z.imag)

        Zt_eff = sum(interpol_Z.*h)
        return nb*(w0/2/pi)*Zt_eff

    def longitudinal_cbi(w, Zl, sigma, nb, w0, nus, eta, E, I_tot, m):
        # Calcula a impedancia longitudinal efetiva dos nb modos de oscilacao,
        # considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        # E calcula as instabilidades de Coupled_bunch a partir dela.
        c = 299792458
        # Calculate Effective Impedance
        pmin = np.ceil( (w[ 0] -    m*nus*w0    )/(w0*nb)) #arredonda em direcao a +infinito
        pmax = np.floor((w[-1] - (nb-1+m*nus)*w0)/(w0*nb)) #arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + np.arange(0,nb)[:,None] + m*nus)

        h = (wp*sigma/c)**(2*abs(m))*np.exp(-(wp*sigma/c)**2)
        interpol_Z = np.interp(wp, w, Zl)
        Zl_eff = np.diag((interpol_Z/wp) * h)

        deltaw = 1j/(2*np.pi*2**m*np.math.factorial(m - 1)) * I_tot*eta/(E*nus*(sigma/c)**2) * Zl_eff
        return deltaw

    def transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot):
        # Calcula a impedância transversal efetiva dos nb modos de oscilacao,
        # considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
        # E calcula as instabilidades de Coupled_bunch a partir dela.
        #
        # deltaw = transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot)

        c = 299792458

        ## Calculate Effective Impedance
        nucro = nut/eta*chrom
        pmin = np.ceil( (w[0]  -        (m*nus + nut)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = np.floor((w[-1] - (nb-1 + nut + m*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(nb*p[None,:] + np.arange(0,nb)[:,None] + nut + m*nus)
        wpcro = wp - nucro*w0

        h = (wpcro*sigma/c)**(2*abs(m))*np.exp(-(wpcro*sigma/c)**2)
        interpol_Z = np.interp(wp,w,Z)
        Zt_eff = np.diag(interpol_Z * h)

        ## Calculate Coupled_bunch Instabilities

        deltaw = -1j*I_tot/(4*np.pi*2**abs(m)*np.math.factorial(abs(m))) / E * w0 * Zt_eff
        return deltaw

    def longitudinal_mode_coupling(w,Z, params):

        def bunch_spectrum(sigW,azi,rad):
            Inl = (1/np.sqrt(np.math.factorial(abs(azi)+rad) * np.math.factorial(rad)) *
                sigW**(abs(azi)+2*rad) * np.exp(-sigW**2))
            return Inl

        def calc_M(interpol_Z, wp, sigma, n_azi, n_rad):
            sigW = wp*sigma/c/np.sqrt(2)
            A = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1])
            M = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1])
            for k in range(n_rad+1):
                for m in range(n_azi+1):
                    Imk = bunch_spectrum(sigW,m,k)
                    A[n_azi+m, k, n_azi+m, k] =  m
                    A[n_azi-m, k, n_azi-m, k] = -m
                    for n in range(m,n_azi+1):
                        Inl =  bunch_spectrum(sigW,n,k)
                        Mmknl = 1j*(1j)**(m-n)*sum((interpol_Z/wp)*Imk*Inl)
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
                        Imk =  bunch_spectrum(sigW,m,k)
                        for n in range(n_azi+1):
                            Inl =  bunch_spectrum(sigW,n,l)
                            Mmknl = 1j*(1j)**(m-n)*sum((interpol_Z/wp)*Imk*Inl)
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


        n_rad = params.n_rad;
        n_azi = params.n_azi;
        sigma = params.sigma;
        I_b   = params.I;
        E     = params.E;
        w0    = params.w0;
        nus   = params.nus;
        eta   = params.eta;
        nb    = params.nb;
        mu    = params.mu;

        c = 299792458;

        pmin = np.ceil( (w[0] -(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = np.floor((w[-1]-(mu + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu + 1*nus);
        interpol_Z = np.interp(wp,w,Z);


        delta = np.zeros([1 + 2*n_azi + n_rad*(2*n_azi+1), len(I_b)]);
        if len(sigma)~=1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wp, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0*eta/(2*np.pi)/(nus*w0)**2/E/(sigma[ii]/c)**2
                B    = A + K*M
                delta[:,ii] = eig(B)
        else:
            A, M = calc_M(interpol_Z, wp, sigma, n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0*eta/(2*np.pi)/(nus*w0)**2/E/(sigma/c)**2
                B = A + K*M
                delta[:,ii] = np.linalg.eigvals(B)


    def transverse_mode_coupling(w,Z, params):

        def bunch_spectrum(sigW,azi,rad):
            Inl = (1/np.sqrt(np.math.factorial(abs(azi)+rad) * np.math.factorial(rad)) *
                sigW**(abs(azi)+2*rad) * np.exp(-sigW**2))
            return Inl

        def calc_M(interpol_Z, wpcro, sigma, n_azi, n_rad):
            sigW = wpcro*sigma/c/np.sqrt(2)
            A = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1])
            M = np.zeros([1 + 2*n_azi, n_rad+1, 1 + 2*n_azi, n_rad+1])
            for k in range(n_rad+1):
                for m in range(n_azi+1):
                    Imk = bunch_spectrum(sigW,m,k)
                    A(n_azi+m, k, n_azi+m, k) =  m
                    A(n_azi-m, k, n_azi-m, k) = -m
                    for n in range(m,n_azi+1):
                        Inl =  bunch_spectrum(sigW,n,k)
                        Mmknl = -1j*(1j)**(m-n)*sum(interpol_Z*Imk*Inl)
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
                        Imk =  bunch_spectrum(sigW,m,k)
                        for n in range(n_azi+1):
                            Inl =  bunch_spectrum(sigW,n,l)
                            Mmknl = -1j*(1j)**(m-n)*sum(interpol_Z*Imk*Inl)
                            M(n_azi+m, k, n_azi+n, l) =             Mmknl
                            M(n_azi-m, k, n_azi+n, l) =             Mmknl
                            M(n_azi+m, k, n_azi-n, l) =             Mmknl
                            M(n_azi-m, k, n_azi-n, l) =             Mmknl
                            M(n_azi+n, l, n_azi+m, k) = (-1)**(m-n)*Mmknl
                            M(n_azi+n, l, n_azi-m, k) = (-1)**(m-n)*Mmknl
                            M(n_azi-n, l, n_azi+m, k) = (-1)**(m-n)*Mmknl
                            M(n_azi-n, l, n_azi-m, k) = (-1)**(m-n)*Mmknl
            return (A.reshape([(1 + 2*n_azi)*(1+n_rad),(1 + 2*n_azi)*(1+n_rad)]).transpose(),
                    M.reshape([(1 + 2*n_azi)*(1+n_rad),(1 + 2*n_azi)*(1+n_rad)]).transpose())


        n_rad = params.n_rad;
        n_azi = params.n_azi;
        sigma = params.sigma;
        I_b   = params.I;
        E     = params.E;
        w0    = params.w0;
        nus   = params.nus;
        nut   = params.nut;
        chrom = params.chrom;
        eta   = params.eta;
        nb    = params.nb;
        mu    = params.mu;

        c = 299792458;

        nucro = nut/eta*chrom
        pmin = np.ceil( (w[0] -(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a +infinito
        pmax = np.floor((w[-1]-(mu + nut + n_azi*nus)*w0)/(w0*nb)) # arredonda em direcao a -infinito

        p = np.arange(pmin,pmax+1)
        wp = w0*(p*nb + mu +nut + 1*nus)
        interpol_Z = np.interp(wp,w,Z)
        wpcro = wp - nucro*w0


        delta = np.zeros([(1 + 2*n_azi)*(1+n_rad), len(I_b)]);
        if len(sigma)~=1:
            for ii in range(len(I_b)):
                A, M = calc_M(interpol_Z, wpcro, sigma[ii], n_azi, n_rad)
                K    = I_b[ii]*nb*w0/(4*np.pi)/(nus*w0)/E
                B    = A + K*M
                delta[:,ii] = np.linalg.eigvals(B)
        else:
            A, M = calc_M(interpol_Z, wpcro, sigma, n_azi, n_rad)
            for ii in range(len(I_b)):
                K = I_b[ii]*nb*w0/(4*np.pi)/(nus*w0)/E
                B = A + K*M
                delta[:,ii] = np.linalg.eigvals(B)
