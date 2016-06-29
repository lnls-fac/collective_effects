#!/usr/bin/env python3

import os as _os
import sh as _sh
import gzip as _gzip
import pickle as _pickle
import numpy as _np
import matplotlib.pyplot as _plt
from matplotlib import rc as _rc

_jnPth = _os.path.sep.join

class EMSimulData:
    def __init__(self):
        self.circumf = 518.396       # Ring circumference [m]
        self.omega0 = 2*_np.pi*299792458/self.circumf   # Revolution frequency [rad/s]
        self.sigma = 2.5e-3          # Nominal Bunch Length[m]
        self.sigmamax = 15e-3        # Maximum Bunch Length of interest[m]
        self.Iavg = 500e-3           # Beam current [A]
        self.h    = 864              # Harmonic Number

        self.path       = _os.path.abspath('.') # Path where the wake file (from any of the softwares) is
        self.datasource = ''    # CST, ACE3P, GdfidL, ECHOz1 ECHOz2, ...
        self.wake_type  = ''    # DX, DY, QX, QY, LL
        self.symmetry   = False # mirror symmetry for transverse?
        self.bunlen     = 0e-3  # Bunch Length Used in simulation[m]
        self.cutoff     = 2     # multiple of the bunch frequency to calculate impedance

        self.s        = _np.array([],dtype=float) # axis: distance from following to drive bunch [m]
        self.W        = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.freq     = _np.array([],dtype=float) # axis: frequency obtained from FFT [GHz]
        self.Z        = _np.array([],dtype=complex) # Imaginary Part of Longitudinal Impedance [Ohm]
        self.klossW   = _np.array([],dtype=float) # Single-bunch Loss Factor Calculated from Wlong [mV/pC]
        self.klossZ   = _np.array([],dtype=float) # Single-bunch Loss Factor Calculated from ReZlong [mV/pC]
        self.sigmak   = _np.array([],dtype=float) # axis: bunch length for kloss|kick integration [mm]
        self.Ploss    = _np.array([],dtype=float) # Power loss from single-bunch loss factor [W]
        self.Plossvec = _np.array([],dtype=float) # Power loss vector for different sigmas [W]

def prepare_struct_for_load(newdir=None, m=0, sigma=5e-4, rootdir=None, code='ECHO'):

    if not rootdir: rootdir = _os.path.abspath(_os.curdir)

    globdata = EMSimulData()
    globdata.simpar.wakepath = rootdir

    if newdir:
        newdir = _jnPth([rootdir,newdir])
        if not _os.path.isdir(newdir): _os.mkdir(newdir)
    else:
        newdir = rootdir

    globdata.simpar.targetdir = newdir
    globdata.simpar.datasource = code
    globdata.simpar.cutoff = 2
    globdata.simpar.bunlen = sigma
    globdata.simpar.m = m
    globdata.simpar.sym = 1
    globdata.simpar.whichaxis = 'y'
    globdata.simpar.units = 1
    return globdata

def _load_data_ACE3P(simpar):

    m       = simpar.m
    sym     = simpar.sym
    whaxis  = simpar.whichaxis
    wdir    = simpar.wakepath
    tardir  = simpar.targetdir
    sigma   = simpar.sigma
    nsigmas = 5
    headerL = 3
    if wdir.startswith(tardir): cpfile = False
    else: cpfile = True

    wakepath = _jnPth([wdir,'wakefield.out'])
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile: _sh.cp(wakepath, tardir)

    spos = loadres[:,0]
    # I know this is correct for ECHO (2015/08/27):
    if m==0: wake = -loadres[:,1]
    else: wake = loadres[:,1]

    spos = spos - nsigmas* sigma   # Performs displacement over s axis

    return spos, wake

def _load_data_GdfidL(globdata):

    m      = simpar.m
    sym    = simpar.sym
    whaxis = simpar.whichaxis
    wdir   = simpar.wakepath
    tardir = simpar.targetdir
    headerL = 11
    if wdir.startswith(tardir): cpfile = False
    else: cpfile = True


    if m==0:
        wakepath = _jnPth([wdir,'Results-Wq_AT_XY.0001'])
        res = _np.loadtxt(wakepath, skiprows=headerL)
        wake = -res[:,1]
    elif m==1:
        if sym:
            shiftpath = _jnPth([wdir,'Results-Wq_AT_XY.0001'])
            wakepath1 = _jnPth([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0001'])
            wakepath2 = _jnPth([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0002'])
            res1 = _np.loadtxt(wakepath1, skiprows=headerL)
            res2 = _np.loadtxt(wakepath2, skiprows=headerL)
            spos = res1[:,0]
            wabs = (res1[:,1] + res2[:,1])/2
        else:
            shiftpath = _jnPth([wdir,'dxdpl','Results-Wq_AT_XY.0001'])
            wakepath1 = _jnPth([wdir,'d'+whaxis+'dpl','Results-W'+whaxis.upper()+'_AT_XY.0001'])
            wakepath2 = _jnPth([wdir,'d'+whaxis+'dpl','Results-W'+whaxis.upper()+'_AT_XY.0002'])
            wakepath3 = _jnPth([wdir,'d'+whaxis+'dmi','Results-W'+whaxis.upper()+'_AT_XY.0001'])
            wakepath4 = _jnPth([wdir,'d'+whaxis+'dmi','Results-W'+whaxis.upper()+'_AT_XY.0002'])
            res1 = _np.loadtxt(wakepath1, skiprows=headerL)
            res2 = _np.loadtxt(wakepath2, skiprows=headerL)
            res3 = _np.loadtxt(wakepath3, skiprows=headerL)
            res4 = _np.loadtxt(wakepath4, skiprows=headerL)
            spos = loadres1[:,0]
            wabs1 = (res1[:,1] + res2[:,1])/2
            wabs2 = (res3[:,1] + res4[:,1])/2
            wabs = (wabs1 - wabs2)/2
            if cpfile:
                _sh.cp(wakepath1, tardir)
                _sh.cp(wakepath2, tardir)
                _sh.cp(wakepath3, tardir)
                _sh.cp(wakepath4, tardir)

        # obtaining shift value
        with open(shiftpath,'r') as fi:
            for _ in range(0,3):
                loadres5 = fi.readline()
        if cpfile: _sh.cp(shiftpath, tardir)

        coord = textscan(loadres5,' %% subtitle= "W_l (x,y)= ( %f, %f ) [m]"')

        if whaxis.startswith('x'):
            shift = coord[1]
        elif whaxis.startswith('y'):
            shift = coord[2]
        wake = wabs/shift

    elif m==2:
        wakepath1 = _jnPth([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0001'])
        if not sym:
            wakepath2 = _jnPth([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0002'])
        loadres1 = _np.loadtxt(wakepath1, skiprows=headerL)
        if cpfile: _sh.cp(wakepath1, tardir)

        spos = loadres1[:,0]
        wabs = loadres1[:,1]

        if ~sym:
            loadres2 = _np.loadtxt(wakepath2, skiprows=headerL)
            if cpfile: _sh.cp(wakepath2, tardir)
            w2 = loadres2[:,1]
            wabs = (wabs - w2)/2

        #obtaining offset value
        with open(wakepath1,'r') as fi:
            for _ in range(0,3):
                loadres5 = fi.readline()
        if cpfile: _sh.cp(wakepath1, tardir)

        if whaxis.startswith('x'):
            coord = textscan(loadres5,' %% subtitle= "integral d/dx W(z) dz, (x,y)=( %f, %f )"')
            shift = coord[1]
        elif whaxis.startswith('y'):
            coord = textscan(loadres5,' %% subtitle= "integral d/dy W(z) dz, (x,y)=( %f, %f )"')
            shift = coord[2]
        wake = -wabs/shift

    return spos, wake

def _load_data_CST(simpar):

    m      = simpar.m
    whaxis = simpar.whichaxis
    wdir   = simpar.wakepath
    tardir = simpar.targetdir
    headerL = 2
    if wdir.startswith(tardir): cpfile = False
    else: cpfile = True

    wakepath = _jnPth([wdir,'wake.txt'])
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile: _sh.cp(wakepath, tardir)

    spos = loadres[:,0]
    wake = loadres[:,1]
    # Adjust s-axis (rescale or shift)
    spos = spos/1000         # Rescaling mm to m
    if m>0: wake = -wake

    return spos, wake

def _load_data_ECHOz1(simpar):

    m      = simpar.m
    sym    = simpar.sym
    whaxis = simpar.whichaxis
    wdir   = simpar.wakepath
    tardir = simpar.targetdir
    headerL = 0
    if wdir.startswith(tardir): cpfile = False
    else: cpfile = True

    if m==0:
        wakepath = _jnPth([wdir,'wake.dat'])
    elif m > 0:
        if sym: wakepath = _jnPth([wdir,'wakeT.dat'])
        else:   wrt = 'symetry error: wrong set'
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile: _sh.cp(wakepath, tardir)

    spos = loadres[:,0]
    # I know this is correct for ECHO (2015/08/27):
    if m==0: wake = -loadres[:,1]
    else: wake = loadres[:,1]
    # Adjust s-axis (rescale or shift)
    spos = spos/100         # Rescaling cm to m

    return spos, wake

def load_wake(globdata):

    codes = {'ECHOz1': _load_data_ECHOz1,
             'ECHOz2': _load_data_ECHOz2,
             'GdfidL': _load_data_GdfidL,
             'ACE3P' : _load_data_ACE3P,
             'CST'   : _load_data_CST
             }

    spos, wake = codes[dsrc](globdata.simpar)

    # Assign to Structure:
    globdata.results.W = wake
    globdata.results.s = spos
    return globdata

def calc_impedance(globdata):
    # Extracts Needed Variables
    m     = globdata.simpar.m
    sigs  = globdata.simpar.bunlen
    wake  = globdata.results.W * 1e12   # rescale W to [V/C]
    saxis = globdata.results.s
    f0    = globdata.ringpar.omega0/(2*_np.pi)

    c     = 299792458
    sigt  = sigs/c

    # frequency scale (Hz):
    dt = (saxis[-1]-saxis[0])/(saxis.shape[0]-1)/c

    # Modified Hanning window proposed by Caiafa. Not sure of mathematical
    # validity:
    window = _np.hanning(2*wake.shape[0])[wake.shape[0]-1:-1]

    # calculates FFT and frequency:
    #fftt == \int exp(-i*2pi*f*t/n) G(t) dt
    # fftt = _np.fft.fft(wake)
    fftt = _np.fft.fft(wake*window)
    freq = _np.fft.fftfreq(wake.shape[0],d=dt)

    # shift the negative frequencies to the correct position
    fftt = _np.roll(fftt,int(_np.floor(fftt.shape[0]/2)))
    freq = _np.roll(freq,int(_np.floor(freq.shape[0]/2)))
    w    = 2*_np.pi*freq

    # Longitudinal position shift to match center of the bunch with zero z:
    shift = _np.exp(-1j*w*saxis[0]/c)

    # Apply correct scale and shift the spectrum:
    VHat = dt * shift * fftt

    # Deconvolve the Transform with a gaussian bunch:
    Jwlist = _np.exp(-(w*sigt)**2/2)
    Z      = VHat/Jwlist

    #Limits the frequency range according to the bunch length
    wmax  = globdata.simpar.cutoff/sigt
    indcs = abs(w) <= wmax
    if m==0:
        globdata.results.freq = freq[indcs]
        # I have to take the conjugate of the fft because:
        #fftt == \int exp(-i*2pi*f*t/n) G(t) dt
        #while impedance, according to Chao and Ng, is given by:
        #Z == \int exp(i*2pi*f*t/n) G(t) dt
        globdata.results.ReZlong =  Z[indcs].real
        globdata.results.ImZlong = -Z[indcs].imag

        #Calc of Z/n
        indcs2 = _np.logical_and(2*_np.pi*abs(globdata.results.freq) < 20e9,
                                globdata.results.freq != 0)
        globdata.results.naxis = globdata.results.freq[indcs2]/f0
        globdata.results.ImZoN = globdata.results.ImZlong[indcs2]/globdata.results.naxis
    elif m>0:
        globdata.results.freq =   freq[indcs]
        #the Transverse impedance, according to Chao and Ng, is given by:
        #Z == i\int exp(i*2pi*f*t/n) G(t) dt
        globdata.results.ReZt = Z[indcs].imag
        globdata.results.ImZt = Z[indcs].real

    return globdata

def calc_loss_factor(globdata):
    # Extracts and Initialize Needed Variables:
    h    = globdata.ringpar.h
    T0   = 2*_np.pi/globdata.ringpar.omega0
    Iavg = globdata.ringpar.Iavg
    sigs = globdata.simpar.bunlen

    wake  = globdata.results.W
    saxis = globdata.results.s
    freq  = globdata.results.freq
    ReZ   = globdata.results.ReZlong

    c = 299792458
    k       = (freq*2*_np.pi)/c
    ksq     = k**2

    # Calculates klossZ vs. sigma:
    sigmax = globdata.ringpar.sigmamax
    sigmin = globdata.simpar.bunlen
    sigi   = _np.linspace(sigmin,sigmax,num=100)

    kZi = _np.zeros(sigi.shape[0])
    for i in range(sigi.shape[0]):
        rhok   = _np.exp(-ksq*sigi[i]**2)
        kZi[i] = _np.trapz(ReZ*rhok, x=k) * c / (2*_np.pi) * 1e-12
    kZ = kZi[0]

    sigvec = _np.array([2.65, 5.3, 2.65, 4, 10, 10],dtype=float)*1e-3  # bunch length scenarios
    Ivec   = _np.array([500, 500, 10, 110, 110, 500],dtype=float)*1e-3 # current scenarios

    kZvec = _np.zeros(sigvec.shape[0])
    for i in range(sigvec.shape[0]):
        rhok     = _np.exp(-ksq*sigvec[i]**2)
        kZvec[i] = _np.trapz(ReZ*rhok, x=k) * c / (2*_np.pi) * 1e-12
    Plossvec = kZvec * Ivec**2 * T0 * 1e12 / h

    globdata.results.klossZ   = kZi
    globdata.results.sigmak   = sigi
    globdata.results.Plossvec = Plossvec

    # Calculates klossW
    ss    = saxis**2
    rhos  = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-ss/(2*sigs**2))
    kW    = _np.trapz(wake*rhos, x=saxis)
    Ploss = kW * Iavg**2 * T0 * 1e12 / h

    globdata.results.klossW = kW
    globdata.results.Ploss  = Ploss

    # Print loss factor calculated in both ways

    print('klossZ = {0:6.5g} mV/pC'.format(kZ*1000))
    print('klossW = {0:6.5g} mV/pC'.format(kW*1000))
    print('Ploss  = {0:6.5g} W     (for {1:5.4g} mA avg current)'.format(Ploss,Iavg*1000))
    return globdata

def calc_kick_factor(globdata):
    # function  globdata = calc_kick_factor(globdata)
    #This function calculates the kick factor according to methodologies:
    # a. using the long. wake data
    # b. using the long impedance data

    # Extracts and Initialize Needed Variables:
    sigs  = globdata.simpar.bunlen
    wake  = globdata.results.W
    saxis = globdata.results.s
    freq  = globdata.results.freq
    ImZ   = globdata.results.ImZt

    c = 299792458

    sigmasq = sigs**2
    w =(freq*2*_np.pi)
    k = w/c
    ksq = k**2

    # Calculates kickZ vs. sigma:
    sigmax = globdata.ringpar.sigmamax
    sigmin = globdata.simpar.bunlen
    sigi = _np.linspace(sigmin,sigmax,num=100)

    rhok  = _np.exp(-ksq*sigs**2)
    kickZ = _np.trapz(ImZ*rhok,x=k) * c / (2*_np.pi) * 1e-12

    kickZi = _np.zeros(sigi.shape[0])
    for i in range(sigi.shape[0]):
        rhok = _np.exp(-ksq*sigi[i]**2)
        kickZi[i] = _np.trapz(ImZ*rhok,x=k) * c / (2*_np.pi) * 1e-12

    # Calculates kickW:
    ss = saxis**2
    rhos = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-ss/(2*sigmasq))
    kickW = _np.trapz(wake*rhos, x=saxis)

    # Assign results to structure:
    globdata.results.kickZ = kickZi
    globdata.results.sigmak = sigi
    globdata.results.kickW = kickW

    # Print kick factor calculated in both ways:
    print('Kick_Z = {0:6.5g} V/pC/m'.format(kickZ))
    print('Kick_W = {0:6.5g} V/pC/m'.format(kickW))
    return globdata

def plot_results(globdata, mostra=False, salva = True):
    _rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #_rc('font',**{'family':'serif','serif':['Palatino']})
    _rc('text', usetex=True)

    # Data info
    tardir= globdata.simpar.targetdir
    dsrc  = globdata.simpar.datasource
    m     = globdata.simpar.m
    waxis = globdata.simpar.whichaxis
    if waxis.startswith('y'):
        wplane = 'Vertical'
    elif waxis.startswith('x'):
        wplane = 'Horizontal'
    taxis = globdata.simpar.whichaxis

    # Wakepotential
    wake  = globdata.results.W
    spos  = globdata.results.s
    sigs  = globdata.simpar.bunlen

    # Impedance
    if m==0:
        rez = globdata.results.ReZlong
        imz = globdata.results.ImZlong
    elif m>0:
        rez = globdata.results.ReZt/1000
        imz = globdata.results.ImZt/1000
    ImZoN = globdata.results.ImZoN
    f     = globdata.results.freq
    naxis = globdata.results.naxis

    # Loss / Kick Factor
    sigi   = globdata.results.sigmak
    kZi    = globdata.results.klossZ
    kW     = globdata.results.klossW
    kickZi = globdata.results.kickZ
    kickW  = globdata.results.kickW

    #% Tick Position # 0: Plot wakepotential
    #% Short Range
    #========= Plot bunch shape =========
    sbun = _np.linspace(-5*sigs,5*sigs,num=1000) # 5 sigma
    bunchshape = wake.max()*_np.exp(-sbun**2/(2*sigs**2))

    _plt.figure(1)
    _plt.plot(sbun*1000,bunchshape,'b',linewidth=2,label='Bunch Shape')
    _plt.plot(spos*1000,wake,'r',linewidth=2,label='Wakepotential')
    _plt.grid(True)
    _plt.xlabel('s [mm]',fontsize=13)
    fname = 'ShortRange'
    if m==0:
        fname += 'LongitWakePot'
        _plt.title ('Longitudinal Wakepotential ('+dsrc+')',fontsize=13)
        _plt.ylabel('W [V]',fontsize=13)
    elif m==1:
        fname += wplane+'DipWakePot'
        _plt.title (wplane+' Dipole Wakepotential ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$W_{{D_{0:s}}}$ [V/m]'.format(waxis),fontsize=13)
    elif m==2:
        fname += wplane+'QuadWakePot'
        _plt.title (wplane+' Quadrupole Wakepotential ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$W_{{Q_{0:s}}}$ [V/m]'.format(waxis),fontsize=13)
    _plt.xlim([spos[0]*1000, 7000*sigs])
    _plt.ylim([wake.min()*1.1, wake.max()*1.1])
    _plt.legend(loc='best')
    if salva: _plt.savefig(_jnPth((tardir,fname+'.svg')))

    #===== Long Range =====
    _plt.figure(2)
    _plt.plot(spos,wake,'r',linewidth=2)
    _plt.grid(True)
    _plt.xlabel('s [m]',fontsize=13)
    fname = 'LongRange'
    if m==0:
        fname += 'LongitWakePot'
        _plt.title ('Longitudinal Wakepotential ('+dsrc+')',fontsize=13)
        _plt.ylabel('W [V]',fontsize=13)
    elif m==1:
        fname += wplane+'DipWakePot'
        _plt.title (wplane+' Dipole Wakepotential ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$W_{{D_{0:s}}}$ [V/m]'.format(waxis),fontsize=13)
    elif m==2:
        fname += wplane+'QuadWakePot'
        _plt.title (wplane+' Quadrupole Wakepotential ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$W_{{Q_{0:s}}}$ [V/m]'.format(waxis),fontsize=13)
    if salva: _plt.savefig(_jnPth((tardir,fname+'.svg')))

    #=========== Plot Impedance ==========================
    _plt.figure(3)
    _plt.plot(f/1e9,rez,'r',linewidth=2,label='Re')
    _plt.plot(f/1e9,imz,'b--',linewidth=2,label='Im')
    _plt.xlabel('Frequency [GHz]',fontsize=13)
    if m==0:
        fname = 'ImpLongit'
        _plt.title('Longitudinal Impedance ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$Z_{||}$ [$\Omega$]',fontsize=13)
    elif m==1:
        fname = 'ImpDip'+wplane
        _plt.title (wplane+' Dipole Impedance ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$Z_{{D_{0:s}}}$ [k$\Omega$/m]'.format(waxis),fontsize=13)
    elif m==2:
        fname = 'ImpQuad'+wplane
        _plt.title (wplane+' Quadrupole Impedance ('+dsrc+')',fontsize=13)
        _plt.ylabel(r'$Z_{{Q_{0:s}}}$ [k$\Omega$/m]'.format(waxis),fontsize=13)
    _plt.grid(True)
    _plt.legend (loc='best')
    _plt.xlim(_np.array(f[[0,-1]],dtype=float)/1e9)
    if salva: _plt.savefig(_jnPth((tardir,fname+'.svg')))

    #===============Plot Loss/Kick Factor vs. Sigma ======================
    if m==0:
        fname = 'LossFactor'
        _plt.figure(4)
        _plt.plot(sigi * 1e3, kZi * 1e3, 'o',markersize=2,label=r'$K_L^Z$')
        _plt.plot(sigs * 1e3, kW * 1e3, '*',markersize=5,linewidth=2,color=[1, 0, 0],label=r'$K_L^W$')
        _plt.xlabel(r'$\sigma$ [mm]')
        _plt.ylabel(r'$K_L$ [mV/pC]')
        _plt.legend(loc='best')
        _plt.grid(True)
        _plt.annotate(r'$K_L^W = {0:5.2f}$ mV/pC'.format(kW*1e3),xy=(sigs*1.1e3, kW*1e3),fontsize=12)
    elif m > 0:
        fname = 'KickFactor'
        if m==1:
            subind = 'D'
        else:
            subind = 'Q'
        _plt.figure(4)
        _plt.plot(sigi * 1e3, kickZi, 'o',markersize=2,label=r"$\kappa_{0:s}^Z$".format(waxis))
        _plt.plot(sigs * 1e3, kickW, '*',markersize=5,linewidth=2,color=[1, 0, 0],label=r"$\kappa_{0:s}^W$".format(waxis))
        _plt.xlabel(r'$\sigma$ [mm]',fontsize=13)
        _plt.ylabel(r'$\kappa_{0:s}$ [V/pC/m]'.format(waxis),fontsize=13)
        _plt.legend(loc='best')
        _plt.grid(True)
        _plt.annotate(r'$\kappa_{0:s}^W = {1:5.2f}$ V/pC/m'.format(waxis,kickW), xy=(sigs * 1.1e3, kickW), fontsize=13)
    if salva: _plt.savefig(_jnPth((tardir,fname+'.svg')))
    if mostra: _plt.show()

def save_results(globdata):
    filesout = globdata.simpar.targetdir
    dsrc     = globdata.simpar.datasource
    m        = globdata.simpar.m

    if m==0:
        wtype = 'long'
    elif m==1:
        wtype = globdata.simpar.whichaxis + 'dip'
    elif m==2:
        wtype = globdata.simpar.whichaxis + 'quad'

    wake = globdata.results.W
    spos = globdata.results.s

    if m==0:
        rez = globdata.results.ReZlong
        imz = globdata.results.ImZlong
    elif m>0:
        rez = globdata.results.ReZt/1000
        imz = globdata.results.ImZt/1000

    ImZoN = globdata.results.ImZoN
    f     = globdata.results.freq
    naxis = globdata.results.naxis
    sigi  = globdata.results.sigmak

    kZi = globdata.results.klossZ
    kW  = globdata.results.klossW

    kickZi = globdata.results.kickZ
    kickW  = globdata.results.kickW

    Ploss = globdata.results.Ploss
    T0    = 2*_np.pi/globdata.ringpar.omega0

    # Tick Position # 2: Export wakepotential
    _np.savetxt(_jnPth((filesout, 'W'+wtype+dsrc+'.txt')),
                _np.array([spos,wake]).transpose(),fmt=['%30.16g','%30.16g'])

    #% Tick Position # 5: Export Impedance
    _np.savetxt(_jnPth((filesout, 'ReZ'+wtype+dsrc+'.txt')),
                _np.array([f,rez]).transpose(),fmt=['%30.16g','%30.16g'])
    _np.savetxt(_jnPth((filesout, 'ImZ'+wtype+dsrc+'.txt')),
                _np.array([f,imz]).transpose(),fmt=['%30.16g','%30.16g'])

    if m==0:
        _np.savetxt(_jnPth((filesout, 'ImZoN'+wtype+dsrc+'.txt')),
                    _np.array([naxis, ImZoN]).transpose(),fmt=['%30.16g','%30.16g'])

    #% Tick Position # 8: Export Loss Factor vs. Sigma and Loss Info
    if m==0:
        with open(_jnPth((filesout,'Loss info_'+dsrc+'.txt')), 'w') as fi:
            fi.writelines('Loss factor Z = {0:10.6f} mV/pC  \n'.format(kZi[0]*1e3))
            fi.writelines('Loss factor W = {0:10.6f} mV/pC  \n'.format(kW*1e3))
            fi.writelines('Power Loss = {0:10.5f} W \n'.format( Ploss))
            fi.writelines('for I = {0:9.4f} mA  h = {1:5.0f}  T0 = {2:8.4f} ns '.format(
                          globdata.ringpar.Iavg*1e3, globdata.ringpar.h, T0*1e9))

        _np.savetxt(_jnPth((filesout, 'Kloss'+dsrc+'.txt')),
                    _np.array([sigi/1e-3, kZi]).transpose(),fmt=['%12.8g','%12.8g'])
    elif m>0:
        with open(_jnPth((filesout,'Kick info_'+wtype+dsrc+'.txt')), 'w') as fi:
            fi.writelines('Kick Z = {0:10.6f} V/pC/m  \n'.format( kickZi[0]))
            fi.writelines('Kick W = {0:10.6f} V/pC/m  \n'.format(kickW))

        _np.savetxt(_jnPth((filesout, 'K'+wtype+dsrc+'.txt')),
                    _np.array([sigi/1e-3, kickZi]).transpose(),fmt=['%12.8g','%12.8g'])


    with _gzip.open(_jnPth((filesout,'globdata'+wtype+dsrc+'.pickle')), 'wb') as f:
        _pickle.dump(globdata,f,_pickle.HIGHEST_PROTOCOL)

def load_results(filename):
    with _gzip.open(filename,'rb') as fh:
        globdata = _pickle.load(fh)
    return globdata

def analysis_example():

    analysis = '''
    #!/usr/bin/env python3
    import import pycolleff.process_wakes as funcs

    newdir = ''
    m = 1
    bunlen = 0.5e-3
    globdata = funcs.prepare_struct_for_load(newdir, m, bunlen)

    # Load wakepotential result from referred software, rescale and save
    #  txt-file on default format
    globdata = funcs.load_wake(globdata)

    # Calculates Impedance Spectrum from Wakepotential Results
    globdata = funcs.calc_impedance(globdata)

    # Calculates Loss Factor
    if m == 0:
        globdata = funcs.calc_loss_factor(globdata)
    elif m > 0:
        globdata = funcs.calc_kick_factor(globdata)

    # Plot Results
    funcs.plot_results(globdata,mostra=True)

    # Export Results
    funcs.save_results(globdata)
    '''
    print(analysis)
    return None
