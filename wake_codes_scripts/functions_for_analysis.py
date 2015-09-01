#!/usr/bin/env python3

import os
import shutil
import gzip
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sirius_structures as si_struc

def prepare_struct_for_load(newdir=None, m=0, sigma=5e-4, rootdir=os.path.abspath(os.curdir)):
    globdata = si_struc.GlobData()

    if newdir is None: raise Exception('must give a folder name')

    globdata.simpar.wakepath = rootdir

    if newdir:
        newdir = os.path.sep.join([rootdir,newdir])
        if not os.path.isdir(newdir): os.mkdir(newdir)
    else:
        newdir = rootdir

    globdata.simpar.targetdir = newdir
    globdata.simpar.datasource = 'ECHO'
    globdata.simpar.cutoff = 2
    globdata.simpar.bunlen = sigma
    globdata.simpar.m = m
    globdata.simpar.sym = 1
    globdata.simpar.whichaxis = 'y'
    globdata.simpar.units = 1
    return globdata


def load_wake(globdata):
    nsigmas = 5 # Only used for ACE3P displacement, standard value!

    dsrc   = globdata.simpar.datasource
    m      = globdata.simpar.m
    sym    = globdata.simpar.sym
    whaxis = globdata.simpar.whichaxis
    wdir   = globdata.simpar.wakepath
    tardir = globdata.simpar.targetdir
    if wdir.startswith(tardir): cpfile = False
    else: cpfile = True

    # Each Software uses a different nomenclature. Setting complete path for loading:
    if dsrc.startswith('ACE3P'):
        wakepath = os.path.sep.join([wdir,'wakefield.out'])
        headerL = 3
    elif dsrc.startswith('GdfidL'):
        headerL = 11

        if m==0:
            wakepath = os.path.sep.join([wdir,'Results-Wq_AT_XY.0001'])
        elif m==1:
            if sym:
                shiftpath = os.path.sep.join([wdir,'Results-Wq_AT_XY.0001'])
                wakepath1 = os.path.sep.join([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0001'])
                wakepath2 = os.path.sep.join([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0002'])
            else:
                shiftpath = os.path.sep.join([wdir,'dxdpl','Results-Wq_AT_XY.0001'])
                wakepath1 = os.path.sep.join([wdir,'d'+whaxis+'dpl','Results-W'+whaxis.upper()+'_AT_XY.0001'])
                wakepath2 = os.path.sep.join([wdir,'d'+whaxis+'dpl','Results-W'+whaxis.upper()+'_AT_XY.0002'])
                wakepath3 = os.path.sep.join([wdir,'d'+whaxis+'dmi','Results-W'+whaxis.upper()+'_AT_XY.0001'])
                wakepath4 = os.path.sep.join([wdir,'d'+whaxis+'dmi','Results-W'+whaxis.upper()+'_AT_XY.0002'])
        elif m==2:
            wakepath1 = os.path.sep.join([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0001'])
            if not sym:
                wakepath2 = os.path.sep.join([wdir,'Results-W'+whaxis.upper()+'_AT_XY.0002'])
    elif dsrc.startswith('CST'):
        wakepath = os.path.sep.join([wdir,'wake.txt'])
        headerL = 2
    elif dsrc.startswith('ECHO'):
        if m==0:
            wakepath = os.path.sep.join([wdir,'wake.dat'])
            headerL = 0
        elif m > 0:
            if sym:
                wakepath = os.path.sep.join([wdir,'wakeT.dat'])
                headerL = 0
            else:
                wrt = 'symetry error: wrong set'

    # Read specified file(s)
    if not dsrc.startswith('GdfidL'):
        loadres = np.loadtxt(wakepath, skiprows=headerL)
        if cpfile: shutil.copy2(wakepath, tardir)

        spos = loadres[:,0]
        # I know this is correct for ECHO (2015/08/27):
        if m==0: wake = -loadres[:,1]
        else: wake = - loadres[:,1]
    else: # I am not sure for GdfidL:
        if m==0:
            wake = -loadres[:,1]
        elif m==1:
            loadres1 = np.loadtxt(wakepath1, skiprows=headerL)
            loadres2 = np.loadtxt(wakepath2, skiprows=headerL)
            if cpfile:
                shutil.copy2(wakepath1, tardir)
                shutil.copy2(wakepath2, tardir)

            spos = loadres1[:,0]
            wabs = (loadres1[:,1]+loadres2[:,1])/2

            if not sym:
                loadres3 = np.loadtxt(wakepath3, skiprows=headerL)
                loadres4 = np.loadtxt(wakepath4, skiprows=headerL)
                if cpfile:
                    shutil.copy2(wakepath3, tardir)
                    shutil.copy2(wakepath4, tardir)

                wabs2 = (loadres3[:,1]+loadres4[:,1])/2
                wabs = (wabs - wabs2)/2

            # obtaining shift value
            with open(shiftpath,'r') as fi:
                for _ in range(0,3):
                    loadres5 = fi.readline()
            if cpfile: shutil.copy2(shiftpath, tardir)

            coord = textscan(loadres5,' %% subtitle= "W_l (x,y)= ( %f, %f ) [m]"')

            if whaxis.startswith('x'):
                shift = coord[1]
            elif whaxis.startswith('y'):
                shift = coord[2]


            wake = wabs/shift
        elif m==2:
            loadres1 = np.loadtxt(wakepath1, skiprows=headerL)
            if cpfile: shutil.copy2(wakepath1, tardir)

            spos = loadres1[:,0]
            wabs = loadres1[:,1]

            if ~sym:
                loadres2 = np.loadtxt(wakepath2, skiprows=headerL)
                if cpfile: shutil.copy2(wakepath2, tardir)

                w2 = loadres2[:,1]
                wabs = (wabs - w2)/2

            #obtaining offset value
            with open(wakepath1,'r') as fi:
                for _ in range(0,3):
                    loadres5 = fi.readline()
            if cpfile: shutil.copy2(wakepath1, tardir)

            if whaxis.startswith('x'):
                coord = textscan(loadres5,' %% subtitle= "integral d/dx W(z) dz, (x,y)=( %f, %f )"')
                shift = coord[1]
            elif whaxis.startswith('y'):
                coord = textscan(loadres5,' %% subtitle= "integral d/dy W(z) dz, (x,y)=( %f, %f )"')
                shift = coord[2]

            wake = -wabs/shift


    # Adjust s-axis (rescale or shift)

    if dsrc.startswith('ACE3P'):
        spos = spos - nsigmas*globdata.simpar.sigma   # Performs displacement over s axis
    elif dsrc.startswith('CST'):
        spos = spos/1000         # Rescaling mm to m
        if m>0:
            wake = -wake
    elif dsrc.startswith('ECHO'):
        spos = spos/100         # Rescaling cm to m


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
    f0    = globdata.ringpar.omega0/(2*np.pi)

    c     = 299792458
    sigt  = sigs/c

    # frequency scale (Hz):
    dt = (saxis[-1]-saxis[0])/(saxis.shape[0]-1)/c

    # Modified Hanning window proposed by Caiafa. Not sure of mathematical
    # validity:
    # window = np.hanning(2*wake.shape[0])[wake.shape[0]-1:-1]
    # fftt = np.fft.rfft(wake*window)

    # calculates FFT and frequency:
    #fftt == \int exp(-i*2pi*f*t/n) G(t) dt
    fftt = np.fft.fft(wake)
    freq = np.fft.fftfreq(wake.shape[0],d=dt)

    # shift the negative frequencies to the correct position
    fftt = np.roll(fftt,int(np.floor(fftt.shape[0]/2)))
    freq = np.roll(freq,int(np.floor(freq.shape[0]/2)))
    w    = 2*np.pi*freq

    # Longitudinal position shift to match center of the bunch with zero z:
    shift = np.exp(-1j*w*saxis[0]/c)

    # Apply correct scale and shift the spectrum:
    VHat = dt * shift * fftt

    # Deconvolve the Transform with a gaussian bunch:
    Jwlist = np.exp(-(w*sigt)**2/2)
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
        indcs2 = np.logical_and(2*np.pi*abs(globdata.results.freq) < 20e9,
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
    T0   = 2*np.pi/globdata.ringpar.omega0
    Iavg = globdata.ringpar.Iavg
    sigs = globdata.simpar.bunlen

    wake  = globdata.results.W
    saxis = globdata.results.s
    freq  = globdata.results.freq
    ReZ   = globdata.results.ReZlong

    c = 299792458
    k       = (freq*2*np.pi)/c
    ksq     = k**2

    # Calculates klossZ vs. sigma:
    sigmax = globdata.ringpar.sigmamax
    sigmin = globdata.simpar.bunlen
    sigi   = np.linspace(sigmin,sigmax,num=100)

    kZi = np.zeros(sigi.shape[0])
    for i in range(sigi.shape[0]):
        rhok   = np.exp(-ksq*sigi[i]**2)
        kZi[i] = np.trapz(ReZ*rhok, x=k) * c / (2*np.pi) * 1e-12
    kZ = kZi[0]

    sigvec = np.array([2.65, 5.3, 2.65, 4, 10, 10],dtype=float)*1e-3  # bunch length scenarios
    Ivec   = np.array([500, 500, 10, 110, 110, 500],dtype=float)*1e-3 # current scenarios

    kZvec = np.zeros(sigvec.shape[0])
    for i in range(sigvec.shape[0]):
        rhok     = np.exp(-ksq*sigvec[i]**2)
        kZvec[i] = np.trapz(ReZ*rhok, x=k) * c / (2*np.pi) * 1e-12
    Plossvec = kZvec * Ivec**2 * T0 * 1e12 / h

    globdata.results.klossZ   = kZi
    globdata.results.sigmak   = sigi
    globdata.results.Plossvec = Plossvec

    # Calculates klossW
    ss    = saxis**2
    rhos  = (1/(sigs*np.sqrt(2*np.pi)))*np.exp(-ss/(2*sigs**2))
    kW    = np.trapz(wake*rhos, x=saxis)
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
    w =(freq*2*np.pi)
    k = w/c
    ksq = k**2

    # Calculates kickZ vs. sigma:
    sigmax = globdata.ringpar.sigmamax
    sigmin = globdata.simpar.bunlen
    sigi = np.linspace(sigmin,sigmax,num=100)

    rhok  = np.exp(-ksq*sigs**2)
    kickZ = np.trapz(ImZ*rhok,x=k) * c / (2*np.pi) * 1e-12

    kickZi = np.zeros(sigi.shape[0])
    for i in range(sigi.shape[0]):
        rhok = np.exp(-ksq*sigi[i]**2)
        kickZi[i] = np.trapz(ImZ*rhok,x=k) * c / (2*np.pi) * 1e-12

    # Calculates kickW:
    ss = saxis**2
    rhos = (1/(sigs*np.sqrt(2*np.pi)))*np.exp(-ss/(2*sigmasq))
    kickW = np.trapz(wake*rhos, x=saxis)

    # Assign results to structure:
    globdata.results.kickZ = kickZi
    globdata.results.sigmak = sigi
    globdata.results.kickW = kickW

    # Print kick factor calculated in both ways:
    print('Kick_Z = {0:6.5g} V/pC/m'.format(kickZ))
    print('Kick_W = {0:6.5g} V/pC/m'.format(kickW))
    return globdata

def plot_results(globdata, mostra=False, salva = True):
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

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
    sbun = np.linspace(-5*sigs,5*sigs,num=1000) # 5 sigma
    bunchshape = wake.max()*np.exp(-sbun**2/(2*sigs**2))

    plt.figure(1)
    plt.plot(sbun*1000,bunchshape,'b',linewidth=2,label='Bunch Shape')
    plt.plot(spos*1000,wake,'r',linewidth=2,label='Wakepotential')
    plt.grid(True)
    plt.xlabel('s [mm]',fontsize=13)
    fname = 'ShortRange'
    if m==0:
        fname += 'LongitWakePot'
        plt.title ('Longitudinal Wakepotential ('+dsrc+')',fontsize=13)
        plt.ylabel('W [V]',fontsize=13)
    elif m==1:
        fname += wplane+'DipWakePot'
        plt.title (wplane+' Dipole Wakepotential ('+dsrc+')',fontsize=13)
        plt.ylabel(r'W_D_{0:s} [V/m]'.format(waxis),fontsize=13)
    elif m==2:
        fname += wplane+'QuadWakePot'
        plt.title (wplane+' Quadrupole Wakepotential ('+dsrc+')',fontsize=13)
        plt.ylabel(r'W_Q_{0:s} [V/m]'.format(waxis),fontsize=13)
    plt.xlim([spos[0]*1000, 7000*sigs])
    plt.ylim([wake.min()*1.1, wake.max()*1.1])
    plt.legend(loc='best')
    if salva: plt.savefig(os.path.sep.join((tardir,fname+'.svg')))
    if mostra: plt.show()

    #===== Long Range =====
    plt.figure(2)
    plt.plot(spos,wake,'r',linewidth=2)
    plt.grid(True)
    plt.xlabel('s [m]',fontsize=13)
    fname = 'LongRange'
    if m==0:
        fname += 'LongitWakePot'
        plt.title ('Longitudinal Wakepotential ('+dsrc+')',fontsize=13)
        plt.ylabel('W [V]',fontsize=13)
    elif m==1:
        fname += wplane+'DipWakePot'
        plt.title (wplane+' Dipole Wakepotential ('+dsrc+')',fontsize=13)
        plt.ylabel(r'W_D_{0:s} [V/m]'.format(waxis),fontsize=13)
    elif m==2:
        fname += wplane+'QuadWakePot'
        plt.title (wplane+' Quadrupole Wakepotential ('+dsrc+')',fontsize=13)
        plt.ylabel(r'W_Q_{0:s} [V/m]'.format(waxis),fontsize=13)
    if salva: plt.savefig(os.path.sep.join((tardir,fname+'.svg')))
    if mostra: plt.show()

    #=========== Plot Impedance ==========================
    plt.figure(3)
    plt.plot(f/1e9,rez,'r',linewidth=2,label='Re')
    plt.plot(f/1e9,imz,'b--',linewidth=2,label='Im')
    plt.xlabel('Frequency [GHz]',fontsize=13)
    if m==0:
        fname = 'ImpLongit'
        plt.title('Longitudinal Impedance ('+dsrc+')',fontsize=13)
        plt.ylabel(r'$\displaystyle Z_{||} [\Omega]$',fontsize=13)
    elif m==1:
        fname = 'ImpDip'+wplane
        plt.title (wplane+' Dipole Impedance ('+dsrc+')',fontsize=13)
        plt.ylabel(r'Z_D_{0:s} [k\Omega/m]'.format(waxis),fontsize=13)
    elif m==2:
        fname = 'ImpQuad'+wplane
        plt.title (wplane+' Quadrupole Impedance ('+dsrc+')',fontsize=13)
        plt.ylabel(r'Z_Q_{0:s} [k\Omega/m]'.format(waxis),fontsize=13)
    plt.grid(True)
    plt.legend (loc='best')
    plt.xlim(np.array(f[[0,-1]],dtype=float)/1e9)
    if salva: plt.savefig(os.path.sep.join((tardir,fname+'.svg')))
    if mostra: plt.show()

    #===============Plot Loss/Kick Factor vs. Sigma ======================
    if m==0:
        fname = 'LossFactor'
        plt.figure(4)
        plt.plot(sigi * 1e3, kZi * 1e3, 'o',markersize=2,label=r'$\displaystyle K_L^Z$')
        plt.plot(sigs * 1e3, kW * 1e3, '*',markersize=5,linewidth=2,color=[1, 0, 0],label=r'$\displaystyle K_L^W$')
        plt.xlabel(r'\sigma [mm]')
        plt.ylabel(r'$\displaystyle K_L [mV/pC]$')
        plt.legend(loc='best')
        plt.grid(True)
        plt.annotate(r'$\displaystyle K_L^W = {0:5.2f} mV/pC$'.format(kW*1e3),xy=(sigs*1.1e3, kW*1e3),fontsize=12)
    elif m > 0:
        fname = 'KickFactor'
        if m==1:
            subind = 'D'
        else:
            subind = 'Q'
        plt.figure(4)
        plt.plot(sigi * 1e3, kickZi, 'o',markersize=2,label=r"$\displaystyle\kappa_{0:s}^Z$".format(waxis))
        plt.plot(sigs * 1e3, kickW, '*',markersize=5,linewidth=2,color=[1, 0, 0],label=r"$\displaystyle\kappa_{0:s}^W$".format(waxis))
        plt.xlabel(r'\sigma [mm]',fontsize=13)
        plt.ylabel(r'$\displaystyle\kappa_{0:s} [V/pC/m]$'.format(waxis),fontsize=13)
        plt.legend(loc='best')
        plt.grid(True)
        plt.annotate(r'$\displaystyle\kappa_{0:s}^W = {1:5.2f} V/pC/m$'.format(waxis,kickW), xy=(sigs * 1.1e3, kickW), fontsize=13)
    if salva: plt.savefig(os.path.sep.join((tardir,fname+'.svg')))
    if mostra: plt.show()


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
    T0    = 2*np.pi/globdata.ringpar.omega0

    # Tick Position # 2: Export wakepotential
    np.savetxt(os.path.sep.join((filesout, 'W'+wtype+dsrc+'.txt')),
                np.array([spos,wake]).transpose(),fmt=['%30.16g','%30.16g'])

    #% Tick Position # 5: Export Impedance
    np.savetxt(os.path.sep.join((filesout, 'ReZ'+wtype+dsrc+'.txt')),
                np.array([f,rez]).transpose(),fmt=['%30.16g','%30.16g'])
    np.savetxt(os.path.sep.join((filesout, 'ImZ'+wtype+dsrc+'.txt')),
                np.array([f,imz]).transpose(),fmt=['%30.16g','%30.16g'])

    if m==0:
        np.savetxt(os.path.sep.join((filesout, 'ImZoN'+wtype+dsrc+'.txt')),
                    np.array([naxis, ImZoN]).transpose(),fmt=['%30.16g','%30.16g'])

    #% Tick Position # 8: Export Loss Factor vs. Sigma and Loss Info
    if m==0:
        with open(os.path.sep.join((filesout,'Loss info_'+dsrc+'.txt')), 'w') as fi:
            fi.writelines('Loss factor Z = {0:10.6f} mV/pC  \n'.format(kZi[0]*1e3))
            fi.writelines('Loss factor W = {0:10.6f} mV/pC  \n'.format(kW*1e3))
            fi.writelines('Power Loss = {0:10.5f} W \n'.format( Ploss))
            fi.writelines('for I = {0:9.4f} mA  h = {1:5.0f}  T0 = {2:8.4f} ns '.format(
                          globdata.ringpar.Iavg*1e3, globdata.ringpar.h, T0*1e9))

        np.savetxt(os.path.sep.join((filesout, 'Kloss'+dsrc+'.txt')),
                    np.array([sigi/1e-3, kZi]).transpose(),fmt=['%12.8g','%12.8g'])
    elif m>0:
        with open(os.path.sep.join((filesout,'Kick info_'+wtype+dsrc+'.txt')), 'w') as fi:
            fi.writelines('Kick Z = {0:10.6f} V/pC/m  \n'.format( kickZi[0]))
            fi.writelines('Kick W = {0:10.6f} V/pC/m  \n'.format(kickW))

        np.savetxt(os.path.sep.join((filesout, 'K'+wtype+dsrc+'.txt')),
                    np.array([sigi/1e-3, kickZi]).transpose(),fmt=['%12.8g','%12.8g'])


    with gzip.open(os.path.sep.join((filesout,'globdata'+wtype+dsrc+'.pickle')), 'wb') as f:
        pickle.dump(globdata,f,pickle.HIGHEST_PROTOCOL)

def load_results(filename):
    with gzip.open(filename,'rb') as fh:
        globdata = pickle.load(fh)
    return globdata
