#!/usr/bin/env python3

import os as _os
import re as _re
import sh as _sh
import gzip as _gzip
import pickle as _pickle
import numpy as _np
import matplotlib.pyplot as _plt
from matplotlib import rc as _rc
from mathphys.constants import light_speed as c
from . import colleff as _colleff
from . import sirius as _sirius

_jnPth = _os.path.sep.join

_si = _sirius.create_ring()

sigvec   = _np.array([2.65, 5.3, 2.65, 4, 10, 10],dtype=float)*1e-3  # bunch length scenarios
Ivec     = _np.array([500, 500, 10, 110, 110, 500],dtype=float)*1e-3 # current scenarios
sigplot  = lambda x:_np.linspace(x,10e-3,num=100)
WAKE_FILENAME_REGEXP = r"[\w-]+W[YXq]{1}_AT_XY.[0-9]{4}"

ANALYSIS_TYPES = {'dx', # horizontal impedance
                   'dy', # vertical impedance
                   'll'  # longitudinal and transverse quadrupolar impedances
                   }

class EMSimulData:
    def __init__(self, path=None, code=None, anal_pl=None, cutoff=None):
        self.path      = path or _os.path.abspath('.')  # Path to the wake files
        self.code      = code      # CST, ACE3P, GdfidL, ECHOz1 ECHOz2, ...
        self.anal_pl   = anal_pl   # dx, dy, ll

        self.bunlen = 0.0                     # Bunch Length Used in simulation[m]
        self.sbun = _np.array([],dtype=float) # positions where the bunch is defined [m]
        self.bun  = _np.array([],dtype=float) # bunch profile used in the simulation [As/m]
        self.s    = _np.array([],dtype=float) # axis: distance from following to drive bunch [m]
        self.Wll  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.Wdx  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.Wdy  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.Wqx  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.Wqy  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.freq = _np.array([],dtype=float) # axis: frequency obtained from FFT [GHz]
        self.Zll  = _np.array([],dtype=complex) # Imaginary Part of Longitudinal Impedance [Ohm]
        self.Zdx  = _np.array([],dtype=complex) # Imaginary Part of Longitudinal Impedance [Ohm]
        self.Zdy  = _np.array([],dtype=complex) # Imaginary Part of Longitudinal Impedance [Ohm]
        self.Zqx  = _np.array([],dtype=complex) # Imaginary Part of Longitudinal Impedance [Ohm]
        self.Zqy  = _np.array([],dtype=complex) # Imaginary Part of Longitudinal Impedance [Ohm]
        self._klossW  = None
        self._kckdxW  = None
        self._kckdyW  = None
        self._kckqxW  = None
        self._kckqyW  = None

    def _kfromW(self,wake):
        T0, sigs, spos = _si.T0, self.bunlen, self.s
        rhos  = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW    = _np.trapz(wake*rhos, x=spos)
        return kW

    def klossW(self):
        if self._klossW: return self._klossW
        kW = self._kfromW(self.Wll)
        self._klossW = kW
        return kW
    def kckdxW(self):
        if self._kckdxW: return self._kckdxW
        kW = self._kfromW(self.Wdx)
        self._kckdxW = kW
        return kW
    def kckdyW(self):
        if self._kckdyW: return self._kckdyW
        kW = self._kfromW(self.Wdy)
        self._kckdyW = kW
        return kW
    def kckqxW(self):
        if self._kckqxW: return self._kckqxW
        kW = self._kfromW(self.Wqx)
        self._kckqxW = kW
        return kW
    def kckqyW(self):
        if self._kckqyW: return self._kckqyW
        kW = self._kfromW(self.Wqy)
        self._kckqyW = kW
        return kW

    def PlossW(self, T0=_si.T0, h=_si.harm_num, Iavg=500e-3):
        kW = self.klossW()
        Ploss = kW * Iavg**2 * T0 * 1e12 / h
        return Ploss

    def klossZ(self,sigma=2.65e-3,n=1):
        si.nbun = n
        klossZ,_ = si.loss_factor(w = self.freq*2*_np.pi, Zl = self.Zll, sigma=sigma)
        return klossZ

    def PlossZ(self,sigma=2.65e-3,n=1):
        ring.nbun = n
        _,PlossZ,*_ = ring.loss_factor(w = self.freq*2*_np.pi, Zl = self.Zll, sigma=sigma)
        return PlossZ


def _load_data_ACE3P(simpar):
    raise NotImplementedError('This function was not tested yet.')
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

def _load_data_CST(simpar):
    raise NotImplementedError('This function was not tested yet.')
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

def _load_data_GdfidL(SimulData,silent=True):

    def _load_dados_info(filename):
        dados, info = [], []
        with open(filename) as fh:
            data = fh.read()
        for line in data.splitlines():
            if not line.startswith((' #',' %',' $')):
                dados.append(line)
            else:
                info.append(line)
        return dados, info

    def _get_charge(info):
        for line in info:
            if line.find('total charge')>=0:
                l = line.split(',')[1]
                charge = float(_re.findall(r'\b[-+]?\d+\.?\d+[eE]?[-+]?\d+\b',l)[0])
                break
        return charge

    def _get_integration_path(info):
        for line in info:
            if line.find('subtitle=')>=0:
                x,y = (float(val) for val in _re.findall(r'\b[-+]?\d+\.?\d+[eE]?[-+]?\d+\b',line))
                break
        return x, y

    def _get_longitudinal_info(path,filelist,pl='ll'):
        if not silent: print('Loading longitunal Wake file:')
        fn = [f for f in filelist if f.find('Wq_AT_XY')>=0]
        if not fn:
            if not silent: print('No longitudinal wake file found. It is needed to have one')
            raise Exception('No longitudinal wake file found. It is needed to have one')
        if len(fn)>1:
            if not silent: print('More than one longitudinal wake file found. It is only allowed 1')
            raise Exception('More than one longitudinal wake file found. It is only allowed 1')
        dados, info = _load_dados_info(_jnPth([path,fn[0]]))
        charge = _get_charge(info)
        if not silent: print('Charge of the driving bunch: {0:7.3g} [pC]'.format(charge*1e12))
        xd, yd = _get_integration_path(info)
        if pl == 'll' and (abs(xd) > 1e-10 or abs(yd) > 1e-10) and not silent:
            print('Driving particle not in the origin. Are you sure this is what you want?')
        elif pl !='ll' and abs(xd) < 1e-10 and abs(yd) < 1e-10 and not silent:
            print('The driving bunch is too close to origin. Are you sure this is what you want?')
        spos,wake = _np.loadtxt(dados,unpack=True) # dados is a list of strings
        a = _np.argmin(_np.diff(spos)) + 1
        sbun   = spos[a:]
        bun    = wake[a:]*charge/_np.trapz(wake[a:],x=sbun) * 1e12 # pC
        wake   = -wake[:a]/charge * 1e-12 # V/pC # minus sign because of convention
        spos   = spos[:a]                # m
        bunlen = -sbun[0]/6            # gdfidl uses a bunch with 6-sigma
        if not silent:
            print('Bunch length of the driving bunch: {0:7.3g} mm'.format(bunlen*1e3))
        return spos, wake, sbun, bun, bunlen, xd, yd

    def _get_transversal_info(path,filelist,pl='qx'):
        stri = 'W{0:s}_AT_XY'.format(pl[1].upper())
        fn = [f for f in f_match if f.find(stri)>=0]
        if not fn and not silent:
            print('No W{0:s} wake file found. Skipping to next'.format(pl))
            return None
        if not silent: print('{0:2d} W{1:s} wake file found.'.format(len(fn),pl))
        dados, info = _load_dados_info(_jnPth([path,fn[0]]))
        charge = _get_charge(info)
        if pl[1] == 'x':
            delta1,_ = _get_integration_path(info)
        else:
            _,delta1 = _get_integration_path(info)
        _, wake1 = _np.loadtxt(dados,unpack=True)
        wake = wake1/delta1 / charge * 1e-12 # V/pC/m
        if len(fn) > 1:
            dados, info = _load_dados_info(_jnPth([path,fn[1]]))
            if pl[1] == 'x':
                delta2,_ = _get_integration_path(info)
            else:
                _,delta2 = _get_integration_path(info)
            _, wake2 = _np.loadtxt(dados,unpack=True)
            if pl[0] == 'd':
                wake = (wake1/delta1 - wake2/delta2)/(1/delta1-1/delta2) / charge * 1e-12 # V/pC
            else:
                wake = (wake1 - wake2)/(delta1-delta2) / charge * 1e-12 # V/pC/m
        return wake

    path     = SimulData.path
    anal_pl  = SimulData.anal_pl

    # list all the files that match the name pattern for wakefields
    f_in_dir = _sh.ls(path).stdout.decode()
    f_match = _re.findall(WAKE_FILENAME_REGEXP,f_in_dir)

    if anal_pl in {'ll'}:
        if not f_match:
            if not silent: print('No files found for longitudinal analysis.')
            raise Exception('No files found for longitudinal analisys')

        #Load longitudinal Wake
        spos, wake, sbun, bun, bunlen, xd, yd = _get_longitudinal_info(path,f_match,pl='ll')
        SimulData.Wll  = wake
        SimulData.s    = spos
        SimulData.bun  = bun
        SimulData.sbun = sbun
        SimulData.bunlen = bunlen

        # And quadrupolar Wakes, if existent:
        if not silent: print('Loading Horizontal Quadrupolar Wake file:')
        wake = _get_transversal_info(path,f_match,pl='qx') # V / pC / m
        if wake: SimulData.Wqx = wake
        if not silent: print('Loading Vertical Quadrupolar Wake file:')
        wake = _get_transversal_info(path,f_match,pl='qy') # V / pC / m
        if wake: SimulData.Wqy = wake
        if not silent: print('Longitudinal Data Loaded.\n')

    elif anal_pl in {'dx','dy'}:
        if not f_match:
            if not silent: print('There is no wake files in this folder. I will assume there is no symmetry.')
            spos,wake,sbun,bun,bunlen,xd,yd = [],[],[],[],[],[],[]
            for sub_fol in ['dpl','dmi']:
                ext_path = _jnPth([path,anal_pl+sul_folder])
                if not silent: print('Looking for '+anal_pl+sub_fol+' subfolder:')
                if not _os.path.isdir(ext_path):
                    if not silent: print('For non-symmetric structures, there must '
                                   'be subfolders {0:s}dpl {0:s}dmi with the data'.format(anal_pl))
                    raise Exception('Files not found')
                # list all the files that match the pattern
                f_in_dir = _sh.ls(ext_path).stdout.decode()
                f_match = _re.findall(WAKE_FILENAME_REGEXP,f_in_dir)
                if not f_match:
                    if not silent: print('No files found for transverse analysis.')
                    raise Exception('No files found for transverse analisys')

                sp, _, sb, bn, bnln, xdi, ydi = _get_longitudinal_info(ext_path,f_match,pl=anal_pl)
                spos.append(sp)
                bun.append(bn)
                sbun.append(sb)
                bunlen.append(bnln)
                xd.append(xdi)
                yd.append(ydi)

                if not silent:
                    print('Loading {0:s} Dipolar Wake file:'.format(
                          'Horizontal' if anal_pl=='dx' else 'Vertical'))
                wk = _get_transversal_info(ext_path,f_match,pl=anal_pl) # V / pC
                if wk:
                    wake.append(wk)
                else:
                    if not silent: print('Actually there is something wrong, these wake files should be here.')
                    raise Exception('Transverse {0:s} dipolar wake files not found'.format(
                                    'Horizontal' if anal_pl=='dx' else 'Vertical'))

            delta = xd if anal_pl=='dx' else yd
            ndel  = yd if anal_pl=='dx' else xd
            if not (_np.allclose(spos[0],spos[1],atol=0) and
                    _np.allclose(sbun[0],sbun[1],atol=0) and
                    _np.allclose( bun[0], bun[1],atol=0) and
                    _np.allclose(ndel[0],ndel[1],atol=0)):
                if not silent: print('There is a mismatch between the paramenters of the'
                            'simulation in the {0:s}dpl and {0:s}dmi folders.'.format(anal_pl))
                raise Exception('Mismatch of the parameters of the simulation in the subfolders.')
            SimulData.s      = spos[0]
            SimulData.bun    = bun[0]
            SimulData.sbun   = sbun[0]
            SimulData.bunlen = bunlen[0]
            setattr(SimulData,'W'+anal_pl, (wake[0]-wake[1])/(delta[0]-delta[1])) # V / pC /m
        else:
            spos, wake, sbun, bun, bunlen, xd, yd = _get_longitudinal_info(path,f_match,pl=anal_pl)
            SimulData.s      = spos
            SimulData.bun    = bun
            SimulData.sbun   = sbun
            SimulData.bunlen = bunlen

            if not silent:
                print('Loading {0:s} Dipolar Wake file:'.format(
                      'Horizontal' if anal_pl=='dx' else 'Vertical'))
            wake = _get_transversal_info(path,f_match,pl=anal_pl) # V / pC
            if wake:
                delta = xd if anal_pl=='dx' else yd
                setattr(SimulData,'W'+anal_pl,wake/delta) # V / pC / m
            else:
                print('Actually there is something wrong, these wake files should be here.')
                raise Exception('Transverse {0:s} dipolar wake files not found'.format(
                                'Horizontal' if anal_pl=='dx' else 'Vertical'))
        if not silent: print('Transverse Data Loaded.\n')

def _load_data_ECHOz1(SimulData,silent=True):

    path     = SimulData.path
    anal_pl  = SimulData.anal_pl

    if anal_pl=='ll':
        if not silent: print('Loading longitudinal Wake file:',end='')
        fname = _jnPth([path,'wake.dat'])
        if _os.path.isfile(fname):
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found')
            Exception('Longitudinal wake file not found')
    else:
        if not silent: print('ECHOz1 only calculates longitudinal wake.')
        raise Exception('ECHOz1 only calculates longitudinal wake.')

    SimulData.spos = loadres[:,0]/100    # Rescaling cm to m
    # I know this is correct for ECHO (2015/08/27):
    SimulData.Wll = -loadres[:,1]       # V / pC / m   the minux sign is due to convention

    # loading driving bunch info
    loadres = _np.loadtxt(_jnPth([path,'bunch.dat']), skiprows=0)
    sbun = loadres[:,1] / 100     # m
    a = _np.argmin(_np.abs(sbun + sbun[0])) # I want the bunch to be symmetric
    SimulData.sbun = sbun[:a]
    SimulData.bun  = loadres[:a,2] # pC
    SimulData.bunlen = abs(sbun[0] + sbun[1])/ 2 / 5 # ECHO uses 5 sigma
    if not silent:
        print('Bunch length of the driving bunch: {0:7.3g} mm'.format(SimulData.bunlen*1e3))
        print('Data Loaded.\n')

def _load_data_ECHOz2(SimulData,silent=True):

    path     = SimulData.path
    anal_pl  = SimulData.anal_pl

    if anal_pl=='ll':
        if not silent: print('Loading longitudinal Wake file:',end='')
        fname = _jnPth([path,'wakeL.dat'])
        if os.path.isfile(fname):
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found')
            Exception('Longitudinal wake file not found')
        SimulData.spos = loadres[:,0]/100    # Rescaling cm to m
        # I know this is correct for ECHO (2015/08/27):
        SimulData.Wll = -loadres[:,1]       # V / pC the minux sign is due to convention

    elif anal_pl in {'dx','dy'}:
        if not silent: print('Loading Transverse Wake file:',end='')
        fname = _jnPth([path,'wakeT.dat'])
        if _os.path.isfile(fname):
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found')
            Exception('Longitudinal wake file not found')
        SimulData.spos = loadres[:,0]/100    # Rescaling cm to m
        setattr(SimulData, 'W'+anal_pl, loadres[:,1])       # V / pC / m  the minux sign is due to convention

    # loading driving bunch info
    loadres = _np.loadtxt(_jnPth([path,'bunch.dat']), skiprows=0)
    sbun = loadres[:,1] / 100     # m
    a = _np.argmin(_np.abs(sbun + sbun[0])) # I want the bunch to be symmetric
    SimulData.sbun = sbun[:a]
    SimulData.bun  = loadres[:a,2] # pC
    SimulData.bunlen = abs(sbun[0] + sbun[1])/ 2 / 5 # ECHO uses 5 sigma
    if not silent:
        print('Bunch length of the driving bunch: {0:7.3g} mm'.format(SimulData.bunlen*1e3))
        print('Data Loaded.\n')


CODES    = {'echoz1': _load_data_ECHOz1,
            'echoz2': _load_data_ECHOz2,
            'gdfidl': _load_data_GdfidL,
            'ace3p' : _load_data_ACE3P,
            'cst'   : _load_data_CST
           }

def load_data(SimulData=None,silent = True):
    if not SimulData: SimulData = EMSimulData()
    path     = SimulData.path
    code     = SimulData.code
    anal_pl  = SimulData.anal_pl

    if not silent: print('*'*60 + '\nLoading Simulation Data')

    #Split the path to try to guess other parameters:
    path_split = set(path.split(_os.path.sep))

    # First try to guess the plane of the analysis if it was not supplied:
    if not anal_pl:
        if not silent: print('Plane of Analysis not supplied, trying to guess from path: ', end='')
        anal_pl_guess = list(ANALYSIS_TYPES & path_split)
        if not anal_pl_guess:
            if not silent: print('could not be guessed.')
            raise Exception('Plane of analysis was not supplied and could not be guessed.')
        else:
            anal_pl = anal_pl_guess[0]
            if not silent: print('ok.')
    if not silent: print('Plane of analysis: {0:s}\n'.format(anal_pl))
    SimulData.anal_pl = anal_pl
    #Now try to guess the code
    if not code:
        if not silent: print('Simulation Code not supplied, trying to guess from path: ', end='')
        code_guess = list(CODES.keys() & path_split)
        if not code_guess:
            if not silent: print('could not be guessed.')
            raise Exception('Simulation Code was not supplied and could not be guessed.')
        else:
            code = code_guess[0]
            if not silent: print('ok.')
    if not silent: print('Code: {0:s}\n'.format(code))
    SimulData.code = code

    CODES[code](SimulData,silent=silent) # changes in SimulData are made implicitly

    print('#'*60+'\n')

def calc_impedance(SimulData, use_win = True, cuttoff = 2, silent = True):

    def _get_impedance(spos,wake,sigt,cutoff):
        dt = (spos[-1]-spos[0]) / (spos.shape[0]-1) / c # frequency scale (Hz):
        VHat = _np.fft.fft(wake) * dt   # fft == \int exp(-i*2pi*f*t/n) G(t) dt
        freq = _np.fft.fftfreq(wake.shape[0],d=dt)
        VHat = _np.fft.fftshift(VHat) # shift the negative frequencies
        freq = _np.fft.fftshift(freq) # to the center of the spectrum
        # Longitudinal position shift to match center of the bunch with zero z:
        w     = 2*_np.pi*freq
        VHat *= _np.exp(-1j*w*spos[0]/c)
        # Deconvolve the Transform with a gaussian bunch:
        Jwlist = _np.exp(-(w*sigt)**2/2)
        Z      = VHat/Jwlist

        wmax  = cutoff/sigt
        indcs = _np.abs(w) <= wmax
        return freq[indcs], Z[indcs]

    # Extracts Needed Variables
    m     = SimulData.anal_pl
    sigt  = SimulData.bunlen / c  # bunch time-length
    spos  = SimulData.spos

    if not silent: print('#'*60 + '\n' + 'Calculating Impedances')

    if use_win:
        if not silent: print('Using Half-Hanning Window')
        # Half Hanning window to zero the end of the signal
        window = _np.hanning(2*spos.shape[0])[spos.shape[0]:]
    else:
        if not silent: print('Not using Window')
        window = _np.ones(spos.shape)

    if not silent: print('Cutoff frequency w = {0:d}/sigmat'.format(cutoff))

    if anal_pl == 'll':
        if not silent: print('Performing FFT of Longitudinal Impedance')
        Wll = SimulData.Wll
        if not Wll or _np.all(Wll == 0):
            if not silent: print('It seems there is no longitudinal wake data.')
            raise Exception('No longitudinal data to perform FFT')
        Wll *= window
        SimulData.freq, Zll = _get_impedance(spos,Wll,sigt,cutoff)
        # I have to take the conjugate of the fft because:
        #fftt == \int exp(-i*2pi*f*t/n) G(t) dt
        #while impedance, according to Chao and Ng, is given by:
        #Z == \int exp(i*2pi*f*t/n) G(t) dt
        SimulData.Zll = Zll.conj()

        if not silent: print('Trying to perform analysis for quadrupolar wakes.')
        for pl in ['qx','qy']:
            print('Horizontal: ' if pl=='qx' else 'Vertical: ', end='')
            Wq = getattr(SimulData,'W'+pl)
            if not Wq or _np.all(Wq == 0):
                if not silent: print('No data found.')
            else:
                if not silent: print('Data Found. ',end='')
                Wq *= window
                _, Zq = _get_impedance(spos,Wq,sigt,cutoff)
                #the Transverse impedance, according to Chao and Ng, is given by:
                #Z == i\int exp(i*2pi*f*t/n) G(t) dt
                setattr(SimulData, 'Z'+pl, 1j*Zq.conj())
                if not silent: print('Impedance Calculated.')

    elif anal_pl in {'dx','dy'}:
        pl = 'Horizontal' if pl=='dx' else 'Vertical'
        if not silent: print('Performing FFT of '+pl+' Impedance.')
        Wd = getattr(SimulData,'W'+anal_pl)
        if not Wd or _np.all(Wd == 0):
            if not silent: print('It seems there is no '+pl+' wake data.')
            raise Exception('No '+pl+' data to perform FFT.')
        else:
            if not silent: print('Data Found. ',end='')
            Wd *= window
            SimulData.freq, Zd = _get_impedance(spos,Wd,sigt,cutoff)
            #the Transverse impedance, according to Chao and Ng, is given by:
            #Z == i\int exp(i*2pi*f*t/n) G(t) dt
            setattr(SimulData, 'Z'+pl, 1j*Zd.conj())
            if not silent: print('Impedance Calculated.')

def calc_loss_factor(SimulData,silent=True):
    # Extracts and Initialize Needed Variables:
    h    = globdata.ringpar.h
    T0   = 2*_np.pi/globdata.ringpar.omega0
    Iavg = globdata.ringpar.Iavg
    sigs = globdata.simpar.bunlen

    wake  = globdata.results.W
    saxis = globdata.results.s
    freq  = globdata.results.freq
    ReZ   = globdata.results.ReZlong

    k   = (freq*2*_np.pi)/c
    ksq = k**2

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
