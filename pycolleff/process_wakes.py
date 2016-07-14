#!/usr/bin/env python3

import os as _os
import re as _re
import sh as _sh
import gzip as _gzip
import pickle as _pickle
import numpy as _np
from scipy import integrate as _scy_int
import matplotlib.pyplot as _plt
from matplotlib import rc as _rc
_rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#_rc('font',**{'family':'serif','serif':['Palatino']})
_rc('text', usetex=True)
from mathphys.constants import light_speed as c
from . import colleff as _colleff
from . import sirius as _sirius

_jnPth = _os.path.sep.join

_si = _sirius.create_ring()

sigvec   = _np.array([2.65, 5.3, 2.65, 4, 10, 10],dtype=float)*1e-3  # bunch length scenarios
Ivec     = _np.array([500, 500, 10, 110, 110, 500],dtype=float)*1e-3 # current scenarios
sigplot  = lambda x:_np.linspace(x,10e-3,num=100)

ANALYSIS_TYPES = {'dx', # horizontal impedance
                  'dy', # vertical impedance
                  'db', # both planes are symmetric
                  'll'  # longitudinal and transverse quadrupolar impedances
                   }

class EMSimulData:
    def __init__(self, path=None, code=None, anal_pl=None):
        self.path      = path or _os.path.abspath('.')  # Path to the wake files
        self.code      = code      # CST, ACE3P, GdfidL, ECHOz1 ECHOz2, ...
        self.anal_pl   = anal_pl   # dx, dy, db, ll

        self.bunlen = 0.0                     # Bunch Length Used in simulation[m]
        self.sbun = _np.array([],dtype=float) # positions where the bunch is defined [m]
        self.bun  = _np.array([],dtype=float) # bunch profile used in the simulation [As/m]
        self.s    = _np.array([],dtype=float) # axis: distance from following to drive bunch [m]
        self.Wll  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/C or V/C/m]
        self.Wdx  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/C or V/C/m]
        self.Wdy  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/C or V/C/m]
        self.Wqx  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/C or V/C/m]
        self.Wqy  = _np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/C or V/C/m]
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
    def klossZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        klossZ,_ = _si.loss_factor(w = self.freq*2*_np.pi, Zl = self.Zll, sigma=sigma)
        return klossZ
    def kckdxZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        kckZ,_ = _si.kick_factor(w = self.freq*2*_np.pi, Z = self.Zdx, sigma=sigma)
        return kckZ
    def kckdyZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        kckZ,_ = _si.kick_factor(w = self.freq*2*_np.pi, Z = self.Zdy, sigma=sigma)
        return kckZ
    def kckqxZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        kckZ,_ = _si.kick_factor(w = self.freq*2*_np.pi, Z = self.Zqx, sigma=sigma)
        return kckZ
    def kckqyZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        kckZ,_ = _si.kick_factor(w = self.freq*2*_np.pi, Z = self.Zqy, sigma=sigma)
        return kckZ

    def PlossW(self, T0=_si.T0, h=_si.harm_num, Iavg=500e-3):
        kW = self.klossW()
        Ploss = kW * Iavg**2 * T0 * 1e12 / h
        return Ploss
    def PlossZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        _,PlossZ,*_ = _si.loss_factor(w = self.freq*2*_np.pi, Zl = self.Zll, sigma=sigma)
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

def _load_data_GdfidL(simul_data,silent=False):

    WAKE_FILENAME_REGEXP = r"[\w-]+W[YXq]{1}_AT_XY.[0-9]{4}"

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
                charge = float(_re.findall(r'[-+]?\d+\.?\d+[eE]?[-+]?\d+',l)[0])
                break
        return charge

    def _get_integration_path(info):
        for line in info:
            if line.find('subtitle=')>=0:
                x,y = (float(val) for val in _re.findall(r'[-+]?\d+\.?\d+[eE]?[-+]?\d+',line))
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
        xd, yd = _get_integration_path(info)
        spos,wake = _np.loadtxt(dados,unpack=True) # dados is a list of strings
        if not silent: print('Charge of the driving bunch: {0:5.3g} pC'.format(charge*1e12))
        if pl == 'll' and (abs(xd) > 1e-10 or abs(yd) > 1e-10) and not silent:
            print('Driving particle not in the origin. Are you sure this is what you want?')
        elif pl !='ll' and abs(xd) < 1e-10 and abs(yd) < 1e-10 and not silent:
            print('The driving bunch is too close to origin. Are you sure this is what you want?')

        a = _np.argmin(_np.diff(spos)) + 1
        sbun   = spos[a:]
        bun    = wake[a:]*charge/_np.trapz(wake[a:],x=sbun) # C
        wake   = -wake[:a]/charge # V/C # minus sign because of convention
        spos   = spos[:a]                # m
        bunlen = -sbun[0]/6            # gdfidl uses a bunch with 6-sigma
        if not silent:
            print('Bunch length of the driving bunch: {0:7.4g} mm'.format(bunlen*1e3))
        return spos, wake, sbun, bun, bunlen, xd, yd

    def _get_transversal_info(path,filelist,pl='qx'):
        stri = 'W{0:s}_AT_XY'.format(pl[1].upper())
        fn = [f for f in f_match if f.find(stri)>=0]
        if not fn:
            if not silent: print('No W{0:s} wake file found. Skipping to next'.format(pl))
            return None
        if not silent: print('{0:2d} W{1:s} wake file found: {2:s}'.format(len(fn),pl,', '.join(fn)))
        dados, info = _load_dados_info(_jnPth([path,fn[0]]))
        charge = _get_charge(info)
        if pl[1] == 'x':
            delta1,_ = _get_integration_path(info)
        else:
            _,delta1 = _get_integration_path(info)
        _, wake1 = _np.loadtxt(dados,unpack=True)
        print('Integration path at {0:s} = {1:8.4g} um '.format(pl[1],delta1*1e6),end='')
        wake = wake1/delta1 / charge # V/C/m
        if len(fn) > 1:
            dados, info = _load_dados_info(_jnPth([path,fn[1]]))
            if pl[1] == 'x':
                delta2,_ = _get_integration_path(info)
            else:
                _,delta2 = _get_integration_path(info)
            _, wake2 = _np.loadtxt(dados,unpack=True)
            print('and {0:8.4g} um'.format(delta2*1e6))
            if pl[0] == 'd':
                wake = (wake1/delta1 - wake2/delta2)/(1/delta1-1/delta2) / charge # V/C
            else:
                wake = (wake1 - wake2)/(delta1-delta2) / charge # V/C/m
        else:
            print()
        return wake

    path     = simul_data.path
    anal_pl  = simul_data.anal_pl

    # list all the files that match the name pattern for wakefields
    f_in_dir = _sh.ls(path).stdout.decode()
    f_match = _re.findall(WAKE_FILENAME_REGEXP,f_in_dir)

    anal_pl_ori = None
    if anal_pl == 'db':
        anal_pl_ori = 'db'
        if not f_match:
            anal_pl = 'dx' if os.path.isdir('dxdpl') else 'dy'
        else:
            anal_pl = 'dx' if [f for f in f_match if f.find('WX_AT_XY')>=0] else 'dy'
        if not silent: print('There is symmetry, calculation performed in the '+anal_pl[1].upper()+' plane.')

    if anal_pl in {'ll'}:
        if not f_match:
            if not silent: print('No files found for longitudinal analysis.')
            raise Exception('No files found for longitudinal analisys')

        #Load longitudinal Wake
        spos, wake, sbun, bun, bunlen, xd, yd = _get_longitudinal_info(path,f_match,pl='ll')
        simul_data.Wll  = wake
        simul_data.s    = spos
        simul_data.bun  = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        # And quadrupolar Wakes, if existent:
        if not silent: print('Loading Horizontal Quadrupolar Wake file:')
        wake = _get_transversal_info(path,f_match,pl='qx') # V/C/m
        if wake is not None: simul_data.Wqx = wake
        if not silent: print('Loading Vertical Quadrupolar Wake file:')
        wake = _get_transversal_info(path,f_match,pl='qy') # V/C/m
        if wake is not None: simul_data.Wqy = wake
        if not silent: print('Longitudinal Data Loaded.')

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
                wk = _get_transversal_info(ext_path,f_match,pl=anal_pl) # V/C
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
            simul_data.s      = spos[0]
            simul_data.bun    = bun[0]
            simul_data.sbun   = sbun[0]
            simul_data.bunlen = bunlen[0]
            setattr(simul_data,'W'+anal_pl, (wake[0]-wake[1])/(delta[0]-delta[1])) # V/C/m
        else:
            spos, wake, sbun, bun, bunlen, xd, yd = _get_longitudinal_info(path,f_match,pl=anal_pl)
            simul_data.s      = spos
            simul_data.bun    = bun
            simul_data.sbun   = sbun
            simul_data.bunlen = bunlen

            if not silent:
                print('Loading {0:s} Dipolar Wake file:'.format(
                      'Horizontal' if anal_pl=='dx' else 'Vertical'))
            wake = _get_transversal_info(path,f_match,pl=anal_pl) # V/C
            if wake is not None:
                delta = xd if anal_pl=='dx' else yd
                setattr(simul_data,'W'+anal_pl,wake/delta) # V/C/m
            else:
                print('Actually there is something wrong, these wake files should be here.')
                raise Exception('Transverse {0:s} dipolar wake files not found'.format(
                                'Horizontal' if anal_pl=='dx' else 'Vertical'))
        if not silent: print('Transverse Data Loaded.')
    else:
        if not silent: print('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))
        raise Exception('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))

    if anal_pl_ori:
        anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
        if not silent: print('There is symmetry. Copying the data from the '+
                    '{0:s} plane to the {1:s} plane'.format(anal_pl[1].upper(),anal_pl_comp[1].upper()))
        setattr(simul_data, 'W'+anal_pl_comp, getattr(simul_data,'W'+anal_pl).copy())

def _load_data_ECHOz1(simul_data,silent=False):

    path     = simul_data.path
    anal_pl  = simul_data.anal_pl

    if anal_pl=='ll':
        if not silent: print('Loading longitudinal Wake file:',end='')
        fname = _jnPth([path,'wake.dat'])
        if _os.path.isfile(fname):
            if not silent: print('Data found.')
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found.')
            Exception('Longitudinal wake file not found.')
    else:
        if not silent: print('ECHOz1 only calculates longitudinal wake.')
        raise Exception('ECHOz1 only calculates longitudinal wake.')

    simul_data.s = loadres[:,0]/100    # Rescaling cm to m
    simul_data.Wll = -loadres[:,1] * 1e12 # V/C/m   the minus sign is due to convention

    # loading driving bunch info
    fname = _jnPth([path,'bunch.dat'])
    if _os.path.isfile(fname):
        if not silent: print('Loading bunch length from file bunch.dat')
        loadres = _np.loadtxt(fname, skiprows=0)
        sbun = loadres[:,0] / 100     # m
        a = _np.argmin(_np.abs(sbun + sbun[0])) + 1 # I want the bunch to be symmetric
        simul_data.sbun = sbun[:a]
        simul_data.bun  = loadres[:a,1] # C
        simul_data.bunlen = abs(sbun[0]) / 5 # ECHO uses 5 sigma
    else:
        if not silent: print('File bunch.dat not found. Loading bunch length from wake.dat')
        spos = simul_data.s.copy()
        ds   = spos[1]-spos[0]
        sbun = _np.hstack([spos[0]-ds,spos])
        bunlen = abs(sbun[0]) / 5
        a    = _np.argmin(_np.abs(sbun + sbun[0])) + 1
        sbun = sbun[:a]
        simul_data.bunlen = bunlen
        simul_data.sbun   = sbun
        simul_data.bun    = _np.exp(-sbun**2/(2*bunlen**2))/(_np.sqrt(2*_np.pi)*bunlen)
    if not silent:
        print('Bunch length of the driving bunch: {0:7.3g} mm'.format(simul_data.bunlen*1e3))
        print('Data Loaded.')

def _load_data_ECHOz2(simul_data,silent=False):

    path     = simul_data.path
    anal_pl  = simul_data.anal_pl

    anal_pl_ori = None
    if anal_pl == 'db':
        anal_pl_ori = 'db'
        anal_pl     = 'dy'
        if not silent: print('Even though there is symmetry, I am loading data to the Y plane.')

    if anal_pl=='ll':
        if not silent: print('Loading longitudinal Wake file:',end='')
        fname = _jnPth([path,'wakeL.dat'])
        if os.path.isfile(fname):
            if not silent: print('Data found.')
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found.')
            Exception('Longitudinal wake file not found.')
        simul_data.s = loadres[:,0]/100    # Rescaling cm to m
        simul_data.Wll = -loadres[:,1] * 1e12 # V/C the minus sign is due to convention

    elif anal_pl in {'dx','dy'}:
        if not silent: print('Loading Transverse Wake file:',end='')
        fname = _jnPth([path,'wakeT.dat'])
        if _os.path.isfile(fname):
            if not silent: print('Data found.')
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found.')
            Exception('Transverse wake file not found.')
        simul_data.s = loadres[:,0]/100    # Rescaling cm to m
        setattr(simul_data, 'W'+anal_pl, loadres[:,1] * 1e12) # V/C/m  the minus sign is due to convention
    else:
        if not silent: print('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))
        raise Exception('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))

    # loading driving bunch info
    fname = _jnPth([path,'bunch.dat'])
    if _os.path.isfile(fname):
        if not silent: print('Loading bunch length from file bunch.dat')
        loadres = _np.loadtxt(fname, skiprows=0)
        sbun = loadres[:,1] / 100     # m
        a = _np.argmin(_np.abs(sbun + sbun[0])) + 1 # I want the bunch to be symmetric
        simul_data.sbun = sbun[:a]
        simul_data.bun  = loadres[:a,2] # C
        simul_data.bunlen = abs(sbun[0]) / 5 # ECHO uses 5 sigma
    else:
        if not silent: print('File bunch.dat not found. Loading bunch length from wakeT.dat')
        spos = simul_data.s.copy()
        ds   = spos[1]-spos[0]
        spos += ds/2
        sbun = _np.hstack([spos[0]-ds,spos])
        bunlen = abs(sbun[0]) / 5
        a    = _np.argmin(_np.abs(sbun + sbun[0])) + 1
        sbun = sbun[:a]
        simul_data.bunlen = bunlen
        simul_data.sbun   = sbun
        simul_data.bun    = _np.exp(-sbun**2/(2*bunlen**2))/(_np.sqrt(2*_np.pi)*bunlen)
    if not silent:
        print('Bunch length of the driving bunch: {0:7.3g} mm'.format(simul_data.bunlen*1e3))
        print('Data Loaded.')

    if anal_pl_ori:
        anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
        if not silent: print('There is symmetry. Copying the data from the '+
                    '{0:s} plane to the {1:s} plane'.format(anal_pl[1].upper(),anal_pl_comp[1].upper()))
        setattr(simul_data, 'W'+anal_pl_comp, getattr(simul_data,'W'+anal_pl).copy())

def _load_data_ECHOzR(simul_data,silent=False):

    WAKE_FILENAME_REGEXP = r"(wakeL_([0-9]{2}).txt)" # the older .dat files are not treated

    def _load_dados_info(fname,mode, bound_cond):
        with open(fname) as f:
            f.readline()
            a = f.readline()
        mstep, offset, wid, bunlen = _np.fromstring(a[1:],sep='\t')
        offset = int(offset)
        dados = _np.loadtxt(fname,skiprows=3)
        spos  = dados[:,0]*1e-2  # cm to m

        y0 = y = mstep*offset * 1e-2 # cm to m
        Kxm = _np.pi/wid*mode * 1e2  # 1/cm to 1/m
        if bound_cond.startswith('elec'): # minus sign is due to convention
            Wm  = -dados[:,1+offset]/_np.sinh(Kxm * y0)/_np.sinh(Kxm * y) * 1e-2 * 1e12 # cm*V/pC to m*V/C
        else: # minus sign is due to convention
            Wm  = -dados[:,1+offset]/_np.cosh(Kxm * y0)/_np.cosh(Kxm * y) * 1e-2 * 1e12 # cm*V/pC to m*V/C
        return spos, Wm, mstep*1e-2, wid*1e-2, bunlen*1e-2

    path     = simul_data.path
    anal_pl  = simul_data.anal_pl

    anal_pl_ori = None
    if anal_pl == 'db':
        if not silent: print('Problem: All geometries simulated with ECHOzR does not have symmetry.')
        raise Exception('Problem: All geometries simulated with ECHOzR does not have symmetry.')

    if anal_pl == 'll':
        if not silent: print('Looking for data files in subfolder magn.')
        bound_cond = 'magnetic'
    elif anal_pl in {'dx','dy'}:
        if not silent: print('Looking for data files in subfolder elec.')
        bound_cond = 'electric'

    pname = _jnPth([path,bound_cond[:4]])
    if not _os.path.isdir(pname):
        pname = path
        if not silent:
            print('Subfolder not found. It would be better to create the subfolder and put the files there...')
            print('Looking for files in the current folder:')

    f_in_dir = _sh.ls(pname).stdout.decode()
    f_match = sorted(_re.findall(WAKE_FILENAME_REGEXP,f_in_dir))
    if not f_match:
        if not silent: print('Files not found. ',end='')
        raise Exception('Files not found.')


    if not silent:
        print('Files found. I am assuming the simulation was performed '+
              'with {0:s} boundary condition.'.format(bound_cond))
        print('Modes found: ' + ', '.join([m for _,m in f_match]))
        print('Loading data from files')

    spos, W, mode, mesh_size, width, bunlen = [], [], [], [], [], []
    for fn, m in f_match:
        s, Wm, ms, wid, bl = _load_dados_info(_jnPth([pname,fn]),int(m),bound_cond)
        mode.append(int(m))
        spos.append(s)
        W.append(Wm)
        mesh_size.append(ms)
        width.append(wid)
        bunlen.append(bl)

    cond = False
    for i in range(1,len(mode)):
        cond |= len(spos[i]) != len(spos[0])
        cond |= not _np.isclose(mesh_size[i], mesh_size[0], rtol=1e-5,atol=0)
        cond |= not _np.isclose(width[i], width[0], rtol=0,atol=1e-7)
        cond |= not _np.isclose(bunlen[i], bunlen[0], rtol=1e-5,atol=0)
        if cond:
            message = 'Parameters of file {0:s} differ from {1:s}.'.format(f_match[i][0],f_match[0][0])
            if not silent: print(message)
            raise Exception(message)


    simul_data.s = spos[0]
    simul_data.bunlen = bunlen[0]
    a = _np.argmin(_np.abs(spos[0] + spos[0][0])) + 1 # I want the bunch to be symmetric
    sbun = spos[0][:a]
    simul_data.sbun = sbun
    simul_data.bun  =_np.exp(-sbun**2/(2*bunlen[0]**2))/(_np.sqrt(2*_np.pi)*bunlen[0])
    if not silent:
        print('Bunch length of the driving bunch: {0:7.4g} mm'.format(simul_data.bunlen*1e3))
        print('Width of the simulated geometry:   {0:7.4g} mm'.format(width[0]*1e3))
        print('Mesh step used in the simulation:  {0:7.4g} um'.format(mesh_size[0]*1e6))
        print('All Data Loaded.')

    if anal_pl=='ll':
        if not silent: print('Calculating longitudinal Wake from data:')
        Wll, frac =  None, 1
        for i in range(len(mode)):
            if mode[i] == 1:
                Wll = W[i].copy()
            elif mode[i] % 2:
                Wll += W[i] #only odd terms
                frac = _np.max(_np.abs(W[i]/Wll))
        if Wll is None:
            if not silent: print('There is none odd mode to calculate Longitudinal Wake.')
        else:
            if not silent: print('Maximum influence of last mode in the final result is: {0:5.2f}%'.format(frac*100))
            Wll *= 1/width[0]/2
            simul_data.Wll = Wll

        if not silent: print('Calculating Quadrupolar Wake from data:')
        Wq, frac =  None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wq = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wq += W[i] * Kxm**2 #only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wq))
        if Wq is None:
            if not silent: print('There is none odd mode to calculate Quadrupolar Wake.')
        else:
            if not silent: print('Maximum influence of last mode in the final result is: {0:5.2f}%'.format(frac*100))
            Wq *= 1/width[0]/2
            Wq  = _scy_int.cumtrapz(Wq,x=spos[0],initial=0)
            simul_data.Wqy =  Wq
            simul_data.Wqx = -Wq

        if not silent: print('Calculating Dipolar Horizontal Wake from data:')
        Wdx, frac =  None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 2:
                Wdx = W[i].copy() * Kxm**2
            elif not mode[i] % 2:
                Wdx += W[i] * Kxm**2 #only even terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wdx))
        if Wdx is None:
            if not silent: print('There is none even mode to calculate Dipolar Horizontal Wake.')
        else:
            if not silent: print('Maximum influence of last mode in the final result is: {0:5.2f}%'.format(frac*100))
            Wdx *= 1/width[0]/2
            Wdx  = _scy_int.cumtrapz(Wdx,x=spos[0],initial=0)
            simul_data.Wdx = -Wdx

    elif anal_pl in {'dx','dy'}:
        if not silent: print('Calculating Dipolar Vertical Wake from data:')
        Wdy, frac =  None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wdy = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wdy += W[i] * Kxm**2 #only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wdy))
        if Wdy is None:
            if not silent: print('There is none even mode to calculate Dipolar Vertical Wake.')
        else:
            if not silent: print('Maximum influence of last mode in the final result is: {0:5.2f}%'.format(frac*100))
            Wdy *= 1/width[0]/2
            Wdy  = _scy_int.cumtrapz(Wdy,x=spos[0],initial=0)
            simul_data.Wdy = -Wdy
    else:
        if not silent: print('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))
        raise Exception('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))


CODES    = {'echoz1': _load_data_ECHOz1,
            'echoz2': _load_data_ECHOz2,
            'echozr': _load_data_ECHOzR,
            'gdfidl': _load_data_GdfidL,
            'ace3p' : _load_data_ACE3P,
            'cst'   : _load_data_CST
           }

def load_raw_data(simul_data=None,silent = False):
    if not simul_data: simul_data = EMSimulData()
    path     = simul_data.path
    code     = simul_data.code
    anal_pl  = simul_data.anal_pl

    if not silent: print('#'*60 + '\nLoading Simulation Data')

    #Split the path to try to guess other parameters:
    path_split = set(path.lower().split(_os.path.sep))

    # First try to guess the plane of the analysis if it was not supplied:
    if not anal_pl:
        if not silent: print('Plane of Analysis not supplied, trying to guess from path: ', end='')
        anal_pl_guess = list(ANALYSIS_TYPES & path_split)
        if not anal_pl_guess:
            if not silent: print('could not be guessed.')
            raise Exception('Plane of analysis was not supplied and could not be guessed.')
        else:
            anal_pl = anal_pl_guess[0]
    if not silent: print(anal_pl)
    simul_data.anal_pl = anal_pl


    #Now try to guess the code
    if not code:
        if not silent: print('Simulation Code not supplied, trying to guess from path: ', end='')
        code_guess = list(CODES.keys() & path_split)
        if not code_guess:
            if not silent: print('could not be guessed.')
            raise Exception('Simulation Code was not supplied and could not be guessed.')
        else:
            code = code_guess[0]
    if not silent: print(code)
    simul_data.code = code

    CODES[code](simul_data,silent=silent) # changes in simul_data are made implicitly

    print('#'*60+'\n')

def calc_impedance(simul_data, use_win = True, cutoff = 2, silent = False):

    def _get_impedance(spos,wake,sigt,cutoff):
        dt = (spos[-1]-spos[0]) / (spos.shape[0]-1) / c # frequency scale (Hz):
        VHat = _np.fft.fft(wake) * dt   # fft == \int exp(-i*2pi*f*t/n) G(t) dt
        freq = _np.fft.fftfreq(wake.shape[0],d=dt)
        VHat = _np.fft.fftshift(VHat) # shift the negative frequencies
        freq = _np.fft.fftshift(freq) # to the center of the spectrum
        # Longitudinal position shift to match center of the bunch with zero z:
        w     = 2*_np.pi*freq
        VHat *= _np.exp(-1j*w*spos[0]/c)
        # find the maximum useable frequency
        wmax  = cutoff/sigt
        indcs = _np.abs(w) <= wmax
        # Deconvolve the Transform with a gaussian bunch:
        Jwlist = _np.exp(-(w*sigt)**2/2)
        Z      = VHat[indcs]/Jwlist[indcs]
        return freq[indcs], Z

    # Extracts Needed Variables
    sigt    = simul_data.bunlen / c  # bunch time-length
    spos    = simul_data.s

    if not silent: print('#'*60 + '\n' + 'Calculating Impedances')
    if use_win:
        if not silent: print('Using Half-Hanning Window')
        # Half Hanning window to zero the end of the signal
        window = _np.hanning(2*spos.shape[0])[spos.shape[0]:]
    else:
        if not silent: print('Not using Window')
        window = _np.ones(spos.shape[0])

    if not silent: print('Cutoff frequency w = {0:d}/sigmat'.format(cutoff))

    for pl in ['ll','dx','dy','qx','qy']:
        if not silent: print('Performing FFT on W{0:s}: '.format(pl),end='')
        Wpl = getattr(simul_data,'W'+pl)
        if Wpl is None or _np.all(Wpl == 0):
            if not silent: print('No Data found.')
            continue
        if not silent: print('Data found. ',end='')
        Wpl *= window
        simul_data.freq, Zpl = _get_impedance(spos,Wpl,sigt,cutoff)
        if pl =='ll':
            # I have to take the conjugate of the fft because:
            #fftt == \int exp(-i*2pi*f*t/n) G(t) dt
            #while impedance, according to Chao and Ng, is given by:
            #Z == \int exp(i*2pi*f*t/n) G(t) dt
            simul_data.Zll = Zpl.conj()
        else:
            #the Transverse impedance, according to Chao and Ng, is given by:
            #Z == i\int exp(i*2pi*f*t/n) G(t) dt
            setattr(simul_data, 'Z'+pl, 1j*Zpl.conj())
        if not silent: print('Impedance Calculated.')

    if not silent: print('#'*60 + '\n')

def predefined_studies(simul_data,silent=False,save_figs=False):
    raise NotImplementedError('Not implemented Yet')
    # First Plot the short range Wake

def plot_short_range_wake(simul_data,silent=False,save_figs=False,pth2sv=None,show=False):

    pls    = ('Wll','Wdx','Wdy','Wqx','Wqy')
    titles = ('Longitudinal',
              'Dipolar Horizontal',
              'Dipolar Vertical',
              'Quadrupolar Horizontal',
              'Quadrupolar Vertical',
             )
    ylabels= (r'$W_l$ [V/pC]',
              r'$W_{{D_x}}$ [V/pC/m]',
              r'$W_{{D_y}}$ [V/pC/m]',
              r'$W_{{Q_x}}$ [V/pC/m]',
              r'$W_{{Q_y}}$ [V/pC/m]')

    sbun = simul_data.sbun
    sigs = simul_data.bunlen
    spos = simul_data.s
    for pl,tit,yl in zip(pls,titles,ylabels):
        wake = getattr(simul_data,pl)*1e-12 # V/C -> V/pC
        if wake is None or _np.all(wake==0): continue
        bunchshape = simul_data.bun * (wake.max()/simul_data.bun.max())

        _plt.figure()
        _plt.plot(sbun*1000,bunchshape,'b',linewidth=2,label='Bunch Shape')
        _plt.plot(spos*1000,wake,'r',linewidth=2,label='Wake Potential')
        _plt.grid(True)
        _plt.xlabel('s [mm]',fontsize=13)
        fname = pl+'_ShortRange'
        _plt.title (tit,fontsize=13)
        _plt.ylabel(yl,fontsize=13)
        _plt.xlim([spos[0]*1000, 7000*sigs])
        _plt.ylim([wake.min()*1.1, wake.max()*1.1])
        _plt.legend(loc='best')
        if save_figs: _plt.savefig(_jnPth((pth2sv,fname+'.svg')))
    if show: _plt.show()

def plot_long_range_wake(simul_data,silent=False,save_figs=False,pth2sv=None,show=False):

    pls    = ('Wll','Wdx','Wdy','Wqx','Wqy')
    titles = ('Longitudinal',
              'Dipolar Horizontal',
              'Dipolar Vertical',
              'Quadrupolar Horizontal',
              'Quadrupolar Vertical',
             )
    ylabels= (r'$W_l$ [V/pC]',
              r'$W_{{D_x}}$ [V/pC/m]',
              r'$W_{{D_y}}$ [V/pC/m]',
              r'$W_{{Q_x}}$ [V/pC/m]',
              r'$W_{{Q_y}}$ [V/pC/m]')

    spos = simul_data.s
    for pl,tit,yl in zip(pls,titles,ylabels):
        wake = getattr(simul_data,pl)*1e-12 # V/C -> V/pC
        if wake is None or _np.all(wake==0): continue

        _plt.figure()
        _plt.plot(spos*1000,wake,'r',linewidth=2)
        _plt.grid(True)
        _plt.xlabel('s [mm]',fontsize=13)
        _plt.title (tit,fontsize=13)
        _plt.ylabel(yl,fontsize=13)
        _plt.ylim([wake.min()*1.1, wake.max()*1.1])
        fname = pl+'_LongRange'
        if save_figs: _plt.savefig(_jnPth((pth2sv,fname+'.svg')))
    if show: _plt.show()

def plot_impedances(simul_data,silent=False,save_figs=False,pth2sv=None,show=False):

    pls    = ('Zll','Zdx','Zdy','Zqx','Zqy')
    titles = ('Longitudinal',
              'Dipolar Horizontal',
              'Dipolar Vertical',
              'Quadrupolar Horizontal',
              'Quadrupolar Vertical',
             )
    ylabels= (r'$Z_l$ [$\Omega$]',
              r'$Z_{{D_x}}$ [$\Omega$/m]',
              r'$Z_{{D_y}}$ [$\Omega$/m]',
              r'$Z_{{Q_x}}$ [$\Omega$/m]',
              r'$Z_{{Q_y}}$ [$\Omega$/m]')

    freq = simul_data.freq
    for pl,tit,yl in zip(pls,titles,ylabels):
        Z = getattr(simul_data,pl)
        if Z is None or _np.all(Z==0): continue

        _plt.figure()
        _plt.plot(freq/1e9,Z.real,'r',linewidth=2,label='Re')
        _plt.plot(freq/1e9,Z.imag,'b--',linewidth=2,label='Im')
        _plt.xlabel('Frequency [GHz]',fontsize=13)
        _plt.grid(True)
        _plt.title (tit,fontsize=13)
        _plt.ylabel(yl,fontsize=13)
        _plt.legend (loc='best')
        _plt.xlim(_np.array(freq[[0,-1]],dtype=float)/1e9)
        if save_figs: _plt.savefig(_jnPth((pth2sv,pl+'.svg')))
    if show: _plt.show()

def calc_plot_losskick_factors(simul_data,silent=False,save_figs=False,pth2sv=None,show=False):
    # Extracts and Initialize Needed Variables:
    h    = globdata.ringpar.h
    T0   = 2*_np.pi/globdata.ringpar.omega0
    Iavg = globdata.ringpar.Iavg
    sigs = globdata.simpar.bunlen

    wake  = globdata.results.W
    saxis = globdata.results.s
    freq  = globdata.results.freq
    ReZ   = globdata.results.ReZlong
    # Extracts and Initialize Needed Variables:
    sigs  = globdata.simpar.bunlen
    wake  = globdata.results.W
    saxis = globdata.results.s
    freq  = globdata.results.freq
    ImZ   = globdata.results.ImZt

    sigmasq = sigs**2
    w   = (freq*2*_np.pi)
    k   = w/c
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

    # Calculates klossW
    ss    = saxis**2
    rhos  = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-ss/(2*sigs**2))
    kW    = _np.trapz(wake*rhos, x=saxis)
    Ploss = kW * Iavg**2 * T0 * 1e12 / h

    # Print loss factor calculated in both ways

    print('klossZ = {0:6.5g} mV/pC'.format(kZ*1000))
    print('klossW = {0:6.5g} mV/pC'.format(kW*1000))
    print('Ploss  = {0:6.5g} W     (for {1:5.4g} mA avg current)'.format(Ploss,Iavg*1000))

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

    # Print kick factor calculated in both ways:
    print('Kick_Z = {0:6.5g} V/pC/m'.format(kickZ))
    print('Kick_W = {0:6.5g} V/pC/m'.format(kickW))


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

show_now = _plt.show

def save_processed_data(simul_data,silent=False,pth2sv=None):

    if not silent: print('#'*60 + '\nSaving Processed data:')
    path = simul_data.path
    spos = simul_data.s
    freq = simul_data.freq

    if pth2sv is None:
        if not silent: print('Saving in the same folder of the raw data')
        pth2sv = path
    elif type(pth2sv) is str:
        if not silent: print('Saving to subfolder: ' + pth2sv)
        pth2sv = _jnPth([path,pth2sv])
        if not _os.path.isdir(pth2sv):
            if not silent: print('Folder does not exist. Creating it...')
            _os.mkdir(pth2sv)
    else:
        if not silent: print('pth2sv must be a string or None object')
        raise Exception('pth2sv must be a string or None')

    #Save wakes
    for par in ['Wll','Wdx','Wdy','Wqx','Wqy']:
        unit = 'V/C' if par == 'Wll' else 'V/C/m'
        header = '{0:30s} {1:30s}'.format('s [m]', '{0:s} [{1:s}]'.format(par,unit))
        fname = _jnPth([pth2sv,par+'.gz'])
        wake  = getattr(simul_data,par)
        if wake is None or _np.all(wake == 0): continue
        if not silent: print('Saving '+ par + ' data to .gz file')
        _np.savetxt(fname,_np.array([spos,wake]).transpose(),
                    fmt=['%30.16g','%30.16g'], header=header)

    #Save Impedances
    for par in ['Zll','Zdx','Zdy','Zqx','Zqy']:
        unit = 'Ohm' if par == 'Zll' else 'Ohm/m'
        header = '{0:30s} {1:30s} {2:30s}'.format('Frequency [GHz]',
                                            'Re{0:s} [{1:s}]'.format(par,unit),
                                            'Im{0:s} [{1:s}]'.format(par,unit))
        fname = _jnPth([pth2sv,par+'.gz'])
        Z  = getattr(simul_data,par)
        if Z is None or _np.all(Z == 0): continue
        if not silent: print('Saving '+ par + ' data to .gz file')
        _np.savetxt(fname,_np.array([freq*1e9,Z.real,Z.imag]).transpose(),
                    fmt=['%30.16g','%30.16g','%30.16g'], header=header)

    if not silent: print('Saving the Complete EMSimulData structure to a .pickle file.')
    with _gzip.open(_jnPth((pth2sv,'SimulData.pickle')), 'wb') as f:
        _pickle.dump(simul_data,f,_pickle.HIGHEST_PROTOCOL)

    if not silent: print('All Data Saved\n' + '#'*60)

def load_processed_data(filename):
    with _gzip.open(filename,'rb') as fh:
        simul_data = _pickle.load(fh)
    return simul_data

def analysis_example():

    analysis = '''
    #!/usr/bin/env python3

    import optparse
    import os
    import pycolleff.process_wakes as proc_wake

    def main(pth2sv='analysis',silent = False):
        simul_data = proc_wake.EMSimulData()
        proc_wake.load_raw_data(simul_data,silent=False)
        proc_wake.calc_impedance(simul_data,silent=False)
        proc_wake.save_processed_data(simul_data,silent=False,pth2sv=pth2sv)
        return simul_data

    if __name__ == '__main__':

        # configuration of the parser for the arguments
        parser = optparse.OptionParser()
        parser.add_option('-p','--noplot',dest='plot',action='store_true',
                          help="Show results", default=False)
        parser.add_option('-s','--silent',dest='silent',action='store_true',
                  help="Print progress", default=False)
        parser.add_option('-c','--calc',dest='calc',action='store_true',
                          help="Calculate results", default=False)
        parser.add_option('--pth2sv',dest='pth2sv',type='str',
                          help="Path to save the data. Relative to the current folder",
                          default = 'analysis')
        (opts, _) = parser.parse_args()

        plot = not opts.plot
        silent = opts.silent
        pth2sv = opts.pth2sv
        file_name = os.path.sep.join([os.path.abspath('.'),
                                      pth2sv,'SimulData.pickle'])

        if opts.calc:
            simul_data = main(pth2sv = pth2sv, silent=silent)
            salva = True
        else:
            simul_data = proc_wake.load_processed_data(file_name)
            salva = False

        proc_wake.plot_short_range_wake(simul_data,silent=silent,
                                    save_figs=salva,pth2sv=pth2sv,show=False)
        proc_wake.plot_long_range_wake(simul_data,silent=silent,
                                    save_figs=salva,pth2sv=pth2sv,show=False)
        proc_wake.plot_impedances(simul_data,silent=silent,
                                save_figs=salva,pth2sv=pth2sv,show=False)
        proc_wake.calc_plot_losskick_factors(simul_data,silent=silent,
                                    save_figs=salva,pth2sv=pth2sv,show=False)
        proc.show_now()
    '''
    print(analysis)
    return None
