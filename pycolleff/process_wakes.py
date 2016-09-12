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
from pyaccel import naff as _naff
from . import colleff as _colleff
from . import sirius as _sirius


## Trought the code I am assuming:
# s positive means particle behind source -->  Wl, Wt = 0 s < 0
# Wl(s) = -c/Q * int El(ct-s,t) dt
# Wx(s) = - int_-inf^s dWl/dx ds'
# Zl =   int exp(i*w*s) Wl(s) ds
# Zx = i*int exp(i*w*s) Wx(s) ds

_jnPth = _os.path.sep.join
_si = _sirius.create_ring()
DEFAULT_FNAME_SAVE = 'SimulData.pickle'
FNAME_ECHOZ1   = r"wake.dat"
FNAME_ECHOZ2   = r"wake[LT]{1}.dat"
FNAME_ECHOZR2D = r"(wakeL_([0-9]{2}).txt)" # the older .dat files are not treated
FNAME_GDFIDL   = r"[\w-]+W[YXq]{1}_AT_XY.[0-9]{4}"

ANALYSIS_TYPES = {'dx', # horizontal impedance
                  'dy', # vertical impedance
                  'db', # both planes are symmetric
                  'll'  # longitudinal and transverse quadrupolar impedances
                  }

PLANES  = ('ll','dx','dy','qx','qy')
TITLES  = {'ll':'Longitudinal',
           'dx':'Dipolar Horizontal',
           'dy':'Dipolar Vertical',
           'qx':'Quadrupolar Horizontal',
           'qy':'Quadrupolar Vertical',
          }
WAKE_YLABELS= {'ll':r'$W_l$ [V/pC]',
               'dx':r'$W_{{D_x}}$ [V/pC/m]',
               'dy':r'$W_{{D_y}}$ [V/pC/m]',
               'qx':r'$W_{{Q_x}}$ [V/pC/m]',
               'qy':r'$W_{{Q_y}}$ [V/pC/m]'
              }
IMPS_YLABELS= {'ll':r'$Z_l$ [$\Omega$]',
               'dx':r'$Z_{{D_x}}$ [$\Omega$/m]',
               'dy':r'$Z_{{D_y}}$ [$\Omega$/m]',
               'qx':r'$Z_{{Q_x}}$ [$\Omega$/m]',
               'qy':r'$Z_{{Q_y}}$ [$\Omega$/m]'
              }

class EMSimulData:
    def __init__(self, code=None):
        self.code      = code      # CST, ACE3P, GdfidL, ECHOz1 ECHOz2, ...

        self.bunlen = 0.0                     # Bunch Length Used in simulation[m]
        self.sbun = _np.array([],dtype=float) # positions where the bunch is defined [m]
        self.bun  = _np.array([],dtype=float) # bunch profile used in the simulation [As/m]
        self.s    = _np.array([],dtype=float) # axis: distance from following to drive bunch [m]
        self.Wll  = _np.array([],dtype=float) # Longitudinal Wakepotential [V/C]
        self.Wdx  = _np.array([],dtype=float) # Dipolar Horizontal Wakepotential [V/C/m]
        self.Wdy  = _np.array([],dtype=float) # Dipolar Vertical Wakepotential [V/C/m]
        self.Wqx  = _np.array([],dtype=float) # Quadrupolar Horizontal Wakepotential [V/C/m]
        self.Wqy  = _np.array([],dtype=float) # Quadrupolar Vertical Wakepotential [V/C/m]
        self.freq = _np.array([],dtype=float) # axis: frequency obtained from FFT [GHz]
        self.Zll  = _np.array([],dtype=complex) # Longitudinal Impedance [Ohm]
        self.Zdx  = _np.array([],dtype=complex) # Dipolar Horizontal Impedance [Ohm]
        self.Zdy  = _np.array([],dtype=complex) # Dipolar Vertical Impedance [Ohm]
        self.Zqx  = _np.array([],dtype=complex) # Quadrupolar Horizontal Impedance [Ohm]
        self.Zqy  = _np.array([],dtype=complex) # Quadrupolar Vertical Impedance [Ohm]
        self._klossW  = None
        self._kckdxW  = None
        self._kckdyW  = None
        self._kckqxW  = None
        self._kckqyW  = None

    def klossW(self):
        if self._klossW: return self._klossW
        T0, sigs, spos = _si.T0, self.bunlen, self.s
        wake = self.Wll
        if wake is None or _np.all(wake==0): return None
        rhos  = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW    = _np.trapz(wake*rhos, x=spos)
        self._klossW = kW
        return kW
    def PlossW(self, T0=_si.T0, h=_si.harm_num, Iavg=500e-3):
        kW = self.klossW()
        Ploss = kW * Iavg**2 * T0 * 1e12 / h
        return Ploss
    def kick_factorW(self,pl='dy'):
        kick = getattr(self,'_kck'+pl+'W')
        if kick: return kick
        T0, sigs, spos = _si.T0, self.bunlen, self.s
        wake = getattr(self,'W'+pl)
        if wake is None or _np.all(wake==0): return None
        rhos  = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW    = _np.trapz(wake*rhos, x=spos)
        setattr(self,'_kck'+pl+'W', kW)
        return kW

    def klossZ(self,sigma=2.65e-3,n=1):
        _si.nbun = n
        klossZ,*_ = _si.loss_factor(w = self.freq*2*_np.pi, Zl = self.Zll, sigma=sigma)
        return klossZ
    def kick_factorZ(self,pl='dy',sigma=2.65e-3,n=1):
        _si.nbun = n
        Z = getattr(self,'Z'+pl)
        if Z is None or _np.all(Z==0): return None
        kckZ,*_ = _si.kick_factor(w = self.freq*2*_np.pi, Z = Z, sigma=sigma)
        return kckZ
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

def _load_data_GdfidL(simul_data,path,anal_pl,silent=False):

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

    # list all the files that match the name pattern for wakefields
    f_in_dir = _sh.ls(path).stdout.decode()
    f_match = _re.findall(FNAME_GDFIDL,f_in_dir)

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
                ext_path = _jnPth([path,anal_pl+sub_fol])
                if not silent: print('Looking for '+anal_pl+sub_fol+' subfolder:')
                if not _os.path.isdir(ext_path):
                    if not silent: print('For non-symmetric structures, there must '
                                   'be subfolders {0:s}dpl {0:s}dmi with the data'.format(anal_pl))
                    raise Exception('Files not found')
                # list all the files that match the pattern
                f_in_dir = _sh.ls(ext_path).stdout.decode()
                f_match = _re.findall(FNAME_GDFIDL,f_in_dir)
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
                if wk is not None:
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

def _load_data_ECHOz1(simul_data,path,anal_pl,silent=False):

    if anal_pl=='ll':
        if not silent: print('Loading longitudinal Wake file:',end='')
        fname = _jnPth([path,FNAME_ECHOZ1])
        if _os.path.isfile(fname):
            if not silent: print('Data found.')
            loadres = _np.loadtxt(fname, skiprows=0)
        else:
            if not silent: print('Not found.')
            raise Exception('Longitudinal wake file not found.')
    else:
        if not silent: print('ECHOz1 only calculates longitudinal wake.')
        raise Exception('ECHOz1 only calculates longitudinal wake.')

    simul_data.s = loadres[:,0]/100    # Rescaling cm to m
    simul_data.Wll = -loadres[:,1] * 1e12 # V/C/m   the minus sign is due to convention

    # loading driving bunch info
    if not silent: print('Loading bunch length from wake.dat')
    sbun = simul_data.s.copy()
    ds   = sbun[1]-sbun[0]
    bunlen = abs(sbun[0]-ds/2) / 5
    a    = _np.argmin(_np.abs(sbun + sbun[0])) + 1
    sbun = sbun[:a]
    simul_data.bunlen = bunlen
    simul_data.sbun   = sbun
    simul_data.bun    = _np.exp(-sbun**2/(2*bunlen**2))/(_np.sqrt(2*_np.pi)*bunlen)
    if not silent:
        print('Bunch length of the driving bunch: {0:7.3g} mm'.format(bunlen*1e3))
        print('Data Loaded.')

def _load_data_ECHOz2(simul_data,path,anal_pl,silent=False):

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
            spos,wl = _np.loadtxt(fname, skiprows=0,usecols=(0,1),unpack=True)
        else:
            if not silent: print('Not found.')
            Exception('Longitudinal wake file not found.')
        simul_data.s = spos/100    # Rescaling cm to m
        simul_data.Wll = -wl * 1e12 # V/C the minus sign is due to convention

    elif anal_pl in {'dx','dy'}:
        fname = _jnPth([path,'wakeL.dat'])
        if _os.path.isfile(fname):
            if not silent: print('Calculating Transverse wake from longitudinal wake file:',end='')
            if not silent: print('Data found.')
            spos,wl = _np.loadtxt(fname, skiprows=0,usecols=(0,1),unpack=True)
            simul_data.s = spos/100    # Rescaling cm to m
            wt = -_scy_int.cumtrapz(-wl,x=spos/100,initial=0) # one minus sign due to convention and the other due to Panofsky-Wenzel
            setattr(simul_data, 'W'+anal_pl, wt * 1e12) # V/C/m
        else:
            if not silent: print('File not found.\nLoading transverse wake from transverse wake file.:',end='')
            fname = _jnPth([path,'wakeT.dat'])
            if _os.path.isfile(fname):
                if not silent: print('Data found.\nDepending on the ECHOz2 program version this may lead to inacurate results.')
                spos,wt = _np.loadtxt(fname, skiprows=0,usecols=(0,1),unpack=True)
            else:
                if not silent: print('Not found.')
                Exception('Transverse wake file not found.')
            simul_data.s = spos/100    # Rescaling cm to m
            # there is an error in the integration of echoz2. It is needed to subtract
            # the first value to correct and offset
            # wt = -_scy_int.cumtrapz(-wl,x=spos/100,initial=0)
            setattr(simul_data, 'W'+anal_pl, (wt-wt[0]) * 1e12) # V/C/m  the minus sign is due to convention
    else:
        if not silent: print('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))
        raise Exception('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))

    # loading driving bunch info
    if not silent: print('Loading bunch length from wake file')
    sbun = simul_data.s.copy()
    ds = sbun[1]-sbun[0]
    bunlen = abs(sbun[0] -ds/2) / 5
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

def _load_data_ECHO_rect(simul_data,code,path,anal_pl,silent):
    def _load_dados(fname,mode, bc, code):
        if code == 'echozr':
            len_unit, charge_unit, header = 1e-2, 1e-12, 3 # cm to m, pC to C
            with open(fname) as f:
                f.readline()
                a = f.readline()
            mstep, offset, wid, bunlen = _np.fromstring(a[1:],sep='\t')
            offset = int(offset)
            arbitrary_factor = 2 # I don't know why I have to divide the echozr data by 2;
            y0 = y = mstep*offset / 100
        elif code == 'echo2d':
            len_unit, charge_unit, header = 1, 1, 5
            with open(fname) as f:
                f.readline()
                mstep, offset = _np.fromstring(f.readline(),sep='\t')
                f.readline()
                wid, bunlen   = _np.fromstring(f.readline(),sep='\t')
            offset = int(offset)
            arbitrary_factor = 1 # But I don't have to do this for the echo2d data.
            y0 = y = mstep*offset
            offset = 0  # has only one column of wake
        spos, Wm = _np.loadtxt(fname,skiprows=header,usecols=(0,1+offset),unpack=True)
        mstep *= len_unit
        wid   *= len_unit
        bunlen*= len_unit
        spos  *= len_unit
        Wm    *= -len_unit/charge_unit / arbitrary_factor  # minus sign is due to convention

        Kxm = _np.pi/wid*mode
        if bc == 'elec': Wm  /= _np.sinh(Kxm*y0)*_np.sinh(Kxm*y)
        else:            Wm  /= _np.cosh(Kxm*y0)*_np.cosh(Kxm*y)
        return spos, Wm, mstep, wid, bunlen


    if anal_pl == 'db':
        if not silent: print('Problem: All rectangular geometries does not have symmetry.')
        raise Exception('Problem: All rectangular geometries does not have symmetry.')

    if anal_pl == 'll':          bc = 'magn'
    elif anal_pl in {'dx','dy'}: bc = 'elec'

    if not silent: print('Looking for data files in subfolder {0:s}.'.format(bc))
    pname = _jnPth([path,bc])
    if not _os.path.isdir(pname):
        pname = path
        if code == 'echozr':
            if not silent:
                print('Subfolder not found. It would be better to'+
                      ' create the subfolder and put the files there...')
                print('Looking for files in the current folder:')
        elif code == 'echo2d':
            if not silent: print('Files not found. ')
            raise Exception('Files not found.')

    f_in_dir = _sh.ls(pname).stdout.decode()
    f_match = sorted(_re.findall(FNAME_ECHOZR2D,f_in_dir))
    if not f_match:
        if not silent: print('Files not found.')
        raise Exception('Files not found.')

    if not silent:
        print('Files found.\n I am assuming the simulation was performed '+
              'with {0:s} boundary condition.'.format('electric' if bc == 'elec' else 'magnetic'))
        print('Modes found: ' + ', '.join([m for _,m in f_match]))
        print('Loading data from files')

    spos, W, mode, mesh_size, width, bunlen = [], [], [], [], [], []
    for fn, m in f_match:
        s, Wm, ms, wid, bl = _load_dados(_jnPth([pname,fn]),int(m),bc,code)
        mode.append(int(m))
        spos.append(s)
        W.append(Wm)
        mesh_size.append(ms)
        width.append(wid)
        bunlen.append(bl)

    cond = False
    for i in range(1,len(mode)):
        cond |= len(spos[i]) != len(spos[0])
        cond |= not _np.isclose(mesh_size[i], mesh_size[0], rtol=1e-5, atol=0)
        cond |= not _np.isclose(width[i], width[0], rtol=0, atol=1e-7)
        cond |= not _np.isclose(bunlen[i], bunlen[0], rtol=1e-5, atol=0)
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
            Wll *= 2/width[0]
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
            Wq *= 2/width[0]
            Wq  = -_scy_int.cumtrapz(Wq,x=spos[0],initial=0) # minus sign is due to Panofsky-Wenzel
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
            Wdx *= 2/width[0]
            Wdx  = -_scy_int.cumtrapz(Wdx,x=spos[0],initial=0) # minus sign is due to Panofsky-Wenzel
            simul_data.Wdx = Wdx

    elif anal_pl in {'dx','dy'}:
        pl = 'Vertical' if anal_pl == 'dy' else 'Horizontal'
        if not silent: print('Calculating Dipolar {0:s} Wake from data:'.format(pl))
        Wd, frac =  None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wd = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wd += W[i] * Kxm**2 #only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wd))
        if Wd is None:
            if not silent: print('There is none even mode to calculate Dipolar {0:s} Wake.'.format(pl))
        else:
            if not silent: print('Maximum influence of last mode in the final result is: {0:5.2f}%'.format(frac*100))
            Wd *= 2/width[0]
            Wd  = -_scy_int.cumtrapz(Wd,x=spos[0],initial=0) # minus sign is due to Panofsky-Wenzel
            setattr(simul_data,'W'+anal_pl,Wd)
    else:
        if not silent: print('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))
        raise Exception('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))

def _load_data_ECHOzR(simul_data,path,anal_pl,silent=False):
    _load_data_ECHO_rect(simul_data,'echozr',path,anal_pl,silent)

def _load_data_ECHO2D(simul_data,path,anal_pl,silent=False):

    if not silent: print('Trying to find out the geometry type: ',end='')

    if (os.path.isdir(_jnPth([path,'magn'])) or
        os.path.isdir(_jnPth([path,'elec']))):
        geo_type = 'rectangular'
    elif (os.path.isfile(_jnPth([path,'wakeL_00.txt'])) or
          os.path.isfile(_jnPth([path,'wakeL_01.txt']))):
        geo_type = 'round'
    else:
        if not silent: print('not ok.\n Could not find out the geometry type.')
        raise Exception('Could not find out the geometry type.')
    if not silent: print(geo_type)

    if geo_type == 'rectangular':
        _load_data_ECHO_rect(simul_data,'echo2d',silent)
    else:
        anal_pl_ori = None
        if anal_pl == 'db':
            anal_pl_ori = 'db'
            anal_pl     = 'dy'
            if not silent: print('Even though there is symmetry, I am loading data to the Y plane.')

        if anal_pl=='ll':
            if not silent: print('Loading longitudinal Wake file:',end='')
            fname = _jnPth([path,'wakeL_00.dat'])
            if os.path.isfile(fname):
                if not silent: print('Data found.')
                with open(fname) as f:
                    f.readline()
                    mstep, offset = _np.fromstring(f.readline(),sep='\t')
                    f.readline()
                    _, bunlen   = _np.fromstring(f.readline(),sep='\t')
                spos, Wm = _np.loadtxt(fname,skiprows=5,unpack=True)
                simul_data.s = spos
                simul_data.Wll = -Wm # V/C the minus sign is due to convention
            else:
                if not silent: print('Not found.')
                Exception('Longitudinal wake file not found.')


        elif anal_pl in {'dx','dy'}:
            if not silent: print('Loading Transverse Wake file:',end='')
            fname = _jnPth([path,'wakeL_01.dat'])
            if _os.path.isfile(fname):
                if not silent: print('Data found.')
                with open(fname) as f:
                    f.readline()
                    mstep, offset = _np.fromstring(f.readline(),sep='\t')
                    f.readline()
                    _, bunlen   = _np.fromstring(f.readline(),sep='\t')
                y0 = mstep*(offset+0.5) # transverse wakes are calculated in the middle of the mesh
                spos, Wm = _np.loadtxt(fname,skiprows=5,unpack=True) # m and V/C/m^2
                simul_data.s = spos
                Wdm = -scy_int.cumtrapz(-Wm/y0^2,x=spos,initial=0) # V/C/m the minus sign is due to convention
                setattr(simul_data, 'W'+anal_pl, Wdm)
            else:
                if not silent: print('File not found.')
                Exception('Transverse wake file not found.')
        else:
            if not silent: print('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))
            raise Exception('Plane of analysis {0:s} does not match any of the possible options'.format(anal_pl))

        if anal_pl_ori:
            anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
            if not silent: print('There is symmetry. Copying the data from the '+
                        '{0:s} plane to the {1:s} plane'.format(anal_pl[1].upper(),anal_pl_comp[1].upper()))
            setattr(simul_data, 'W'+anal_pl_comp, getattr(simul_data,'W'+anal_pl).copy())

        a    = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        simul_data.bunlen = bunlen
        simul_data.sbun   = sbun
        simul_data.bun    = _np.exp(-sbun**2/(2*bunlen**2))/(_np.sqrt(2*_np.pi)*bunlen)
        if not silent:
            print('Bunch length of the driving bunch: {0:7.3g} mm'.format(simul_data.bunlen*1e3))
            print('Mesh size used in the simulation:  {0:7.4g} um'.format(mstep*1e6))
            print('Data Loaded.')


CODES    = {'echoz1': _load_data_ECHOz1,
            'echoz2': _load_data_ECHOz2,
            'echo2d': _load_data_ECHO2D,
            'echozr': _load_data_ECHOzR,
            'gdfidl': _load_data_GdfidL,
            'ace3p' : _load_data_ACE3P,
            'cst'   : _load_data_CST
           }

def load_raw_data(simul_data=None, code=None, path=None, anal_pl=None, silent=False):
    if not simul_data: simul_data = EMSimulData()

    if path is None: path = _os.path.abspath('.')


    if not silent: print('#'*60 + '\nLoading Simulation Data')

    #Split the path to try to guess other parameters:
    path_split = set(path.lower().split(_os.path.sep))

    #First try to guess the code used in simulation, if not supplied:
    if code is None:
        if not silent: print('Simulation Code not supplied, trying to guess from path: ', end='')
        code_guess = list(CODES.keys() & path_split)
        if code_guess:    code = code_guess[0]
        else:
            if not silent: print('could not be guessed.')
            if not silent: print('Trying to guess from files in folder: ', end='')
            f_mat = None
            f_in_dir = _sh.ls(path).stdout.decode()
            if len(_re.findall(FNAME_GDFIDL,f_in_dir)):     code = 'gdfidl'
            elif len(_re.findall(FNAME_ECHOZ1,f_in_dir)):   code = 'echoz1'
            elif len(_re.findall(FNAME_ECHOZ2,f_in_dir)):   code = 'echoz2'
            elif len(_re.findall(FNAME_ECHOZR2D,f_in_dir)):
                fol = path
                f_mat = _re.findall(FNAME_ECHOZR2D,f_in_dir)
            elif _os.path.isdir(_jnPth([path,'elec'])):
                fol = _jnPth([path,'elec'])
                f_in_fol = _sh.ls(fol).stdout.decode()
                f_mat = _re.findall(FNAME_ECHOZR2D,f_in_fol)
            elif _os.path.isdir(_jnPth([path,'magn'])):
                fol = _jnPth([path,'magn'])
                f_in_fol = _sh.ls(fol).stdout.decode()
                f_mat = _re.findall(FNAME_ECHOZR2D,f_in_fol)
            else:   raise Exception('Simulation Code was not supplied and could not be guessed.')
            if f_mat is not None:
                if _os.path.isfile(_jnPth([fol,f_mat[0][0]])):
                    with open(_jnPth([fol,f_mat[0][0]])) as f:
                        code = 'echozr' if f.readline().find('[cm]')> 0 else 'echo2d'
                else:  raise Exception('Simulation Code was not supplied and could not be guessed.')
    if not silent: print(code)
    simul_data.code = code

    # Now try to guess the plane of the analysis:
    if anal_pl is None:
        if not silent: print('Plane of Analysis not supplied, trying to guess from path: ', end='')
        anal_pl_guess = list(ANALYSIS_TYPES & path_split)
        if anal_pl_guess:    anal_pl = anal_pl_guess[0]
        else:
            if not silent: print('could not be guessed.')
            if not silent: print('Trying to guess from files in folder and code: ', end='')
            if  code == 'echoz1':   anal_pl = 'll'
            elif code == 'echoz2':  anal_pl = 'dy' if _os.path.isfile('wakeT.dat') else 'll'
            elif code == 'gdfidl':
                f_mat = _re.findall(r"[\w-]+W([YXq]{2})_AT_XY.[0-9]{4}",f_in_dir)
                if len(f_mat) > 0:  anal_pl = 'd'+f_mat[0][0].lower()
                else:               anal_pl = 'll'
            elif code == 'echozr':
                if _os.path.isdir(_jnPth([path,'magn'])):    anal_pl = 'll'
                elif _os.path.isdir(_jnPth([path,'elec'])):  anal_pl = 'dy'
                elif _os.path.isfile('wakeL_01.txt'):
                    w = _np.loadtxt('wakeL_01.txt',skiprows=3,usecols=(1,),unpack=True)
                    if _np.all(w==0):  anal_pl = 'dy'
                    else:              anal_pl = 'll'
                else:  raise Exception('Plane of analysis was not supplied and could not be guessed.')
            elif code == 'echo2d':
                if _os.path.isdir(_jnPth([path,'magn'])):     anal_pl = 'll'
                elif _os.path.isdir(_jnPth([path,'elec'])):   anal_pl = 'dy'
                elif _os.path.isfile('wakeL_00.txt'):         anal_pl = 'll'
                else:                                         anal_pl = 'dy'
            else:    raise Exception('Plane of analysis was not supplied and could not be guessed.')
    if not silent: print(anal_pl)


    CODES[code](simul_data,silent=silent,path=path,anal_pl=anal_pl) # changes in simul_data are made implicitly

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
        VHat *= _np.exp(-1j*w*(spos[0])/c)
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

    for pl in PLANES:
        if not silent: print('Performing FFT on W{0:s}: '.format(pl),end='')
        Wpl = getattr(simul_data,'W'+pl).copy()
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


def calc_impedance_naff(simul_data, pl='ll', s_min = None, win = 1, nr_ff = 20):

    if pl not in PLANES:
        raise Exception('Value of variable pl not accepted. Must be one of these: '+', '.join(PLANES))

    # Extracts Needed Variables
    sigt    = simul_data.bunlen / c  # bunch time-length
    spos    = simul_data.s.copy()
    W    = getattr(simul_data,'W' + pl).copy()

    if W is None or not len(W) or _np.all(W == 0):
        raise Exception('No Data found.')

    if s_min is None: s_min = spos[0]

    inds = spos >= s_min
    spos = spos[inds]
    W    = W[inds]
    dt   = (spos[1]-spos[0])/c
    leng = len(W) + 1 - (len(W)%6)
    spos = spos[-leng:]
    W    = W[-leng:]
    if win == 1/2:
        W  *= _np.hanning(2*spos.shape[0])[spos.shape[0]:]
        tu,a = _naff.naff_general(W,use_win=0, is_real=False, nr_ff=nr_ff)
    elif isinstance(win,int):
        tu,a = _naff.naff_general(W,use_win=win, is_real=False, nr_ff=nr_ff)
    else:
        raise Exception('Win must be 1/2 for half-hanning window or an integer for other windows(0 --> no window).')

    freq = tu/dt
    w    = 2*_np.pi*freq
    # Longitudinal position shift to match center of the bunch with zero z:
    a   *= _np.exp(-1j*w*(spos[0])/c)
    # Deconvolve the Transform with a gaussian bunch:
    a   /= _np.exp(-(w*sigt)**2/2)

    # Must multiply by the vector length due to difference in the meaning of the
    # amplitune in the NAFF transform and the fourier transform
    Z    = a*dt*leng

    if pl =='ll':
        # I have to take the conjugate of the fft because:
        #fftt == \int exp(-i*2pi*f*t/n) G(t) dt
        #while impedance, according to Chao and Ng, is given by:
        #Z == \int exp(i*2pi*f*t/n) G(t) dt
        Z = Z.conj()
    else:
        #the Transverse impedance, according to Chao and Ng, is given by:
        #Z == i\int exp(i*2pi*f*t/n) G(t) dt
        Z = 1j*Z.conj()

    return freq, Z, leng

def plot_wakes(simul_data,save_figs=False,pth2sv=None,show=False,pls=None):

    sbun = simul_data.sbun
    sigs = simul_data.bunlen
    spos = simul_data.s
    if pls is None: pls = PLANES
    for pl in pls:
        wake = getattr(simul_data,'W'+pl)*1e-12 # V/C -> V/pC
        if wake is None or _np.all(wake==0): continue
        max_wake = wake[_np.abs(wake).argmax()]
        bunchshape = simul_data.bun * (max_wake/simul_data.bun.max())

        f,axs = _plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(14,6))
        ax = axs[0]
        b = ax.get_position()
        b.x0, b.x1 = 0.05, 0.45
        ax.set_position(b)
        ax.plot(sbun*1000,bunchshape,'b',linewidth=2,label='Bunch Shape')
        ax.plot(spos*1000,wake,'r',linewidth=2,label='Wake Potential')
        ax.grid(True)
        ax.set_ylabel(WAKE_YLABELS[pl],fontsize=13)
        ax.set_xlim([spos[0]*1000, 8000*sigs])
        ax.set_ylim([wake.min()*1.1, wake.max()*1.1])
        ax.legend(loc='best')

        ax = axs[1]
        b = ax.get_position()
        b.x0, b.x1 = 0.45, 0.95
        ax.set_position(b)
        ax.plot(spos*1000,wake,'r',linewidth=2)
        ax.grid(True)
        tit = ax.set_title(TITLES[pl],fontsize=13)
        tit.set_x(0.1)
        xl  = ax.set_xlabel('s [mm]',fontsize=13)
        xl.set_x(0.1)
        ax.set_xlim([8000*sigs,spos[-1]*1000])
        ax.set_ylim([wake.min()*1.1, wake.max()*1.1])

        if save_figs: f.savefig(_jnPth((pth2sv,'W'+pl+'.svg')))
    if show: _plt.show()

def plot_impedances(simul_data,save_figs=False,pth2sv=None,show=False,pls=None):

    freq = simul_data.freq
    if pls is None: pls = PLANES
    for pl in pls:
        Z = getattr(simul_data,'Z'+pl)
        if Z is None or _np.all(Z==0): continue

        _plt.figure()
        _plt.plot(freq/1e9,Z.real,'r',linewidth=2,label='Re')
        _plt.plot(freq/1e9,Z.imag,'b--',linewidth=2,label='Im')
        _plt.xlabel('Frequency [GHz]',fontsize=13)
        _plt.grid(True)
        _plt.title (TITLES[pl],fontsize=13)
        _plt.ylabel(IMPS_YLABELS[pl],fontsize=13)
        _plt.legend (loc='best')
        _plt.xlim(freq[[0,-1]]/1e9)
        if save_figs: _plt.savefig(_jnPth((pth2sv,'Z'+pl+'.svg')))
    if show: _plt.show()

def plot_losskick_factors(simul_data,save_figs=False,pth2sv=None,show=False,pls=None):
    # Extracts and Initialize Needed Variables:
    _si.nom_cur = 500e-3
    sigvec = _np.array([2.65, 5, 8, 10, 15],dtype=float)*1e-3  # bunch length scenarios
    Ivec   = _np.linspace(10e-3,_si.nom_cur,num=50) # current scenarios


    bunlen = simul_data.bunlen
    sigi   = _np.linspace(bunlen,18e-3,num=50)
    fill_pat = _np.array([1,864,864/2,864/4],dtype=int)
    if pls is None: pls = PLANES
    pls2 = []
    for pl in pls:
        W = getattr(simul_data,'W'+pl)
        if W is None or _np.all(W==0): continue
        Z = getattr(simul_data,'Z'+pl)
        if Z is None or _np.all(Z==0): continue
        pls2.append(pl)

    for pl in pls2:
        if pl == 'll':
            f,axs = _plt.subplots(nrows=1, ncols=2, figsize=(12,6),gridspec_kw=dict(left=0.08,right=0.97))
            ax = axs[0]
            fname = 'Loss_factor'
            for i in range(fill_pat.shape[0]):
                kZi = _np.zeros(sigi.shape[0])
                for j in range(sigi.shape[0]):
                    kZi[j] = simul_data.klossZ(sigma=sigi[j],n=fill_pat[i]) * 1e-12 #V/pC
                ax.semilogy(sigi * 1e3, kZi * 1e3, 'o',markersize=4,label=r'$n = {0:03d}$'.format(fill_pat[i]))
                if not i: kZ = kZi[0]
            # Calculates klossW
            kW    = simul_data.klossW() * 1e-12
            # Print loss factor calculated in both ways
            ax.semilogy(bunlen * 1e3, kW * 1e3, '*',markersize=7,color=[1, 0, 0],label=r'$K_L^W$')
            ax.set_title('Loss Factor for $n$ equally spaced bunches.')
            ax.set_xlabel(r'$\sigma$ [mm]')
            ax.set_ylabel(r'$K_L$ [mV/pC]')
            ax.legend(loc='best')
            ax.grid(True)
            ax.annotate(r'$K_L^W = {0:5.2f}$ mV/pC'.format(kW*1e3),xy=(bunlen*1.1e3, kW*1e3),fontsize=12)

            ax = axs[1]
            kZvec = _np.zeros(sigvec.shape[0])
            labels = []
            for i in range(sigvec.shape[0]):
                kZvec[i] = simul_data.klossZ(sigma=sigvec[i], n=_si.harm_num) #V/C
                labels.append(r'$\sigma = {0:05.2f}$ mm'.format(sigvec[i]*1e3))
            Plossvec = kZvec[None,:] * Ivec[:,None]**2 * _si.T0/_si.harm_num
            ax.semilogy(Ivec*1e3, Plossvec,markersize=4)
            ax.set_title('Power Loss for ${0:d}$ equally spaced bunches.'.format(_si.harm_num))
            ax.set_xlabel(r'$I_{{avg}}$ [mA]')
            ax.set_ylabel(r'Power [W]')
            ax.legend(labels,loc='best')
            ax.grid(True)
        else:
            f  = _plt.figure(figsize=(6,6))
            ax = _plt.axes()
            fname = 'Kck'+pl+'_factor'
            for i in range(fill_pat.shape[0]):
                kZi = _np.zeros(sigi.shape[0])
                for j in range(sigi.shape[0]):
                    kZi[j] = simul_data.kick_factorZ(pl=pl,sigma=sigi[j],n=fill_pat[i]) * 1e-12 #V/pC/m
                ax.plot(sigi*1e3, kZi, 'o',markersize=4,label=r'n = {0:03d}'.format(fill_pat[i]))
                if not i: kickZ = kZi[0]
            # Calculates kickW:
            kickW = simul_data.kick_factorW(pl=pl) * 1e-12
            # Print loss factor calculated in both ways
            ax.plot(bunlen*1e3, kickW, '*',markersize=7,color=[1, 0, 0],label=r'$K_{{{0:s}_{1:s}}}^W$'.format(pl[0].upper(),pl[1]))
            ax.set_title('Kick Factor for $n$ equally spaced bunches.')
            ax.set_xlabel(r'$\sigma$ [mm]',fontsize=13)
            ax.set_ylabel(r'$K_{{{0:s}_{1:s}}}$ [V/pC/m]'.format(pl[0].upper(),pl[1]),fontsize=13)
            ax.legend(loc='best')
            ax.grid(True)
            stri = r'$K_{{{0:s}_{1:s}}}^W = {2:5.2f}$ V/pC/m'.format(pl[0].upper(),pl[1],kickW)
            ax.annotate(stri,xy=(bunlen*1.1e3, kickW),fontsize=13)
        if save_figs: _plt.savefig(_jnPth((pth2sv,fname+'.svg')))
    if show: _plt.show()

def show_now():
    _plt.show()

def save_processed_data(simul_data,silent=False,pth2sv=None):

    if not silent: print('#'*60 + '\nSaving Processed data:')
    spos = simul_data.s
    freq = simul_data.freq

    if pth2sv is None:
        if not silent: print('Saving in the same folder of the raw data')
        pth2sv = os.path.abspath('.')
    elif type(pth2sv) is str:
        if not silent: print('Saving to subfolder: ' + pth2sv)
        if not _os.path.isdir(pth2sv):
            if not silent: print('Folder does not exist. Creating it...')
            _os.mkdir(pth2sv)
    else:
        if not silent: print('pth2sv must be a string or None object')
        raise Exception('pth2sv must be a string or None')

    #Save wakes
    for par in PLANES:
        unit = 'V/C' if par == 'll' else 'V/C/m'
        header = '{0:30s} {1:30s}'.format('s [m]', 'W{0:s} [{1:s}]'.format(par,unit))
        fname = _jnPth([pth2sv,'W'+par+'.gz'])
        wake  = getattr(simul_data,'W'+par)
        if wake is None or _np.all(wake == 0): continue
        if not silent: print('Saving W'+ par + ' data to .gz file')
        _np.savetxt(fname,_np.array([spos,wake]).transpose(),
                    fmt=['%30.16g','%30.16g'], header=header)

    #Save Impedances
    for par in PLANES:
        unit = 'Ohm' if par == 'll' else 'Ohm/m'
        header = '{0:30s} {1:30s} {2:30s}'.format('Frequency [GHz]',
                                            'ReZ{0:s} [{1:s}]'.format(par,unit),
                                            'ImZ{0:s} [{1:s}]'.format(par,unit))
        fname = _jnPth([pth2sv,'Z'+par+'.gz'])
        Z  = getattr(simul_data,'Z'+par)
        if Z is None or _np.all(Z == 0): continue
        if not silent: print('Saving Z'+ par + ' data to .gz file')
        _np.savetxt(fname,_np.array([freq/1e9,Z.real,Z.imag]).transpose(),
                    fmt=['%30.16g','%30.16g','%30.16g'], header=header)

    if not silent: print('Saving the Complete EMSimulData structure to a .pickle file.')
    with _gzip.open(_jnPth((pth2sv,DEFAULT_FNAME_SAVE)), 'wb') as f:
        _pickle.dump(simul_data,f,_pickle.HIGHEST_PROTOCOL)

    if not silent: print('All Data Saved\n' + '#'*60)

def load_processed_data(filename):
    with _gzip.open(filename,'rb') as fh:
        simul_data = _pickle.load(fh)
    return simul_data

def create_make_fig_file(path = None):
    if path is None: path = os.path.abspath('.')
    fname = _jnPth([path,'create_figs.py'])
    analysis  = '#!/usr/bin/env python3\n\n'
    analysis += 'import os\n'
    analysis += 'import pycolleff.process_wakes as ems\n\n'
    analysis += 'opts = dict(save_figs=False,show=False)\n'
    analysis += 'path = os.path.abspath(__file__).rpartition(os.path.sep)[0]\n'
    analysis += "file_name = os.path.sep.join([path,'{0:s}'])\n".format(DEFAULT_FNAME_SAVE)
    analysis += 'simul_data = ems.load_processed_data(file_name)\n'
    analysis += 'ems.plot_wakes(simul_data,**opts)\n'
    analysis += 'ems.plot_impedances(simul_data,**opts)\n'
    analysis += 'ems.plot_losskick_factors(simul_data,**opts)\n'
    analysis += 'ems.show_now()\n'
    with open(fname,'w') as f:
        f.writelines(analysis)
    _sh.chmod('+x',fname)
