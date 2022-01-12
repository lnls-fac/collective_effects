#!/usr/bin/env python3

import os as _os
import re as _re
import sh as _sh
import gzip as _gzip
import pickle as _pickle
import logging as _log
import json as _json

import numpy as _np
from scipy import integrate as _scy_int
import matplotlib.pyplot as _plt
from matplotlib import rc as _rc

from mathphys.constants import light_speed as c

from . import sirius as _sirius

try:
    from pyaccel import naff as _naff
    bool_pyaccel = True
except Exception:
    bool_pyaccel = False

_rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# for Palatino and other serif fonts use:
# _rc('font',**{'family':'serif','serif':['Palatino']})
_rc('text', usetex=True)

# # Troughout the code I am assuming:
# s positive means particle behind source -->  Wl, Wt = 0 s < 0
# Wl(s) = -c/Q * int El(ct-s,t) dt
# Wx(s) = - int_-inf^s dWl/dx ds'
# Zl =   int exp(i*w*s) Wl(s) ds
# Zx = i*int exp(i*w*s) Wx(s) ds

_jnpth = _os.path.sep.join
_si = _sirius.create_ring()

DEFAULT_FNAME_SAVE = 'SimulData.pickle'
# the older .dat files are not treated:
FNAME_ECHOZR2D = r"(wakeL_([0-9]{2}).txt)"
FNAME_GDFIDL = r"[\w-]+W[YXq]{1}_AT_XY.[0-9]{4}"

ANALYSIS_TYPES = {
    'dx',  # horizontal impedance
    'dy',  # vertical impedance
    'db',  # both planes are symmetric
    'll'   # longitudinal and transverse quadrupolar impedances
    }

PLANES = ('ll', 'dx', 'dy', 'qx', 'qy')
TITLES = {
    'll': 'Longitudinal',
    'dx': 'Dipolar Horizontal',
    'dy': 'Dipolar Vertical',
    'qx': 'Quadrupolar Horizontal',
    'qy': 'Quadrupolar Vertical'}
WAKE_YLABELS = {
    'll': r'$W_l$ [V/pC]',
    'dx': r'$W_{{D_x}}$ [V/pC/m]',
    'dy': r'$W_{{D_y}}$ [V/pC/m]',
    'qx': r'$W_{{Q_x}}$ [V/pC/m]',
    'qy': r'$W_{{Q_y}}$ [V/pC/m]'}
IMPS_YLABELS = {
    'll': r'$Z_l$ [$\Omega$]',
    'dx': r'$Z_{{D_x}}$ [$\Omega$/m]',
    'dy': r'$Z_{{D_y}}$ [$\Omega$/m]',
    'qx': r'$Z_{{Q_x}}$ [$\Omega$/m]',
    'qy': r'$Z_{{Q_y}}$ [$\Omega$/m]'}


def load_raw_data(
        simul_data=None, code=None, path=None, anal_pl=None, silent=False):
    """."""
    if not simul_data:
        simul_data = EMSimulData()

    _log.basicConfig(level=_log.CRITICAL if silent else _log.INFO)

    if path is None:
        path = _os.path.abspath('.')

    _log.info('#'*60 + '\nLoading Simulation Data')

    # First try to guess the code used in simulation, if not supplied:
    if code is None:
        _log.info('Simulation Code not supplied.')
        code = _get_code(path)
    _log.info(code)
    simul_data.code = code

    # Now try to guess the plane of the analysis:
    if anal_pl is None:
        _log.info('Plane of Analysis not supplied.')
        anal_pl = _get_plane_of_analysis(path, code)
    _log.info(anal_pl)

    # changes in simul_data are made implicitly
    CODES[code](simul_data, path=path, anal_pl=anal_pl)

    _log.info('#'*60+'\n')
    return simul_data


# ########################## Auxiliary Methods ##########################

def _get_code(path):
    # Split the path to try to guess other parameters:

    _log.info('Trying to guess from path: ')
    path_split = set(path.lower().split(_os.path.sep))
    code_guess = list(CODES.keys() & path_split)
    if code_guess:
        return code_guess[0]
    _log.info('could not be guessed by path.')

    _log.info('Trying to guess from files in folder: ')
    f_in_dir = ' '.join(_os.listdir(path))
    if len(_re.findall(FNAME_GDFIDL, f_in_dir)):
        return 'gdfidl'

    if len(_re.findall(FNAME_ECHOZ1, f_in_dir)):
        return 'echoz1'

    if len(_re.findall(FNAME_ECHOZ2, f_in_dir)):
        return 'echoz2'

    if len(_re.findall(FNAME_ECHO3D, f_in_dir)):
        return 'echo3d'

    f_mat = None
    if len(_re.findall(FNAME_ECHOZR2D, f_in_dir)):
        fol = path
        f_mat = _re.findall(FNAME_ECHOZR2D, f_in_dir)
    elif _os.path.isdir(_jnpth([path, 'elec'])):
        fol = _jnpth([path, 'elec'])
        f_in_fol = ' '.join(_os.listdir(fol))
        f_mat = _re.findall(FNAME_ECHOZR2D, f_in_fol)
    elif _os.path.isdir(_jnpth([path, 'magn'])):
        fol = _jnpth([path, 'magn'])
        f_in_fol = ' '.join(_os.listdir(fol))
        f_mat = _re.findall(FNAME_ECHOZR2D, f_in_fol)

    if f_mat is not None and _os.path.isfile(_jnpth([fol, f_mat[0][0]])):
        with open(_jnpth([fol, f_mat[0][0]])) as f:
            code = 'echozr'
            if f.readline().find('[cm]') <= 0:
                code = 'echo2d'
            return code

    msg = 'Simulation Code was not supplied and could not be guessed.'
    _log.info(msg)
    raise Exception(msg)


def _get_plane_of_analysis(path, code):
    # Split the path to try to guess other parameters:
    path_split = set(path.lower().split(_os.path.sep))

    _log.info('Trying to guess from path: ')
    anal_pl_guess = list(ANALYSIS_TYPES & path_split)
    if anal_pl_guess:
        return anal_pl_guess[0]
    _log.info('could not be guessed by path.')

    _log.info('Trying to guess from files in folder and code: ')
    if code == 'echoz1':
        return 'll'

    if code == 'echoz2':
        return 'dy' if _os.path.isfile('wakeT.dat') else 'll'

    if code == 'gdfidl':
        f_in_dir = ' '.join(_os.listdir(path))
        f_mat = _re.findall(
            r"[\w-]+W([YXq]{2})_AT_XY.[0-9]{4}", f_in_dir)
        if len(f_mat) > 0:
            # f_mat = _re.findall(
            #     r"[\w-]+W([YX]{1})_AT_XY.[0-9]{4}",f_in_dir)
            # countx = [x for x in f_mat if x=='X']
            # county = [y for y in f_mat if y=='Y']
            # anal_pl = 'dy' if len(county) >= len(county) else 'dx'
            anal_pl = 'd'+f_mat[0][0].lower()
        else:
            anal_pl = 'll'
        return anal_pl

    if code == 'echozr':
        if _os.path.isdir(_jnpth([path, 'magn'])):
            return 'll'
        if _os.path.isdir(_jnpth([path, 'elec'])):
            return 'dy'
        if _os.path.isfile('wakeL_01.txt'):
            w = _np.loadtxt(
                'wakeL_01.txt', skiprows=3, usecols=(1, ), unpack=True)
            anal_pl = 'll'
            if _np.allclose(w, 0, atol=0):
                anal_pl = 'dy'
            return anal_pl
        msg = 'Plane of analysis was not supplied '
        msg += 'and could not be guessed.'
        _log.info(msg)
        raise Exception(msg)

    if code == 'echo2d':
        if _os.path.isdir(_jnpth([path, 'magn'])):
            return 'll'
        if _os.path.isdir(_jnpth([path, 'elec'])):
            return 'dy'

        anal_pl = 'dy'
        if _os.path.isfile('wakeL_00.txt'):
            anal_pl = 'll'
        return anal_pl

    if code == 'echo3d':
        if _os.path.isdir(_jnpth([path, 'magn'])):
            return 'll'
        if _os.path.isdir(_jnpth([path, 'elec'])):
            return 'dy'

    msg = 'Plane of analysis was not supplied '
    msg += 'and could not be guessed.'
    _log.info(msg)
    raise Exception(msg)


def _ACE3P_load_data(simpar):
    raise NotImplementedError('This function was not tested yet.')
    nsigmas = 5
    headerL = 3
    if wdir.startswith(tardir):
        cpfile = False
    else:
        cpfile = True

    wakepath = _jnpth([wdir, 'wakefield.out'])
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile:
        _sh.cp(wakepath, tardir)

    spos = loadres[:, 0]
    # I know this is correct for ECHO (2015/08/27):
    if m == 0:
        wake = -loadres[:, 1]
    else:
        wake = loadres[:, 1]

    spos = spos - nsigmas * bunlen  # Performs displacement over s axis
    return spos, wake


def _CST_load_data(simpar):
    raise NotImplementedError('This function was not tested yet.')
    headerL = 2
    if wdir.startswith(tardir):
        cpfile = False
    else:
        cpfile = True

    wakepath = _jnpth([wdir, 'wake.txt'])
    loadres = _np.loadtxt(wakepath, skiprows=headerL)
    if cpfile:
        _sh.cp(wakepath, tardir)

    spos = loadres[:, 0]
    wake = loadres[:, 1]
    # Adjust s-axis (rescale or shift)
    spos = spos/1000  # Rescaling mm to m
    if m > 0:
        wake = -wake
    return spos, wake


def _ECHO2D_load_data(simul_data, path, anal_pl):

    _log.info('Trying to find out the geometry type: ')

    if (_os.path.isdir(_jnpth([path, 'magn'])) or
            _os.path.isdir(_jnpth([path, 'elec']))):
        geo_type = 'rectangular'
    elif (_os.path.isfile(_jnpth([path, 'wakeL_00.txt'])) or
            _os.path.isfile(_jnpth([path, 'wakeL_01.txt']))):
        geo_type = 'round'
    else:
        msg = 'Could not find out the geometry type.'
        _log.info(msg)
        raise Exception(msg)
    _log.info(geo_type)

    if geo_type == 'rectangular':
        _ECHO_rect_load_data(simul_data, 'echo2d', path, anal_pl)
    else:
        anal_pl_ori = None
        if anal_pl == 'db':
            anal_pl_ori = 'db'
            anal_pl = 'dy'
            _log.info(
                'Even though there is symmetry, '
                'I am loading data to the Y plane.')

        if anal_pl == 'll':
            _log.info('Loading longitudinal Wake file:')
            fname = _jnpth([path, 'wakeL_00.txt'])
            if _os.path.isfile(fname):
                _log.info('Data found.')
                with open(fname) as f:
                    f.readline()
                    mstep, offset = _np.fromstring(f.readline(), sep='\t')
                    f.readline()
                    _, bunlen = _np.fromstring(f.readline(), sep='\t')
                spos, Wm = _np.loadtxt(fname, skiprows=6, unpack=True)
                simul_data.s = spos
                simul_data.Wll = -Wm  # V/C the minus sign is due to convention
            else:
                _log.info('Not found.')
                Exception('Longitudinal wake file not found.')
        elif anal_pl in {'dx', 'dy'}:
            _log.info('Loading Transverse Wake file:')
            fname = _jnpth([path, 'wakeL_01.txt'])
            if _os.path.isfile(fname):
                _log.info('Data found.')
                with open(fname) as f:
                    f.readline()
                    mstep, offset = _np.fromstring(f.readline(), sep='\t')
                    f.readline()
                    _, bunlen = _np.fromstring(f.readline(), sep='\t')
                # transverse wakes are calculated in the middle of the mesh
                y0 = mstep*(offset+0.5)
                # m and V/C/m^2
                spos, Wm = _np.loadtxt(fname, skiprows=6, unpack=True)
                simul_data.s = spos
                # V/C/m the minus sign is due to convention
                Wdm = -_scy_int.cumtrapz(-Wm/(y0*y0), x=spos, initial=0)
                setattr(simul_data, 'W'+anal_pl, Wdm)
            else:
                _log.info('File not found.')
                Exception('Transverse wake file not found.')
        else:
            _log.info(
                f'Plane of analysis {anal_pl:s} does not match '
                'any of the possible options')
            raise Exception(
                f'Plane of analysis {anal_pl:s} does not match '
                'any of the possible options')

        if anal_pl_ori:
            anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
            _log.info(
                'There is symmetry. Copying the data from the ' +
                f'{anal_pl[1].upper():s} plane to' +
                f' the {anal_pl_comp[1].upper():s} plane')
            setattr(
                simul_data, 'W'+anal_pl_comp,
                getattr(simul_data, 'W'+anal_pl).copy())

        a = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        simul_data.bunlen = bunlen
        simul_data.sbun = sbun
        simul_data.bun = _np.exp(-sbun**2/(2*bunlen**2))
        simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen
        _log.info(
            'Bunch length of the driving bunch: ' +
            f'{simul_data.bunlen*1e3:7.3g} mm')
        _log.info(f'Mesh size used in the simulation: {mstep*1e6:7.4g} um')
        _log.info('Data Loaded.')


def _GdfidL_load_data(simul_data, path, anal_pl):
    # list all the files that match the name pattern for wakefields
    f_in_dir = ' '.join(_os.listdir(path))
    f_match = _re.findall(FNAME_GDFIDL, f_in_dir)

    if anal_pl == 'll':
        if not f_match:
            msg = 'No files found for longitudinal analysis.'
            _log.info(msg)
            raise Exception(msg)

        # Load longitudinal Wake
        spos, wake, sbun, bun, bunlen, xd, yd = _GdfidL_get_longitudinal_info(
            path, f_match, pl='ll')
        simul_data.Wll = wake
        simul_data.s = spos
        simul_data.bun = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        # And quadrupolar Wakes, if existent:
        _log.info('Loading Horizontal Quadrupolar Wake file:')
        wake = _GdfidL_get_transversal_info(
            path, f_match, pl='qx')  # V/C/m
        if wake is not None:
            simul_data.Wqx = wake

        _log.info('Loading Vertical Quadrupolar Wake file:')
        wake = _GdfidL_get_transversal_info(
            path, f_match, pl='qy')  # V/C/m
        if wake is not None:
            simul_data.Wqy = wake
        _log.info('Longitudinal Data Loaded.')
        return

    anal_pl_ori = None
    if anal_pl == 'db':
        anal_pl_ori = 'db'
        if not f_match:
            anal_pl = 'dx' if _os.path.isdir('dxdpl') else 'dy'
        else:
            cond = [f for f in f_match if f.find('WX_AT_XY') >= 0]
            anal_pl = 'dx' if cond else 'dy'
        _log.info(
            'There is symmetry y=x, calculation performed in the ' +
            anal_pl[1].upper() + ' plane.')

    if anal_pl not in {'dx', 'dy'}:
        _log.info(
            f'Plane of analysis {anal_pl:s} does not match any of '
            'the possible options')
        raise Exception(
            f'Plane of analysis {anal_pl:s} does not match any of '
            'the possible options')

    elec_symm = False
    if not f_match:
        _log.info('There is no wake files in this folder.')
        elec_fol = _jnpth([path, 'elec'])
        if _os.path.isdir(elec_fol):
            _log.info(
                ' I found a folder named "elec". I will assume the '
                'simulation has this symmetry.')
            f_in_dir = ' '.join(_os.listdir(elec_fol))
            f_match = _re.findall(FNAME_GDFIDL, f_in_dir)
            elec_symm = True

    if not f_match:
        _log.info(' I will assume there is no symmetry.')

        spos, wake, sbun, bun, bunlen, xd, yd = [], [], [], [], [], [], []
        for sub_fol in ['dpl', 'dmi']:
            ext_path = _jnpth([path, anal_pl+sub_fol])
            _log.info('Looking for '+anal_pl+sub_fol+' subfolder:')
            if not _os.path.isdir(ext_path):
                _log.info(
                    'For non-symmetric structures, there must '
                    f'be subfolders {anal_pl:s}dpl {anal_pl:s}dmi '
                    'with the data')
                raise Exception('Files not found')
            # list all the files that match the pattern
            f_in_dir = ' '.join(_os.listdir(ext_path))
            f_match = _re.findall(FNAME_GDFIDL, f_in_dir)
            if not f_match:
                msg = 'No files found for transverse analysis.'
                _log.info(msg)
                raise Exception(msg)

            sp, _, sb, bn, bnln, xdi, ydi = _GdfidL_get_longitudinal_info(
                ext_path, f_match, pl=anal_pl)
            spos.append(sp)
            bun.append(bn)
            sbun.append(sb)
            bunlen.append(bnln)
            xd.append(xdi)
            yd.append(ydi)

            _log.info('Loading {0:s} Dipolar Wake file:'.format(
                'Horizontal' if anal_pl == 'dx' else 'Vertical'))
            # V/C
            wk = _GdfidL_get_transversal_info(
                ext_path, f_match, pl=anal_pl)
            if wk is not None:
                wake.append(wk)
            else:
                _log.info(
                    'Actually there is something wrong, '
                    'these wake files should be here.')
                raise Exception(
                    'Transverse {0:s} dipolar wake files not found'.format(
                        'Horizontal' if anal_pl == 'dx' else 'Vertical'))

        # If the simulation is not ready yet the lenghts may differ.
        # This line is used to truncate the longer wake in the length of
        # the shorter:
        l1 = min(len(spos[0]), len(spos[1]))

        delta = xd if anal_pl == 'dx' else yd
        ndel = yd if anal_pl == 'dx' else xd
        if not (_np.allclose(spos[0][:l1], spos[1][:l1], atol=0) and
                _np.allclose(sbun[0], sbun[1], atol=0) and
                _np.allclose(bun[0], bun[1], atol=0) and
                _np.allclose(ndel[0], ndel[1], atol=0)):
            msg = 'There is a mismatch between the paramenters of the'
            msg += f'simulation in the {anal_pl:s}dpl and {anal_pl:s}dmi '
            msg += 'folders.'
            _log.info(msg)
            raise Exception(msg)
        simul_data.s = spos[0][:l1]
        simul_data.bun = bun[0]
        simul_data.sbun = sbun[0]
        simul_data.bunlen = bunlen[0]
        setattr(
            simul_data, 'W'+anal_pl,
            (wake[0][:l1]-wake[1][:l1])/(delta[0]-delta[1]))  # V/C/m
    else:
        if elec_symm:
            path = elec_fol
        spos, wake, sbun, bun, bunlen, xd, yd = \
            _GdfidL_get_longitudinal_info(
                path, f_match, pl=anal_pl)
        simul_data.s = spos
        simul_data.bun = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        _log.info('Loading {0:s} Dipolar Wake file:'.format(
            'Horizontal' if anal_pl == 'dx' else 'Vertical'))
        wake = _GdfidL_get_transversal_info(
            path, f_match, pl=anal_pl)  # V/C
        if wake is not None:
            delta = xd if anal_pl == 'dx' else yd
            delta *= 2 if elec_symm else 1
            setattr(simul_data, 'W'+anal_pl, wake/delta)  # V/C/m
        else:
            _log.info(
                'Actually there is something wrong, '
                'these wake files should be here.')
            raise Exception(
                'Transverse {0:s} dipolar wake files not found'.format(
                    'Horizontal' if anal_pl == 'dx' else 'Vertical'))
    _log.info('Transverse Data Loaded.')

    if anal_pl_ori:
        anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
        _log.info(
            'There is symmetry. Copying the data from the ' +
            f'{anal_pl[1].upper():s} plane to the ' +
            f'{anal_pl_comp[1].upper():s} plane')
        setattr(
            simul_data, 'W'+anal_pl_comp,
            getattr(simul_data, 'W'+anal_pl).copy())


CODES = {
    'echoz1': _ECHOz1_load_data,
    'echoz2': _ECHOz2_load_data,
    'echo2d': _ECHO2D_load_data,
    'echo3d': _ECHO3D_load_data,
    'echozr': _ECHOzR_load_data,
    'gdfidl': _GdfidL_load_data,
    'ace3p': _ACE3P_load_data,
    'cst': _CST_load_data}


def _GdfidL_load_dados_info(filename):
    dados, info = [], []
    with open(filename) as fh:
        data = fh.read()
    for line in data.splitlines():
        if not line.startswith((' #', ' %', ' $')):
            dados.append(line)
        else:
            info.append(line)
    return dados, info


def _GdfidL_get_charge(info):
    for line in info:
        if line.find('total charge') >= 0:
            lin = line.split(',')[1]
            charge = float(_re.findall(r'[-+]?\d+\.?\d+[eE]?[-+]?\d+', lin)[0])
            break
    return charge


def _GdfidL_get_integration_path(info):
    for line in info:
        if line.find('subtitle=') >= 0:
            x, y = (float(val) for val in _re.findall(
                r'[-+]?\d+\.?\d+[eE]?[-+]?\d+', line))
            break
    return x, y


def _GdfidL_get_longitudinal_info(path, filelist, pl='ll'):
    _log.info('Loading longitunal Wake file:')
    fn = [f for f in filelist if f.find('Wq_AT_XY') >= 0]
    if not fn:
        msg = 'No longitudinal wake file found. It is needed to have one'
        _log.info(msg)
        raise Exception(msg)
    if len(fn) > 1:
        msg = 'More than one longitudinal wake file found. Only 1 is allowed'
        _log.info(msg)
        raise Exception(msg)

    dados, info = _GdfidL_load_dados_info(_jnpth([path, fn[0]]))
    charge = _GdfidL_get_charge(info)
    xd, yd = _GdfidL_get_integration_path(info)
    spos, wake = _np.loadtxt(dados, unpack=True)  # dados is a list of strings
    _log.info(f'Charge of the driving bunch: {charge*1e12:5.3g} pC')
    if pl == 'll' and (abs(xd) > 1e-10 or abs(yd) > 1e-10):
        _log.info(
            'Driving particle not in the origin. '
            'Are you sure this is what you want?')
    elif pl != 'll' and abs(xd) < 1e-10 and abs(yd) < 1e-10:
        _log.info(
            'The driving bunch is too close to origin. '
            'Are you sure this is what you want?')

    a = _np.argmin(_np.diff(spos)) + 1
    sbun = spos[a:]
    bun = wake[a:]*charge/_np.trapz(wake[a:], x=sbun)  # C
    wake = -wake[:a]/charge  # V/C (minus sign because of convention)
    spos = spos[:a]  # m
    bunlen = -sbun[0]/6  # gdfidl uses a bunch with 6-sigma
    _log.info(f'Bunch length of the driving bunch: {bunlen*1e3:7.4g} mm')
    return spos, wake, sbun, bun, bunlen, xd, yd


def _GdfidL_get_transversal_info(path, filelist, pl='qx'):
    stri = 'W{0:s}_AT_XY'.format(pl[1].upper())
    fn = [f for f in filelist if f.find(stri) >= 0]
    if not fn:
        _log.info(f'No W{pl:s} wake file found. Skipping to next')
        return None
    _log.info(f"{len(fn):2d} W{pl:s} wake file found: {', '.join(fn):s}")

    dados, info = _GdfidL_load_dados_info(_jnpth([path, fn[0]]))
    charge = _GdfidL_get_charge(info)
    if pl[1] == 'x':
        delta1, _ = _GdfidL_get_integration_path(info)
    else:
        _, delta1 = _GdfidL_get_integration_path(info)
    _, wake1 = _np.loadtxt(dados, unpack=True)
    _log.info(f'Integration path at {pl[1]:s} = {delta1*1e6:8.4g} um ')

    wake = wake1/delta1/charge  # V/C/m
    if len(fn) > 1:
        dados, info = _GdfidL_load_dados_info(_jnpth([path, fn[1]]))
        if pl[1] == 'x':
            delta2, _ = _GdfidL_get_integration_path(info)
        else:
            _, delta2 = _GdfidL_get_integration_path(info)
        _, wake2 = _np.loadtxt(dados, unpack=True)
        _log.info(f'and {delta2*1e6:8.4g} um')
        if pl[0] == 'd':
            wake = (wake1/delta1 - wake2/delta2)
            wake /= (1/delta1-1/delta2) * charge  # V/C
        else:
            wake = (wake1 - wake2)/(delta1-delta2)/charge  # V/C/m
    else:
        _log.info('')
    return wake


def _ECHO_rect_load_data(simul_data, code, path, anal_pl):
    def _load_dados(fname, mode, bc, code):
        if code == 'echozr':
            len_unit, charge_unit, header = 1e-2, 1e-12, 3  # cm to m, pC to C
            with open(fname) as f:
                f.readline()
                a = f.readline()
            mstep, offset, wid, bunlen = _np.fromstring(a[1:], sep='\t')
            offset = int(offset)
            # I don't know why I have to divide the echozr data by 2;
            arbitrary_factor = 2
            y0 = y = mstep*offset / 100
        elif code == 'echo2d':
            len_unit, charge_unit, header = 1, 1, 6
            with open(fname) as f:
                f.readline()
                mstep, offset = _np.fromstring(f.readline(), sep='\t')
                f.readline()
                wid, bunlen = _np.fromstring(f.readline(), sep='\t')
            offset = int(offset)
            # But I don't have to do this for the echo2d data.
            arbitrary_factor = 1
            y0 = y = mstep*offset
            offset = 0  # has only one column of wake
        spos, Wm = _np.loadtxt(
            fname, skiprows=header, usecols=(0, 1+offset), unpack=True)
        mstep *= len_unit
        wid *= len_unit
        bunlen *= len_unit
        spos *= len_unit
        # minus sign is due to convention:
        Wm *= -len_unit/charge_unit/arbitrary_factor

        Kxm = _np.pi/wid*mode
        if bc == 'elec':
            Wm /= _np.sinh(Kxm*y0)*_np.sinh(Kxm*y)
        else:
            Wm /= _np.cosh(Kxm*y0)*_np.cosh(Kxm*y)
        return spos, Wm, mstep, wid, bunlen

    if anal_pl == 'db':
        msg = 'Problem: All rectangular geometries does not have symmetry.'
        _log.info(msg)
        raise Exception(msg)

    bc = 'magn' if anal_pl == 'll' else 'elec'

    _log.info(f'Looking for data files in subfolder {bc:s}.')
    pname = _jnpth([path, bc])
    if not _os.path.isdir(pname):
        pname = path
        if code == 'echozr':
            _log.info(
                'Subfolder not found. It would be better to ' +
                'create the subfolder and put the files there...')
            _log.info('Looking for files in the current folder:')
        elif code == 'echo2d':
            msg = 'Files not found.'
            _log.info(msg)
            raise Exception(msg)

    f_in_dir = ' '.join(_os.listdir(pname))
    f_match = sorted(_re.findall(FNAME_ECHOZR2D, f_in_dir))
    if not f_match:
        msg = 'Files not found.'
        _log.info(msg)
        raise Exception(msg)

    _log.info(
        'Files found.\n I am assuming the simulation was performed ' +
        'with {0:s} boundary condition.'.format(
            'electric' if bc == 'elec' else 'magnetic'))
    _log.info('Modes found: ' + ', '.join([m for _, m in f_match]))
    _log.info('Loading data from files')

    spos, W, mode, mesh_size, width, bunlen = [], [], [], [], [], []
    for fn, m in f_match:
        if int(m) == 0:
            continue
        s, Wm, ms, wid, bl = _load_dados(_jnpth([pname, fn]), int(m), bc, code)
        mode.append(int(m))
        spos.append(s)
        W.append(Wm)
        mesh_size.append(ms)
        width.append(wid)
        bunlen.append(bl)

    cond = False
    for i in range(1, len(mode)):
        cond |= len(spos[i]) != len(spos[0])
        cond |= not _np.isclose(mesh_size[i], mesh_size[0], rtol=1e-5, atol=0)
        cond |= not _np.isclose(width[i], width[0], rtol=0, atol=1e-7)
        cond |= not _np.isclose(bunlen[i], bunlen[0], rtol=1e-5, atol=0)
        if cond:
            message = 'Parameters of file {0:s} differ from {1:s}.'.format(
                f_match[i][0], f_match[0][0])
            _log.info(message)
            raise Exception(message)

    simul_data.s = spos[0]
    simul_data.bunlen = bunlen[0]
    # I want the bunch to be symmetric:
    a = _np.argmin(_np.abs(spos[0] + spos[0][0])) + 1
    sbun = spos[0][:a]
    simul_data.sbun = sbun
    simul_data.bun = _np.exp(-sbun**2/(2*bunlen[0]**2))
    simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen[0]
    _log.info(f'Length of the driving bunch: {simul_data.bunlen*1e3:7.4g} mm')
    _log.info(f'Width of the simulated geometry: {width[0]*1e3:7.4g} mm')
    _log.info(f'Mesh step used in the simulation: {mesh_size[0]*1e6:7.4g} um')
    _log.info('All Data Loaded.')

    if anal_pl == 'll':
        _log.info('Calculating longitudinal Wake from data:')
        Wll, frac = None, 1
        for i in range(len(mode)):
            if mode[i] == 1:
                Wll = W[i].copy()
            elif mode[i] % 2:
                Wll += W[i]  # only odd terms
                frac = _np.max(_np.abs(W[i]/Wll))
        if Wll is None:
            _log.info('There is none odd mode to calculate Longitudinal Wake.')
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wll *= 2/width[0]
            simul_data.Wll = Wll

        _log.info('Calculating Quadrupolar Wake from data:')
        Wq, frac = None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wq = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wq += W[i] * Kxm**2  # only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wq))
        if Wq is None:
            _log.info('There is none odd mode to calculate Quadrupolar Wake.')
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wq *= 2/width[0]
            # minus sign is due to Panofsky-Wenzel
            Wq = -_scy_int.cumtrapz(Wq, x=spos[0], initial=0)
            simul_data.Wqy = Wq
            simul_data.Wqx = -Wq

        _log.info('Calculating Dipolar Horizontal Wake from data:')
        Wdx, frac = None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 2:
                Wdx = W[i].copy() * Kxm**2
            elif not mode[i] % 2:
                Wdx += W[i] * Kxm**2  # only even terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wdx))
        if Wdx is None:
            _log.info(
                "There's no even mode to calculate Dip. Horizontal Wake.")
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wdx *= 2/width[0]
            # minus sign is due to Panofsky-Wenzel
            Wdx = -_scy_int.cumtrapz(Wdx, x=spos[0], initial=0)
            simul_data.Wdx = Wdx

    elif anal_pl in {'dx', 'dy'}:
        pl = 'Vertical' if anal_pl == 'dy' else 'Horizontal'
        _log.info(f'Calculating Dipolar {pl:s} Wake from data:')
        Wd, frac = None, 1
        for i in range(len(mode)):
            Kxm = _np.pi/width[0]*mode[i]
            if mode[i] == 1:
                Wd = W[i].copy() * Kxm**2
            elif mode[i] % 2:
                Wd += W[i] * Kxm**2  # only odd terms
                frac = _np.max(_np.abs(W[i]*Kxm**2/Wd))
        if Wd is None:
            _log.info(
                f"There's no even mode to calculate Dipolar {pl:s} Wake.")
        else:
            _log.info(
                'Maximum influence of last mode in the final result '
                f'is: {frac*100:5.2f}%')
            Wd *= 2/width[0]
            # minus sign is due to Panofsky-Wenzel
            Wd = -_scy_int.cumtrapz(Wd, x=spos[0], initial=0)
            setattr(simul_data, 'W'+anal_pl, Wd)
    else:
        msg = f'Plane of analysis {anal_pl:s} does not match any '
        msg += 'of the possible options'
        _log.info(msg)
        raise Exception(msg)
