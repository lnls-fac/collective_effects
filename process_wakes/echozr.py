"""Load and process ECHOZ2 simulation results.

In case the plane of analysis is the transverse plane ('dx' or 'dy')
it is preferable to get the transverse wake from the longitudinal one
because of an integration error in ECHOZ2 when it calculates itself
the transverse wake.

Troughout the code I am assuming:
    s positive means particle behind source -->  Wl, Wt = 0 s < 0


Definitions:
    - Longitudinal Wake:
        Wl(s) = -c/Q * int El(ct-s,t) dt

    - Dipolar wake:
        Wdx(s) = - 1/xs * int_-inf^s dWl(s)/dx ds'
    where xs is source particle displacement

    - Quadrupolar wake:
        Wqx(s) = - int_-inf^s dWl(s)/dx^2 ds'

    - Longitudinal impedance:
        Zl = int exp(i*w*s) Wl(s) ds

    - Dipolar and Quadrupolar impedances:
        Zx = i*int exp(i*w*s) Wx(s) ds

"""

from os import listdir as _listdir
from os.path import join as _jnpth, isfile as _isfile
import re as _re
import logging as _log

import numpy as _np
from scipy import integrate as _scy_int

from .base_load import LoadWakes


class ECHOZ2(LoadWakes):
    """."""

    FNAME = r"wake[LT]{1}.dat"

    def load_and_process_raw_data(self, path: str, anal_pl: str):
        """Load raw simulaton data and process it.

        Args:
            path (str): path to the folder where data is placed.
            anal_pl (str): plane where the analysis was performed.

        """
        anal_pl_ori = None
        if anal_pl == 'db':
            anal_pl_ori = 'db'
            anal_pl = 'dy'
            _log.info(
                'Even though there is symmetry, '
                'I am loading data to the Y plane.')

        if anal_pl == 'll':
            _log.info('Loading longitudinal Wake file:')
            fname = _jnpth([path, 'wakeL.dat'])
            if _isfile(fname):
                _log.info('Data found.')
                spos, wl = _np.loadtxt(
                    fname, skiprows=0, usecols=(0, 1), unpack=True)
            else:
                _log.info('Not found.')
                Exception('Longitudinal wake file not found.')

            self.simul_data.s = spos/100  # Rescaling cm to m
            # V/C (minus sign is due to convention)
            self.simul_data.Wll = -wl * 1e12
        elif anal_pl in {'dx', 'dy'}:
            fname = _jnpth([path, 'wakeL.dat'])
            if _isfile(fname):
                _log.info(
                    'Calculating Transverse wake from longitudinal wake file:')
                _log.info('Data found.')
                spos, wl = _np.loadtxt(
                    fname, skiprows=0, usecols=(0, 1), unpack=True)
                self.simul_data.s = spos/100  # Rescaling cm to m
                # one minus sign due to convention and
                # the other due to Panofsky-Wenzel:
                wt = -_scy_int.cumtrapz(-wl, x=spos/100, initial=0)
                setattr(self.simul_data, 'W'+anal_pl, wt * 1e12)  # V/C/m
            else:
                _log.info('File not found.')
                _log.info(
                    'Loading transverse wake from transverse wake file.:')
                fname = _jnpth([path, 'wakeT.dat'])
                if _isfile(fname):
                    _log.info('Data found.')
                    _log.info(
                        'Depending on the ECHOz2 program version this '
                        'may lead to inacurate results.')
                    spos, wt = _np.loadtxt(
                        fname, skiprows=0, usecols=(0, 1), unpack=True)
                else:
                    _log.info('Not found.')
                    Exception('Transverse wake file not found.')
                self.simul_data.s = spos/100  # Rescaling cm to m
                # there is an error in the integration of echoz2.
                # It is needed to subtract the first value to correct
                # an offset
                # wt = -_scy_int.cumtrapz(-wl, x=spos/100, initial=0)
                setattr(self.simul_data, 'W'+anal_pl, (wt - wt[0]) * 1e12)
                # V/C/m (minus sign is due to convention)
        else:
            msg = f'Plane of analysis {anal_pl:s} does not match any '
            msg += 'of the possible options'
            _log.info(msg)
            raise Exception(msg)

        # loading driving bunch info
        _log.info('Loading bunch length from wake file')
        sbun = self.simul_data.s.copy()
        ds = sbun[1]-sbun[0]
        bunlen = abs(sbun[0] - ds/2) / 5
        a = _np.argmin(_np.abs(sbun + sbun[0])) + 1
        sbun = sbun[:a]
        self.simul_data.bunlen = bunlen
        self.simul_data.sbun = sbun
        self.simul_data.bun = _np.exp(-sbun**2/(2*bunlen**2))
        self.simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen
        _log.info(
            'Bunch length of the driving bunch: '
            f'{self.simul_data.bunlen*1e3:7.3g} mm')
        _log.info('Data Loaded.')

        if anal_pl_ori:
            anal_pl_comp = 'dx' if anal_pl == 'dy' else 'dy'
            _log.info(
                'There is symmetry. Copying the data from the ' +
                f'{anal_pl[1].upper():s} plane to ' +
                f'the {anal_pl_comp[1].upper():s} plane')
            setattr(
                self.simul_data, 'W'+anal_pl_comp,
                getattr(self.simul_data, 'W' + anal_pl).copy())

    def _get_analysis_plane(self, path, anal_pl):
        # list all files that match the name pattern for wakefields
        f_match = self._find_files(path)
        return anal_pl

    @classmethod
    def _find_files(cls, path):
        f_in_dir = ' '.join(_listdir(path))
        f_match = _re.findall(cls.FNAME, f_in_dir)
        return f_match


def _ECHOzR_load_data(simul_data, path, anal_pl):
    _ECHO_rect_load_data(simul_data, 'echozr', path, anal_pl)

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
