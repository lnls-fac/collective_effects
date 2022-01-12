"""Load and process ECHO3D simulation results.

The coordinate system used in ECHO3D is not the same as sirius:
      ECHO3D      --->       Sirius
    (x, y, z)     --->     (z, x, y)

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
from os.path import join as _jnpth, isdir as _isdir
import re as _re
import logging as _log
import json as _json

import numpy as _np
from scipy import integrate as _scy_int

from .base_load import LoadWakes


class ECHO3D(LoadWakes):
    """."""

    FNAME = r'wake3Dindirect.bin'

    def load_and_process_raw_data(self, path: str, anal_pl: str):
        """Load raw simulaton data and process it.

        Args:
            path (str): path to the folder where data is placed.
            anal_pl (str): plane where the analysis was performed.

        """
        anal_pl, anal_pl_comp, elec_symm = self._get_analysis_plane(
            path, anal_pl)

        if anal_pl == 'll':
            self._load_longitudinal_data(path)
            return

        if elec_symm:
            path = _jnpth(path, 'elec')
            self._load_transverse_symetric_data(
                path, anal_pl, elec_symm=elec_symm)
        else:
            self._load_transverse_assymetric_data(path, anal_pl)
        _log.info('Transverse Data Loaded.')

        if anal_pl_comp:
            _log.info(
                'There is symmetry. Copying the data from the ' +
                f'{anal_pl[1].upper():s} plane to the ' +
                f'{anal_pl_comp[1].upper():s} plane')
            setattr(
                self.simul_data, 'W'+anal_pl_comp,
                getattr(self.simul_data, 'W'+anal_pl).copy())

    def _get_analysis_plane(self, path, anal_pl):
        # list all files that match the name pattern for wakefields
        f_match = self._find_files(path)

        if anal_pl == 'll':
            if not f_match:
                msg = 'No files found for longitudinal analysis.'
                _log.info(msg)
                raise Exception(msg)
            return anal_pl, None, False

        elec_symm = _isdir(_jnpth(path, 'elec'))
        if anal_pl == 'db':
            if f_match:
                msg = (
                    "I know this simulation works for both planes, but I " +
                    "couldn't find out which plane was in fact used.")
                _log.info(msg)
                raise Exception(msg)

            _log.info(
                'There is symmetry y=x, calculation performed in the ' +
                anal_pl[1].upper() + ' plane.')

            anal_pl = 'dx' if _isdir('dxdpl') else 'dy'
            anal_pl_comp = ({'dx', 'dy'} - {anal_pl}).pop()
            return anal_pl, anal_pl_comp, elec_symm

        if anal_pl not in {'dx', 'dy'}:
            msg = f'Plane of analysis {anal_pl:s} does not match any of '
            msg += 'the possible options'
            _log.info(msg)
            raise Exception(msg)

        if not f_match:
            _log.info('There is no wake files in this folder.')
            if elec_symm:
                _log.info(
                    ' I found a folder named "elec". I will assume the '
                    'simulation has this symmetry.')
            else:
                _log.info(' I will assume there is no symmetry.')

        return anal_pl, None, elec_symm

    def _load_transverse_symetric_data(self, path, anal_pl, elec_symm=False):
        simul_data = self.simul_data

        wake, origx, origy = self._load_wakes(path)
        nz = wake.shape[0]

        dic = self._load_info(path)
        bunlen = dic['BunchSigma']*1e-3  # in m
        mstepz = dic['Steps'][0]*1e-3  # in m
        mstepx = dic['Steps'][1]*1e-3  # in m
        mstepy = dic['Steps'][2]*1e-3  # in m
        bunx = dic['BunchPosition'][0]  # in step index
        buny = dic['BunchPosition'][1]  # in step index

        spos = mstepz*_np.arange(nz, dtype=float)
        spos -= 5 * bunlen
        a = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        bun = _np.exp(-sbun**2/(2*bunlen*bunlen))
        bun /= _np.sqrt(2*_np.pi)*bunlen

        simul_data.s = spos
        simul_data.bun = bun
        simul_data.sbun = sbun
        simul_data.bunlen = bunlen

        # Panofsky-Wenzel:
        # dWt/dz = -dWl(t, z)/dt
        # Wt = -int(dWl(t, zp)/dt, zp=-inf, z)
        # with t in {x, y}
        orig = origx if anal_pl == 'dx' else origy
        origo = origy if anal_pl == 'dx' else origx
        mstep = mstepx if anal_pl == 'dx' else mstepy
        bunp = bunx if anal_pl == 'dx' else buny

        axiso = 2 if anal_pl == 'dx' else 1

        waket = wake.take(origo, axis=axiso)
        waket = _np.gradient(waket, mstep, axis=1)[:, orig]
        waket = -_scy_int.cumtrapz(waket, x=spos, initial=0, axis=0)

        mstep *= 2 if elec_symm else 1
        waket /= mstep * (bunp-orig)
        setattr(simul_data, 'W'+anal_pl, waket)  # V/C/m

    def _load_transverse_assymetric_data(self, path, anal_pl):
        simul_data = self.simul_data
        sposs, wakes, sbuns, buns, bunlens, xd, yd = [], [], [], [], [], [], []
        origs, bunsp, msteps = [], [], []

        for sub_fol in ['dpl', 'dmi']:
            ext_path = _jnpth(path, anal_pl+sub_fol)
            _log.info('Looking for '+anal_pl+sub_fol+' subfolder:')
            if not _isdir(ext_path):
                _log.info(
                    'For non-symmetric structures, there must '
                    f'be subfolders {anal_pl:s}dpl {anal_pl:s}dmi '
                    'with the data')
                raise Exception('Files not found')
            # list all files that match the pattern
            f_match = self._find_files(ext_path)
            if not f_match:
                msg = 'No files found for transverse analysis.'
                _log.info(msg)
                raise Exception(msg)

            _log.info('Loading wake file and calulating {0:s} wake:'.format(
                'Horizontal' if anal_pl == 'dx' else 'Vertical'))
            wake, origx, origy = self._load_wakes(ext_path)
            nz = wake.shape[0]

            dic = self._load_info(ext_path)
            bunlen = dic['BunchSigma']*1e-3  # in m
            mstepz = dic['Steps'][0]*1e-3  # in m
            mstepx = dic['Steps'][1]*1e-3  # in m
            mstepy = dic['Steps'][2]*1e-3  # in m
            bunx = dic['BunchPosition'][0]  # in step index
            buny = dic['BunchPosition'][1]  # in step index

            spos = mstepz*_np.arange(nz, dtype=float)
            spos -= 5 * bunlen
            a = _np.argmin(_np.abs(spos + spos[0])) + 1
            sbun = spos[:a]
            bun = _np.exp(-sbun**2/(2*bunlen*bunlen))
            bun /= _np.sqrt(2*_np.pi)*bunlen

            sposs.append(spos)
            buns.append(bun)
            sbuns.append(sbun)
            bunlens.append(bunlen)
            xd.append(bunx)
            yd.append(buny)

            # Panofsky-Wenzel:
            # dWt/dz = -dWl(t, z)/dt
            # Wt = -int(dWl(t, zp)/dt, zp=-inf, z)
            # with t in {x, y}
            orig = origx if anal_pl == 'dx' else origy
            origo = origy if anal_pl == 'dx' else origx
            mstep = mstepx if anal_pl == 'dx' else mstepy
            bunp = bunx if anal_pl == 'dx' else buny
            origs.append(orig)
            bunsp.append(bunp)
            msteps.append(mstep)

            waket = wake.take(origo, axis=2 if anal_pl == 'dx' else 1)
            waket = _np.gradient(waket, mstep, axis=1)
            waket = -_scy_int.cumtrapz(waket, x=spos, initial=0, axis=0)
            wakes.append(waket[:, orig])

        # If the simulation is not ready yet the lenghts may differ.
        # This line is used to truncate the longer wake in the length of
        # the shorter:
        l1 = min(len(sposs[0]), len(sposs[1]))

        bunp = xd if anal_pl == 'dx' else yd
        ndel = yd if anal_pl == 'dx' else xd
        if not (_np.allclose(sposs[0][:l1], sposs[1][:l1], atol=0) and
                _np.allclose(sbuns[0], sbuns[1], atol=0) and
                _np.allclose(buns[0], buns[1], atol=0) and
                _np.allclose(ndel[0], ndel[1], atol=0)):
            msg = 'There is a mismatch between the parameters of the'
            msg += f'simulation in the {anal_pl:s}dpl and {anal_pl:s}dmi '
            msg += 'folders.'
            _log.info(msg)
            raise Exception(msg)
        simul_data.s = sposs[0][:l1]
        simul_data.bun = buns[0]
        simul_data.sbun = sbuns[0]
        simul_data.bunlen = bunlens[0]

        waket = wakes[0][:l1] - wakes[1][:l1]
        waket /= (bunsp[0]-origs[0])*msteps[0] - (bunsp[1]-origs[1])*msteps[1]
        setattr(simul_data, 'W'+anal_pl, waket)  # V / pC / m

    def _load_longitudinal_data(self, path):
        simul_data = self.simul_data

        wake, origx, origy = self._load_wakes(path)
        nz = wake.shape[0]

        dic = self._load_info(path)
        bunlen = dic['BunchSigma']*1e-3  # in m
        mstepz = dic['Steps'][0]*1e-3  # in m
        mstepx = dic['Steps'][1]*1e-3  # in m
        mstepy = dic['Steps'][2]*1e-3  # in m
        bunx = dic['BunchPosition'][0]  # in step index
        buny = dic['BunchPosition'][1]  # in step index

        if bunx != origx:
            _log.warning(
                'Attention, Bunch was not passed at origin X for '
                'longitudinal wake simulation.'
                'Results might not be accurate.')
        if buny != origy:
            _log.warning(
                'Attention, Bunch was not passed at origin Y for '
                'longitudinal wake simulation.'
                'Results might not be accurate.')

        spos = mstepz*_np.arange(nz, dtype=float)
        spos -= 5 * bunlen
        a = _np.argmin(_np.abs(spos + spos[0])) + 1
        sbun = spos[:a]
        bun = _np.exp(-sbun**2/(2*bunlen*bunlen))
        bun /= _np.sqrt(2*_np.pi)*bunlen

        simul_data.s = spos
        simul_data.bunlen = bunlen
        simul_data.sbun = sbun
        simul_data.bun = bun
        _log.info(
            'Bunch length of the driving bunch: ' +
            f'{simul_data.bunlen*1e3:7.3g} mm')
        _log.info(
            f'Mesh size used in the simulation: {mstepz*1e6:7.4g} um')

        # Load longitudinal Wake
        simul_data.Wll = wake[:, origx, origy]

        # Horizontal quadrupolar wake
        _log.info('Loading Horizontal Quadrupolar Wake file:')
        # Panofsky-Wenzel:
        # dWx/dz = -dWl(x, y, z)/dx
        # Wx = -int(dWl(x, y, zp)/dx, zp=-inf, z)
        wakex = _np.gradient(wake[:, :, origy], mstepx, axis=1)
        wakex = -_scy_int.cumtrapz(wakex, x=spos, initial=0, axis=0)
        # Isolate the quadrupolar wake:
        # Wqx = dWx(z)/dx
        simul_data.Wqx = _np.gradient(wakex, mstepx, axis=1)[:, origx]  # V/C/m

        # Vertical quadrupolar wake
        _log.info('Loading Vertical Quadrupolar Wake file:')
        # Panofsky-Wenzel:
        # dWy/dz = -dWl(x, y, z)/dy
        # Wy = -int(dWl(x, y, zp)/dy, zp=-inf, z)
        wakey = _np.gradient(wake[:, origx, :], mstepy, axis=1)
        wakey = -_scy_int.cumtrapz(wakey, x=spos, initial=0, axis=0)
        # Isolate the quadrupolar wake:
        # Wqy = dWy(z)/dy
        simul_data.Wqy = _np.gradient(wakey, mstepy, axis=1)[:, origy]  # V/C/m
        _log.info('Longitudinal Data Loaded.')

    @classmethod
    def _load_wakes(cls, path):
        fname = _jnpth(path, 'Results', cls.FNAME)
        nz, nx, ny = _np.fromfile(fname, dtype=_np.int32, count=3)
        data = _np.fromfile(fname, dtype=_np.float64, offset=4*3)
        x = data[:nx] * 1e-2  # from cm to m
        y = data[nx:nx+ny] * 1e-2  # from cm to m
        wake = data[nx+ny:]
        wake = wake.reshape(nz, nx, ny)
        wake *= -1  # the minus sign is due to convention
        wake *= 1e12  # from V/pC to V/C

        # find the index with the origin
        origx = _np.isclose(x, 0.0).nonzero()[0]
        if not origx.size:
            _log.warning(
                'Attention, wake was not calculated at origin X. '
                'Results might not be accurate.')
            origx = _np.argmin(_np.abs(x))
        else:
            origx = origx[0]

        origy = _np.isclose(y, 0.0).nonzero()[0]
        if not origy.size:
            _log.warning(
                'Attention, wake was not calculated at origin Y. '
                'Results might not be accurate.')
            origy = _np.argmin(_np.abs(y))
        else:
            origy = origy[0]
        return wake, origx, origy

    @staticmethod
    def _load_info(path):
        with open(_jnpth(path, 'input.txt'), 'r') as fil:
            data = fil.readlines()
        dic = dict()
        for line in data:
            if '=' not in line:
                continue
            lin = line.split('%')[0]
            var, val = lin.split('=')
            var = var.replace(' ', '')
            val = val.strip()
            if "'" in val:
                dic[var] = val.replace("'", '')
            else:
                valo = ''
                while val != valo:
                    valo = val
                    val = valo.replace('  ', ' ')
                val = val.replace(' ', ',')
                dic[var] = _json.loads(val)
        return dic

    @classmethod
    def _find_files(cls, path):
        path_ = _jnpth(path, 'Results')
        f_match = None
        if _isdir(path_):
            f_in_dir = ' '.join(_listdir(path_))
            f_match = _re.findall(cls.FNAME, f_in_dir)
        return f_match
