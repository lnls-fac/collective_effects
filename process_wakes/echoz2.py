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
