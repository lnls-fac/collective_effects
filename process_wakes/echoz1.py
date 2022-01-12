"""Load and process ECHOZ1 simulation results.

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

from .base_load import LoadWakes


class ECHOZ1(LoadWakes):
    """."""

    FNAME = r"wake.dat"

    def load_and_process_raw_data(self, path: str, anal_pl: str):
        """Load raw simulaton data and process it.

        Args:
            path (str): path to the folder where data is placed.
            anal_pl (str): plane where the analysis was performed.

        """
        if anal_pl == 'll':
            _log.info('Loading longitudinal Wake file:')
            fname = _jnpth([path, self.FNAME])
            if _isfile(fname):
                _log.info('Data found.')
                loadres = _np.loadtxt(fname, skiprows=0)
            else:
                _log.info('Not found.')
                raise Exception('Longitudinal wake file not found.')
        else:
            msg = 'ECHOz1 only calculates longitudinal wake.'
            _log.info(msg)
            raise Exception(msg)

        self.simul_data.s = loadres[:, 0]/100  # Rescaling cm to m
        # V/C/m (the minus sign is due to convention):
        self.simul_data.Wll = -loadres[:, 1] * 1e12

        # loading driving bunch info
        _log.info('Loading bunch length from wake.dat')
        sbun = self.simul_data.s.copy()
        ds = sbun[1] - sbun[0]
        bunlen = abs(sbun[0]-ds/2) / 5
        a = _np.argmin(_np.abs(sbun + sbun[0])) + 1
        sbun = sbun[:a]
        self.simul_data.bunlen = bunlen
        self.simul_data.sbun = sbun
        self.simul_data.bun = _np.exp(-sbun**2/(2*bunlen**2))
        self.simul_data.bun /= _np.sqrt(2*_np.pi)*bunlen
        _log.info(f'Bunch length of the driving bunch: {bunlen*1e3:7.3g} mm')
        _log.info('Data Loaded.')

    def _get_analysis_plane(self, path, anal_pl):
        # list all files that match the name pattern for wakefields
        f_match = self._find_files(path)
        return anal_pl

    @classmethod
    def _find_files(cls, path):
        f_in_dir = ' '.join(_listdir(path))
        f_match = _re.findall(cls.FNAME, f_in_dir)
        return f_match
