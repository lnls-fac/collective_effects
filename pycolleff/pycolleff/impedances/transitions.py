"""Implement some impedance models for transitions."""

import numpy as _np
import scipy.constants as _constants


_LSPEED = _constants.speed_of_light
_Z0 = _constants.physical_constants['characteristic impedance of vacuum'][0]


def taper(w, r1, r2, t, wid=0, geom='round'):
    r"""Calculate symmetric taper impedance.

    Geometry:
                        L2
        _____|- - - - - - - - - - - - |_____
              \                      /    :
               \                    /     :
                \                  /      :
                 \_______L1_______/       : R2
                    :                     :
                    : R1                  :
                    :                     :
        - - - - - - - - - - - - - - - - - - -

    """
    diff = _np.abs(r2-r1)
    summ = r2 + r1
    prod = r2*r1
    ums = _np.ones(w.shape)
    if geom == 'round':
        Zll = 2 * -1j*w*_Z0/4/_np.pi/_LSPEED * diff/t
        Zdx = 2 * -1j*_Z0/2/_np.pi * diff/t/prod * ums
        Zdy = 1*Zdx
        Zqx = 0*w
    else:
        Zll = 2 * -1j*0.43*w*_Z0/_np.pi/_LSPEED * diff/t
        Zdx = 2 * -1j*_Z0/4/_np.pi * diff/t/prod * ums
        Zdy = 2 * -1j*_Z0/2 * wid/t*summ*diff/prod/prod * ums
        Zqx = -1*Zdx
    return Zll, Zdx, Zdy, Zqx
