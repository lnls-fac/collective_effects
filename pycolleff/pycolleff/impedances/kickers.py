"""Calculate Coupled and Uncoupled Ferrite Kicker impedances."""

import numpy as _np

import mathphys as _mp

_c = _mp.constants.light_speed
_mu0 = _mp.constants.vacuum_permeability
_Z0 = _mp.constants.vacuum_impedance
E0 = _mp.constants.electron_rest_energy * _mp.units.joule_2_eV


def kicker_coupled_flux(w, h, W, t, L, mur, Zg):
    r"""Calculate Impedances for a ferrite kicker.

      - For the Coupled Flux, it uses Davino-Hahn model.

      DAVINO-HAHN MODEL:

     #######################################    |
     ###############FERRITE#################    t
     #######################################    |
     ###**                             **###  |
     ###**  VACUUM                     **###  |      ______|  |_________|
     ###**                             **###  |            |  |         |
     ###**             +   .           **###  w            |  |         |
     ###**                             **###  |            )||(         \
     ###**             |_D_|           **###  |      Zk  L1)||(L2     Zg/
     ###**                             **###  |            )||(         \
     #######################################               | M|         /
     #######################################               |  |         |
     #######################################         ______|  |_________|
         |______________h______________|

    Bibliografias:

    - Davino_D Hahn_H - Improved Analytical Model of the transverse coupling
      impedance of ferrite kicker magnets - Phys. Rev. ST-AB v6 012001 2003

    - Nassibian G Sacherer F - Methods for measuring tranverse coupling
      impedances in circular Accelerators,
      Nucl Inst and Meth. 159 21-27 1979

    """
    # Equivalent Circuit model.
    D = 0.5e-3
    M = L*D*_mu0/W
    L2 = L*h*_mu0/W
    # L2 = L*h*_mu0/W*(mur*t/(mur*t+h*(h/W+1)))

    Zk = (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg) * w
    Zx = (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg) * _c/D**2

    return Zk.conj(), Zx.conj()  # conjugate to adapt impedance convention


def kicker_tsutsui_model(w, epr, mur, a, b, d, L, n):
    """Calculate uncoupled flux impedance of a kicker magnet.

    - For the Uncoupled Flux, we can choose between three models:

      TSUTSUI MODEL:

      ******************PEC**********************
      *******************************************
      **#######################################**      |
      **################FERRITE################**      |
      **#######################################**      |
      **                                       **  |   d
      **                                       **  b   |
      **                                       **  |   |
      **     VACUUM        .                   **  |   |
      **                                       **
      **                                       **
      **                                       **
      **#######################################**
      **#######################################**
      **#######################################**
      *******************************************
      *******************************************
                           |__________a________|

    Inputs:

    w   = vector of angular frequencies to evaluate impedances [rad/s]
    epr = vector with real and imaginary electric permeability of ferrite for
          the frequency range of interest
    mur = vector with real and imaginary magnetic permissivity of ferrite for
          the frequency range of interest
    n   = max order of terms to sum
    L   = length of the structure [m]

    Outputs:

    Zl = Impedancia Longitudinal [Ohm]
    Zh = Impedancia Horizontal [Ohm/m]
    Zv = Impedancia Vertical   [Ohm/m]

    Bibliografias:

    - Tsutsui_H - Some Simplified Models of Ferrite Kicker Magnet for
      Calculation of longitudinal Coupling Impedance - CERN-SL-2000-004

    - Tsutsui_H - Transverse Coupling Impedance of a Simplified Ferrite
      Kicker Magnet Model - LHC Project Note 234 - 2000

    - Salvant, B. et al - Quadrupolar Impedance of Simple
    Models of Kickers, Proceedings of IPAC 2010, pp. 2054-2057

    """
    # Valores do artigo do Wang et al para testar a implementacao das formulas
    # do modelo do Tsutui.
    # a = 103e-3
    # b = 67e-3
    # d = 122e-3
    # L = 1.0
    # Valores do artigo do Tsutsui para testar a implementacao das formulas
    # a = 20e-3
    # b = 16e-3
    # d = 66e-3
    # L = 0.6

    # Terms for the infinite sum:
    n = _np.arange(0, n+1)[:, None]

    k = _np.ones(n.shape)*w/_c
    epr = _np.ones(n.shape)*epr
    mur = _np.ones(n.shape)*mur

    kxn = _np.repeat((2*n + 1)*_np.pi/2/a, w.shape[0], axis=1)
    kyn = _np.sqrt((epr*mur-1)*k**2 - kxn**2)
    sh = _np.sinh(kxn*b)
    ch = _np.cosh(kxn*b)
    tn = _np.tan(kyn*(b-d))
    ct = 1/_np.tan(kyn*(b-d))

    Zl = 1j*_Z0*L/2/a / (
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(
            mur*sh**2*tn - epr*ch**2*ct)
         )/(epr*mur-1) - k/kxn*sh*ch)
    Zq = -Zl*kxn*kxn/k
    Zq = Zq.sum(0)
    Zl = Zl.sum(0)

    Zv = 1j*_Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(
            mur*ch**2*tn - epr*sh**2*ct)
         )/(epr*mur-1) - k/kxn*sh*ch)
    Zv = Zv.sum(0)

    kxn = _np.repeat(2*(n + 1)*_np.pi/2/a, w.shape[0], axis=1)
    kyn = _np.sqrt((epr*mur-1)*k**2 - kxn**2)
    sh = _np.sinh(kxn*b)
    ch = _np.cosh(kxn*b)
    tn = _np.tan(kyn*(b-d))
    ct = 1/_np.tan(kyn*(b-d))

    Zh = 1j*_Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(
            mur*sh**2*tn - epr*ch**2*ct)
         )/(epr*mur-1) - k/kxn*sh*ch)
    Zh = Zh.sum(0)

    return Zl.conj(), Zh.conj(), Zv.conj(), Zq.conj()  # impedance convention
