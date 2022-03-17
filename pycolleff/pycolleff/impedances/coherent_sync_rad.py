"""."""

import numpy as _np
import scipy.signal as _scysig
from scipy.special import hyp1f1 as _hyp1f1, gamma as _gamma, airy as _airy, \
    airye as _airye

import mathphys as _mp

_LSPEED = _mp.constants.light_speed
_Z0 = _mp.constants.vacuum_impedance


class CSRElement:
    _X = _np.linspace(-900, 900, 10001)
    _Y = None

    def __init__(
            self, rho=17.2, h=12e-3, bl=2.5e-3, nus=4.6e-3,
            espread=8.5e-4, rev_time=1.73e-6):
        self.bl = bl
        self.rho = rho
        self.h = h
        self.nus = nus
        self.espread = espread
        self.rev_time = rev_time

    @property
    def shielding(self):
        return self.rho**(1/2) * self.bl / self.h**(3/2)

    @property
    def threshold(self):
        S = self.calc_normalized_strength(1e-3) / 1e-3
        Ith = (0.5 + 0.12*self.shielding)/S
        return Ith

    def calc_normalized_current(self, I0):
        return (120/4*_LSPEED*I0*self.rev_time /
                (2*_np.pi*self.nus*3e9*self.espread))

    def calc_normalized_strength(self, I0):
        curr = self.calc_normalized_current(I0)
        return curr*self.rho**(1/3)/self.bl**(4/3)

    def calc_formation_length(self, w):
        return (24*self.rho**2/w/_LSPEED)**(1/3)

    def wake(self, z, maxi=25, L=None, bunlen=80e-6, convolved=True):
        gaus_to_si = 120*_LSPEED/4  # _Z0 _LSPEED / (4 pi)
        L = L or (2 * _np.pi * self.rho)

        # Free Space Term
        inds = z < 0
        W0 = _np.zeros(len(z))
        W0[inds] = (-2/3**(4/3)/self.rho**(2/3) /
                    _np.power((-z[inds]), 4/3)*gaus_to_si*L)

        # Shielding Term
        W1 = _np.zeros(len(z))
        zshield = self.shielding / self.bl * z
        for i in range(1, maxi):
            uai = self._getY(3*zshield/i**(3/2))
            W1 += 8*_np.pi*(-1)**(i+1)/i/i*uai*(3 - uai)/(1 + uai)**3
        W1 *= -1/self.h**2 / (2 * _np.pi) * gaus_to_si * L

        # If want the convolved values
        if convolved:
            # For the free space used analytical formulas of convolution
            bl = bunlen
            C = _Z0*_LSPEED*L
            C /= 2**(13/6) * _np.pi**(3/2) * (3*self.rho**2*bl**10)**(1/3)
            W0 = C*(2**(1/2)*_gamma(5/6)*(
                        bl**2*_hyp1f1(-1/3, 1/2, -z*z/2/bl**2) -
                        z**2*_hyp1f1(2/3, 3/2, -z*z/2/bl**2)) +
                    z*bl*_gamma(4/3)*(
                        3*_hyp1f1(1/6, 1/2, -z*z/2/bl**2) -
                        2*_hyp1f1(1/6, 3/2, -z*z/2/bl**2)))
            # For shielding perform convolution numerically
            bunch = _np.exp(-(z*z/bl**2)/2)/_np.sqrt(2*_np.pi)/bl  # gaussian
            W1 = _scysig.fftconvolve(W1, bunch, mode='same') * (z[1]-z[0])
        return W0, W1

    def impedance(self, w, L=None, imax=3, free=False):
        L = L or (2*_np.pi * self.rho)
        Lfrac = L / (2*_np.pi * self.rho)
        n = w/_LSPEED * self.rho

        if free:
            Z = (120/2 * 1.354/3**(1/3)*_np.exp(1j*_np.pi/6) *
                 (w/_LSPEED/self.rho**2)**(1/3) * L)
            return Z

        u0 = _np.pi**2/2**(2/3) * (n * (2*self.h/self.rho)**(3/2))**(-4/3)
        Z = (1 + 0*1j)*_Z0 * 16 * u0

        F = _np.zeros(len(w), dtype=complex)
        for p in range(0, imax):
            up = u0*(2*p + 1)**2
            Ai, Ail, Bi, Bil = _airy(up)
            Ri = Ail*Ail + up * Ai*Ai
            Ai, Ail, Bi, Bil = _airye(up)
            Im = Ail*Bil + up * Ai*Bi
            F += Ri - 1j * Im
        Z *= F
        return Z * n * self.h / self.rho * Lfrac

    def _getY(self, values):
        if self._Y is None:
            self._Y = _np.zeros(self._X.shape)
            for i, x in enumerate(self._X):
                sol = _np.roots([1, 0, 0, -x, -3])
                sol2 = []
                for s in sol:
                    if _np.isclose(s.imag, 0):
                        sol2.append(s.real)
                if len(sol2) != 2:
                    print('war')
                self._Y[i] = min(sol2)
            self._Y = self._Y*self._Y*self._Y*self._Y
        return _np.interp(values, self._X, self._Y, left=None, right=None)
