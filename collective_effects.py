#!/usr/bin/env python3

import numpy as np

class Element:
    def __init__(self,name='',betax=7.2,betay=11,quantity=1):
        self.name     = name
        self.betax    = betax
        self.betay    = betay
        self.w        = np.array([],dtype=float)
        self.Zl       = np.array([],dtype=float)
        self.Zdv      = np.array([],dtype=float)
        self.Zdh      = np.array([],dtype=float)
        self.Zqv      = np.array([],dtype=float)
        self.Zqh      = np.array([],dtype=float)
        self.Wqh      = np.array([],dtype=float)
        self.Wqh      = np.array([],dtype=float)
        self.Wqh      = np.array([],dtype=float)
        self.Wqh      = np.array([],dtype=float)
        self.Wqh      = np.array([],dtype=float)

class Budget(list):
    def __init__(self):
        super.__init__(self)


class Ring:
    def __init__(self,phase=None):
        self.version   = 'SI.v10.c01'
        self.circ      = 518.396
        self.T0        = self.circ/299792458
        self.f0        = 1/self.T0
        self.w0        = 2*np.pi*self.f0 # revolution angular frequency [Hz]
        self.mom_cmpct = 1.7e-4       # momentum compaction factor
        self.E         = 3              # energy [GeV]
        self.nuy       = 13.116
        self.nux       = 48.131
        self.harm_num  = 864           # harmonic Number

        I = np.linspace(0,4,num=40)
        self.I   = I*1e-3
        if phase.startswith('commissioning'):
            self.phase   = 'Commissioning'
            self.nom_cur = 0.100       # total current [A]
            self.nus     = 0.00435    # synchrotron tune
            self.espread = 7.64e-4 + 0*I
            self.emitx   = 271e-12 + 0*I
            self.emity   = 2.71e-12 + 0*I
            # damping times [s]
            self.damptx = 17.1e-3
            self.dampty = 22.7e-3
            self.dampte = 13.6e-3
            self.en_lost_rad = 456740.6 #eV
        elif phase.startswith('phase_1'):
            self.phase   = 'Phase 1 with IDS'
            self.I_tot = 0.10         # total current [A]
            self.nus = 0.00435        # synchrotron tune
            self.espread = 1e-02*(0.093995 + 0.038011*I -0.018279*I.^2 + 0.0047843*I.^3 -0.00047294*I.^4)
            self.emitx   = 1e-09*(0.23011  + 0.15699 *I -0.063581*I.^2 + 0.015965* I.^3 -0.0015505* I.^4)
            self.emity   = 1e-12*(2.1496   + 1.8725*  I -0.84932 *I.^2 + 0.22507*  I.^3 -0.022538 * I.^4)
            # damping times [s]
            self.taux = 12.4e-3
            self.tauy = 15.1e-3
            self.taue =  8.5e-3
            self.U0   = 685374.1 #eV
        elif phase.startswith('phase_2'):
            self.stage   = 'Phase 2 with IDS'
            self.I_tot   = 0.35        # total current [A]
            self.nus     = 0.00435    # synchrotron tune
            self.espread = 1e-02*(0.088704 + 0.015765*I -0.005477*I.^2 + 0.0012452*I.^3 -0.00011434*I.^4)
            self.emitx   = 1e-09*(0.18859  + 0.056781*I -0.015909*I.^2 + 0.003445* I.^3 -0.00031039*I.^4)
            self.emity   = 1e-12*(1.6497   + 1.0422*  I -0.51454 *I.^2 + 0.14498*  I.^3 -0.015059 * I.^4)
            # damping times [s]
            self.taux = 10.6e-3
            self.tauy = 12.5e-3
            self.taue =  6.9e-3
            self.U0   = 829761.9 #eV
        elif phase.startswith('phase_2_HC'):
            self.stage   = 'Phase 2 with IDS High Current'
            self.I_tot   = 0.5        # total current [A]
            self.nus     = 0.00435    # synchrotron tune
            self.espread = 1e-02*(0.088704 + 0.015765*I -0.005477*I.^2 + 0.0012452*I.^3 -0.00011434*I.^4)
            self.emitx   = 1e-09*(0.18859  + 0.056781*I -0.015909*I.^2 + 0.003445* I.^3 -0.00031039*I.^4)
            self.emity   = 1e-12*(1.6497   + 1.0422*  I -0.51454 *I.^2 + 0.14498*  I.^3 -0.015059 * I.^4)
            # damping times [s]
            self.taux = 10.6e-3
            self.tauy = 12.5e-3
            self.taue =  6.9e-3
            self.U0   = 829761.9 #eV
