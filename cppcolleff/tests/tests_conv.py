#!/usr/bin/env python-sirius

import time
import math
import numpy as np
import cppcolleff as coll
import pycolleff.sirius as si


def perform_time_domain_simulation(rin, wake, currents,
                                   init_espread, name, niter=100):
    espread = np.zeros(len(currents))
    s0 = np.zeros(len(currents))
    bl = np.zeros(len(currents))
    nus_ave = np.zeros(len(currents))
    nus_std = np.zeros(len(currents))
    Vs, dists = [], []
    string = '{0:^17s} {1:^17s} {2:^17s} {3:^17s} {4:^17s} {5:^17s}'.format(
              'Ib [mA]', 'spread x 1000', 'bun len [mm]', 'sync ph [mm]',
              'ave(nus) x 1000', 'std(nus) x 1000\n')
    print(string, end='')
    with open(name+'_data.txt', 'w') as fi:
        fi.write(string)
    dist = rin.get_distribution()
    for cur, i in zip(currents, range(len(currents))):
        rin.espread = max(init_espread[i], rin.espread)
        espread[i] = coll.find_equilibrium_energy_spread(wake, rin, cur,
                                                         niter, dist)
        rin.espread = espread[i]
        V = coll.solve_Haissinski_get_potential(wake, rin, cur,
                                                niter*2, dist)
        dist = rin.get_distribution(V)
        dists.append(dist)
        Vs.append(V)

        ss = rin.cav.ref_to_xi()
        # calc average and standard deviation of synchrotron tune
        der = np.abs(np.diff(V))/(ss[1]-ss[0]) * rin.mom_comp * rin.circum
        nus_ave[i] = np.trapz(np.sqrt(der)*np.array(dist)[:-1],
                              x=ss[:-1]) / 2 / np.pi
        x = np.trapz(der*np.array(dist)[:-1],
                     x=ss[:-1]) / (2*np.pi)**2 - nus_ave[i]**2
        nus_std[i] = np.sqrt(x if x >= 0 else 0.0)

        # phase shift:
        s0[i] = np.trapz(np.array(dist)*np.array(ss), x=ss)

        # bunch length
        bl[i] = np.sqrt(np.trapz(np.array(dist)*np.array(ss)*np.array(ss),
                                 x=ss) - s0[i]*s0[i])

        string = \
         '{0:^17.3g} {1:^17.3g} {2:^17.3g} {3:^17.3g} {4:^17.3g} {5:^17.3g}\n'\
         .format(cur*1e3, espread[i]*1e3, bl[i]*1e3, s0[i]*1e3,
                 nus_ave[i]*1e3, nus_std[i]*1e3)
        print(string, end='')
        with open(name+'_data.txt', 'a') as fi:
            fi.write(string)
    rin.espread = init_espread[0]

    return currents, espread, bl, s0, nus_ave, nus_std, dists


coll.set_num_threads(8)

ring = coll.Ring_t()
ring.energy = 3e9
ring.en_lost_rad = 470e3/ring.energy
ring.tunex = 0.13
ring.chromx = 0
ring.emitx = 250e-12
ring.espread = 8e-4
ring.circum = 518.396
ring.T0 = 518.396/coll.light_speed
ring.mom_comp = 1.7e-4
ring.harm_num = 864
ring.betax = 19
ring.etax = 0
ring.damp_nume = 1.7
ring.damp_numx = 1.3

V0 = 3e6/ring.energy
phi0 = coll.TWOPI/2 - math.asin(ring.en_lost_rad/V0)
krf = coll.TWOPI*ring.harm_num/ring.circum

ss = coll.my_Dvector()
V = coll.my_Dvector()
for i in range(-10000, 10001):
    s = 1e-4 * 15e-2 * i
    ss.push_back(s)
    V.push_back(V0*math.sin(phi0 + krf*s) - ring.en_lost_rad)
ring.cav.set_xy(ss, V)

num_part = 100000
nturns = 10000
bun = coll.Bunch_t(num_part, 1e-3)
coll.generate_bunch(ring, bun)
bun.sort()


resonators = [  # Rs    Q      wr
          [2000,  1,    1.0e11],
          [2500,  3,    2.2e11],
          [2500,  4.5,  3.6e11],
          [2000,  1.0,  5.0e11],
          [2000,  4.0,  8.7e11],
          [6500,  1.3, 13.0e11],
          [30000, 0.7,  45.0e11],
]
wake = coll.Wake_t()
wake.Wl.resonator = True
for Rsi, Qi, wri in resonators:
    wake.Wl.wr.push_back(wri)
    wake.Wl.Rs.push_back(Rsi)
    wake.Wl.Q.push_back(Qi)

ring_aux = si.create_ring(phase=0)
ring_aux.cur_bun = np.linspace(0, 1, 21)*1e-3
currents = ring_aux.cur_bun
init_espread = ring_aux.espread(currents)
name = 'test_conv'

t0 = time.time()
currents, espread, bl, s0, nus_ave, nus_std, dists = \
    perform_time_domain_simulation(rin=ring,
                                   wake=wake,
                                   currents=currents,
                                   init_espread=init_espread,
                                   name=name,
                                   niter=100)
print('elapsed time: {0:15.4f}'.format(time.time()-t0))
