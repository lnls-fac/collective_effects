#!/usr/bin/env python-sirius

import time
import math
import cppcolleff as coll
import matplotlib.pyplot as plt
import numpy as np

coll.set_num_threads(32)
coll.set_seed_num(5004930)

ring = coll.Ring_t()
ring.energy = 3e9
ring.en_lost_rad = 470e3/ring.energy
ring.tunex = 0.13
ring.chromx = 0
ring.emitx = 250e-12
ring.espread = 8e-4
ring.circum = 518.396
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

num_part = 10000  # takes 21 seconds with 32 processors.
nturns = 100
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
# resonators = [  # Rs    Q      wr
#            [4000,  1,    1.0e11]]

wake = coll.Wake_t()
wake.Wl.resonator = False
for Rsi, Qi, wri in resonators:
    wake.Wl.wr.push_back(wri)
    wake.Wl.Rs.push_back(Rsi)
    wake.Wl.Q.push_back(Qi)

# x = coll.my_Dvector()
# for i in range(50000):
#     x.push_back(2e-5 * 5e-2 * i)
# y = coll.my_Dvector(wake.Wl.get_wake_at_points(x, 1))
# wake.Wl.WF.set_xy(x, y)
# wake.Wl.resonator = True
# wake.Wl.wake_function = False

fb = coll.Feedback_t()

results = coll.Results_t(nturns, 1)
results.set_keepWl(False)
results.set_save_distributions_every(10)
results.save_distribution_de = False
results.save_distribution_ss = False
results.set_nparticles_to_track(2)
indcs = coll.my_Ivector()
indcs.push_back(900)
indcs.push_back(1)
bun.set_track_indcs(indcs)
results.bins[2] = 5000
results.bins[3] = 5000

# bun.scale_longitudinal(0.1)
# bun.scale_transverse(0.1)
t0 = time.time()
coll.single_bunch_tracking(ring, wake, fb, bun, results)
print('elapsed time: {0:15.4f}'.format(time.time()-t0))
results.to_file('results.txt')
bun.to_file('bunch.txt')

# t0 = time.time()
# dists = []
# for i in range(10):
#     results = coll.Results_t(10, 1)
#     coll.single_bunch_tracking(ring, wake, fb, bun, results)
#     dists.append(bun.calc_particles_distribution(ss))
#     plt.plot(ss, dists[-1])
#     plt.xlim(-1.5e-2, 1.5e-2)
#     plt.show()
# print('elapsed time: {0:15.4f}'.format(time.time()-t0))
