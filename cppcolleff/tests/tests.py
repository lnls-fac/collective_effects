#!/usr/bin/env python-sirius

import time
import math
import cppcolleff as coll

coll.set_num_threads(32)

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

wake = coll.Wake_t()
wake.Wl.resonator = False
wake.Wl.wr.push_back(30e9*coll.TWOPI)
wake.Wl.Rs.push_back(1e4)
wake.Wl.Q.push_back(1)

fb = coll.Feedback_t()

results = coll.Results_t(nturns, 100)

# bun.scale_longitudinal(0.1)
# bun.scale_transverse(0.1)
t0 = time.time()
coll.single_bunch_tracking(ring, wake, fb, bun, results)
print('elapsed time: {0:15.4f}'.format(time.time()-t0))
