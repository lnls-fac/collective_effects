#!/usr/bin/env python-sirius


# In this example the bunch is tracked for 40000 turns 
# for each current in a range from 0.1 to 4 mA. 
# The final bunch of one tracking is used as input to
# the next step; only the current of the bunch is increased.
# The script will save all the parameters of the simulation
# in a way that is possible to reload them again with 
# functions of the cppcolleff library.
# It will create folder for each current and save the average
# and spread of each of the four dimensions of the bunch
# for each turn there. The longitudinal and energy distributions
# will also be saved as well as the initial and final bunch.
#
# After the simulation the results can be analised with
# the scripts analyse_results.py and analyse_distributions.py.

import time
import os
import math
import cppcolleff as coll
import numpy as np
import sys


# global parameters:
coll.set_num_threads(24)
coll.set_seed_num(1)
num_part = 1000000  
nturns = 40000  #five damping times
currents = [x*1e-4 for x in range(1,40)] # up to 4mA

# setup the ring
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
ring.to_file('ring.txt')


# setup the bunch
bun = coll.Bunch_t(num_part, 1e-3)
coll.generate_bunch(ring, bun)
bun.sort()


# setup the wake
resonators = [  # Rs    Q      wr
          [2000,  1,    1.0e11],
          [2500,  3,    2.2e11],
          [2500,  4.5,  3.6e11],
          [2000,  1.0,  5.0e11],
          [2000,  4.0,  8.7e11],
          [6500,  1.3, 13.0e11],
          [30000, 0.7,  45.0e11]]
wake = coll.Wake_t()
wake.Wl.resonator = True
for Rsi, Qi, wri in resonators:
    wake.Wl.wr.push_back(wri)
    wake.Wl.Rs.push_back(Rsi)
    wake.Wl.Q.push_back(Qi)
wake.to_file('wake.txt')


# setup the feedback (not used here.)
fb = coll.Feedback_t()


for cur in currents:
    print('\n\n' + 50*'#')
    print('Starting calculation for Ib = {0:5.3f} mA'.format(cur*1e3))
    print(50*'#')
    sys.stdout.flush()
    dname = "curr_{0:05.3f}mA".format(cur*1e3)
    os.mkdir(dname)
    os.chdir(dname)

    # setup results:
    results = coll.Results_t(nturns, 1)
    results.set_keepWl(True)
    results.set_print_every(100)
    results.set_save_distributions_every(215)  # one synchrotron period
    results.save_distribution_de = True
    results.save_distribution_ss = True
    results.bins[2] = 5000
    results.bins[3] = 5000

    # setup current, the bunch will be reused from one simulation to another:
    bun.Ib = cur
    bun.to_file('bunch_initial.txt')

    t0 = time.time()
    coll.single_bunch_tracking(ring, wake, fb, bun, results)
    print('elapsed time: {0:15.4f}'.format(time.time()-t0))
    print(50*'#')
    sys.stdout.flush()
    results.to_file('results.txt')
    bun.to_file('bunch_final.txt')
    os.chdir('../')
