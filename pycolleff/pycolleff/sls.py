#!/usr/bin/env python-sirius
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import mathphys as _mp
import simulate_landau as landau_cav

c = _mp.constants.light_speed


def main():

    ring = landau_cav.Ring()

    # SLS
    ring.E0 = 2.4e9
    ring.en_spread = 9e-4
    ring.frf = 499.632e6
    ring.harm_num = 480
    ring.mom_cmpct = 6e-4
    ring.en_lost_rad = 530e3
    ring.peak_rf = 2.1e6
    ring.nom_cur = 400e-3

    r = ring.en_lost_rad/ring.peak_rf

    wrf = ring.wrf
    krf = wrf/c
    h = ring.harm_num
    It = ring.nom_cur
    sigz = ring.bunlen

    Q = 2e8
    hc = landau_cav.HarmCav(wrf, n=3, Q=Q, Rs=88.4*Q, r=r)  # SLS

    F = np.exp(-(sigz*hc.num*wrf/c)**2)

    hc.calc_flat_potential(ring, F=F)

    hc.dw = 70.251e3 * 2 * np.pi

    zlim = 30*sigz
    npoints = 1501
    z = np.linspace(-zlim, zlim, npoints)

    epsilon = 1e-5
    param_conv = 10
    n_iters = 5000

    Ib = np.zeros(h, dtype=float)
    s_fill = h
    n_trains = 1
    s_gap = (h - n_trains * s_fill) // n_trains
    Ib[:] = It / s_fill / n_trains

    for j in range(n_trains):
        Ib[j * (s_fill + s_gap) + s_fill:(j + 1) * (s_fill+s_gap)] = 0

    lamb = landau_cav.Lambda(z, Ib, ring)

    Vrf = ring.get_Vrf(z)
    dist_old = lamb.get_gauss_dist(hc.length_shift, sigz)
    lamb.dist = np.array(dist_old)

    V_i = landau_cav.get_potentials_imp(ring, hc, lamb)
    Vt_imp = Vrf[None, :] + V_i
    dist_old = lamb.get_dist(Vt_imp, ring)

    # V_w = get_potentials_wake(ring, hc, lamb)
    # Vt_wake = Vrf[None, :] + V_w
    # dist_old = lamb.get_dist(Vt_wake, ring)

    plt.figure(figsize=(10, 14))
    gs = gridspec.GridSpec(4, 1)
    gs.update(left=0.10, right=0.95, bottom=0.10,
              top=0.97, wspace=0.35, hspace=0.25)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0], sharex=ax2)
    ax4 = plt.subplot(gs[3, 0], sharex=ax2)

    dist_new = np.zeros(dist_old.shape)

    for i in range(n_iters):
        lamb.dist = np.array(dist_old)
        lamb.form_fact = np.array(F)
        F_new = np.trapz(lamb.dist*np.exp(1j*hc.num*wrf*z/c), z)
        F_new /= np.trapz(lamb.dist, z)
        F_mod = np.absolute(np.mean(F_new))

        V_i = landau_cav.get_potentials_imp(ring, hc, lamb)
        Vt_imp = Vrf[None, :] + V_i
        dist_new = lamb.get_dist(Vt_imp, ring)

        # V_w = get_potentials_wake(ring, hc, lamb)
        # Vt_wake = Vrf[None, :] + V_w
        # dist_new = lamb.get_dist(Vt_wake, ring)

        dif = dist_new - dist_old
        res = np.trapz(np.absolute(dif), z)
        conv = np.mean(res)
        print('{0:03d}: {1:f}, F: {2:f}'.format(i, conv, F_mod))

        if conv < epsilon:
            break

        dist_old *= param_conv
        dist_old += dist_new
        dist_old /= param_conv + 1

    lamb.dist = np.array(dist_new)
    sigma_z_imp = lamb.get_bun_lens()
    bl_imp = np.mean(sigma_z_imp)
    z_ave_i = lamb.get_synch_phases()
    z_ave_ave_i = np.mean(z_ave_i)

    print('IMPEDANCE')
    hc.print_param(ring)

    V_i = landau_cav.get_potentials_imp(ring, hc, lamb)
    Vt_imp = Vrf + V_i[0, :]

    # V_w = get_potentials_wake(ring, hc, lamb)
    # Vt_wake = Vrf + V_w[0, :]

    # fwhm = 2*np.sqrt(2*np.log(2))   # To calculate FWHM Bunch Length for Gaussian distribution

    print('sync phase: {0:7.3f} mm'.format(z_ave_ave_i*1e3))
    print('bun length: {0:7.3f} mm ({1:7.3f} ps)'.format(bl_imp*1e3,
                                                         bl_imp*1e12/c))
    print('=============================')

    # print('ANALYTICAL')
    #
    # # hc.shunt_imp = 2.017e6
    # # hc.psi = 103.717*np.pi/180
    # hc.form_fact = F_mod
    #
    # V_a = get_potential_analytic_uni_fill(ring, hc, z)
    # Vt_a = Vrf + V_a
    #
    # lamb.dist = lamb.get_dist(Vt_a, ring)
    # sigma_z_a = lamb.get_bun_lens()
    # bl_a = np.mean(sigma_z_a)
    # z_ave_a = lamb.get_synch_phases()
    # z_ave_ave_a = np.mean(z_ave_a)
    #
    # print('sync phase analy: {0:7.3f}'.format(z_ave_ave_a*1e3))
    # print('bun length analy: {0:7.3f} mm ({1:7.3f} ps)'.format(bl_a*1e3,
    #                                                            bl_a*1e12/c))

    ph = z*krf
    ax1.plot(ph, dist_new[0, :], label='Distribution 1st bunch')

    ax2.plot(lamb.cur, label='Current')

    mask = lamb.cur < 1e-5
    sigma_z_imp[mask] = np.nan
    z_ave_i[mask] = np.nan

    ax3.plot(sigma_z_imp/ring.bunlen, label='Bunch Lengthening factor')
    ax4.plot(ring.synch_phase + z_ave_i*krf, label='Synch Phase')

    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax3.legend(loc='best')
    ax4.legend(loc='best')
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax4.grid(True)
    plt.show()


if __name__ == "__main__":
    landau_cav.memory_limit()   # Limitates maximum memory usage to half
    try:
        main()
    except MemoryError:
        sys.stderr.write('\n\nERROR: Memory Exception\n')
sys.exit(1)
