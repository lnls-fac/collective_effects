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

    Ib = np.zeros(h, dtype=float)
    s_fill = 390
    n_trains = 1
    s_gap = (h - n_trains * s_fill) // n_trains
    Ib[:] = It / s_fill / n_trains

    for j in range(n_trains):
        Ib[j * (s_fill + s_gap) + s_fill:(j + 1) * (s_fill+s_gap)] = 0

    Ib[409] = 4e-3  # camshaft

    lamb = landau_cav.Lambda(z, Ib, ring)

    _, _, dist_new = landau_cav.calc_equilibrium_potential(ring, lamb, hc, z,
                                                           epsilon=1e-5,
                                                           param_conv=15,
                                                           n_iters=1000)
    lamb.dist = np.array(dist_new)
    sigma_z_imp = lamb.get_bun_lens()
    sigma_z_imp_fwhm = lamb.get_bun_lens_fwhm() / 2.35

    z_ave_i = lamb.get_synch_phases()
    bl_imp = np.mean(sigma_z_imp)
    z_ave_ave_i = np.mean(z_ave_i)

    print('IMPEDANCE')
    hc.print_param(ring)

    # fwhm = 2*np.sqrt(2*np.log(2))

    mask = lamb.cur < 1e-8
    sigma_z_imp[mask] = np.nan
    z_ave_i[mask] = np.nan
    ind_max = np.nanargmax(sigma_z_imp)

    print('sync phase: {0:7.3f} mm'.format(z_ave_ave_i*1e3))
    print('bun length: {0:7.3f} mm ({1:7.3f} ps)'.format(bl_imp*1e3,
                                                         bl_imp*1e12/c))
    print('max bun length factor: {0:1.3f} (Bunch number {1:1g})'.format(np.nanmax(sigma_z_imp)/sigz, ind_max))
    print('max bun length factor FWHM: {0:1.3f} (Bunch number {1:1g})'.format(np.nanmax(sigma_z_imp_fwhm)/sigz, ind_max))
    plt.figure(figsize=(10, 14))
    gs = gridspec.GridSpec(4, 1)
    gs.update(left=0.10, right=0.95, bottom=0.10,
              top=0.97, wspace=0.35, hspace=0.25)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0], sharex=ax2)
    ax4 = plt.subplot(gs[3, 0], sharex=ax2)

    ph = z*krf

    ax1.plot(ph, dist_new[0, :], label='Distribution 1st bunch')
    ax2.plot(lamb.cur, label='Current - [mA]')

    mask = lamb.cur < 1e-6
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
