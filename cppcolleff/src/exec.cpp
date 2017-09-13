// #include <string>
// #include <iostream>
// #include <cstdio>
// #include <cppcolleff/Bunch.h>
// #include <cppcolleff/essentials.h>
#include <chrono>  //std::high_resolution_clock
#include <iostream> //std::
#include <cppcolleff/cppcolleff.h>

int main()
{
    set_num_threads(32);
    set_seed_num(5004930);

    Ring_t ring = Ring_t();
    ring.energy = 3e9;
    ring.en_lost_rad = 470e3/ring.energy;
    ring.tunex = 0.13;
    ring.chromx = 0;
    ring.emitx = 250e-12;
    ring.espread = 8e-4;
    ring.circum = 518.396;
    ring.mom_comp = 1.7e-4;
    ring.harm_num = 864;
    ring.betax = 19;
    ring.etax = 0;
    ring.damp_nume = 1.7;
    ring.damp_numx = 1.3;


    double&& V0 = 3e6/ring.energy;
    double&& phi0 = TWOPI/2 - asin(ring.en_lost_rad/V0);
    double&& krf = TWOPI*ring.harm_num/ring.circum;

    my_Dvector ss;
    my_Dvector V;
    for (long i=-10000; i<=10000; ++i){
        double&& s = 1e-4 * 10e-2 * i;
        ss.push_back(s);
        V.push_back(V0*sin(phi0 + krf*s) - ring.en_lost_rad);
    }
    ring.cav.set_xy(ss, V);

    long num_part = 500000;  // takes 21 seconds with 32 processors.
    long nturns = 100;
    Bunch_t bun (num_part, 3e-3);
    bun.generate_particles(ring);
    bun.sort();

    double Rs[7] = {2000, 2500, 2500, 2000, 2000, 6500, 30000};
    double Q[7] = {1, 3, 4.5, 1.0, 4.0, 1.3, 0.7};
    double wr[7] = {1.0e11, 2.2e11, 3.6e11, 5.0e11, 8.7e11, 13.0e11, 45.0e11};

    Wake_t wake;
    for (long i=0; i<7; ++i){
        wake.Wl.wr.push_back(wr[i]);
        wake.Wl.Rs.push_back(Rs[i]);
        wake.Wl.Q.push_back(Q[i]);
    }
    my_Dvector x;
    int npt (100000);
    for (long i=-npt/2; i<npt/2; ++i)
        x.push_back(15e-2 * i/(npt/2));
    my_Dvector y (wake.Wl.get_wake_at_points(x, 1));
    wake.Wl.WF.set_xy(x, y);
    wake.Wl.WFC.prepare(x, y, false);
    wake.Wl.resonator = false;
    wake.Wl.wake_function = !wake.Wl.resonator;

    Feedback_t fb;

    Results_t results (nturns, 1);
    // results.set_print_every(2);
    results.set_keepWl(false);
    results.set_save_distributions_every(10);
    results.save_distribution_de = false;
    results.save_distribution_ss = false;
    results.bins[2] = 5000;
    results.bins[3] = 5000;

    // keep track of some particles
    my_Ivector indcs;
    indcs.push_back(900);
    indcs.push_back(1);
    bun.set_track_indcs(indcs);
    results.set_nparticles_to_track(indcs.size());

    typedef chrono::high_resolution_clock clock_;
    typedef chrono::duration<double, ratio<1> > s_;
    chrono::time_point<clock_> beg_ = clock_::now();
    single_bunch_tracking(ring,wake,fb,bun,results);
    cout << "ET: " << chrono::duration_cast<s_> (clock_::now()-beg_).count() << " s" << endl;
    return 0;
}
