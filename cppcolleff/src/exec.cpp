// #include <string>
// #include <iostream>
// #include <cstdio>
// #include <cppcolleff/Bunch.h>
// #include <cppcolleff/essentials.h>
#include <chrono>  //std::high_resolution_clock
#include <iostream> //std::
#include <cppcolleff/cppcolleff.h>

// int main()
// {
//     const long num_part = 10000000;
//     const long nturns   = 100;
//
//     Ring_t ring;
//     ring.energy      = 3e9;
//     ring.tunex       = 0.13;
//     ring.emitx       = 250e-12;
//     ring.espread     = 8e-4;
//     ring.circum      = 518.396;
//     ring.T0          = 518.396/light_speed;
//     ring.mom_comp    = 1.7e-4;
//     ring.harm_num    = 864;
//     ring.betax       = 19;
//
//     my_Dvector ss ;
//     my_Dvector V  ;
//     double phi0 (171.24/360*TWOPI), V0 (3e6), krf (TWOPI*ring.harm_num/ring.circum);
//     double U0 (V0*sin(phi0));
//     for (long i=-10000; i<=10000; i++){
//         double&& s = 1e-4 * 10e-2 * i;
//         ss.push_back(s);
//         V.push_back(V0*sin(phi0 + krf*s) - U0);
//     }
//     ring.cav.set_xy(ss,V);
//
//     Bunch_t bun (num_part,1e-3); //number of particles and current in A;
//
//     set_num_threads(8);
//     generate_bunch(ring, bun);
//     typedef std::chrono::high_resolution_clock clock_;
//     typedef std::chrono::duration<double, std::ratio<1> > s_;
//     std::chrono::time_point<clock_> beg1_ = clock_::now();
//     bun.sort();
//     cout << chrono::duration_cast<s_> (clock_::now()-beg1_).count() << endl;
// }

int main()
{

    Ring_t ring;
    ring.energy      = 3e9;
    ring.tunex       = 0.13;
    ring.emitx       = 250e-12;
    ring.espread     = 8e-4;
    ring.circum      = 518.396;
    ring.T0          = 518.396/light_speed;
    ring.mom_comp    = 1.7e-4;
    ring.harm_num    = 864;
    ring.betax       = 19;

    my_Dvector ss ;
    my_Dvector V  ;
    double phi0 (171.24/360*TWOPI), V0 (3e6/ring.energy), krf (TWOPI*ring.harm_num/ring.circum);
    double U0 (V0*sin(phi0));
    for (long i=-10000; i<=10000; i++){
        double&& s = 1e-4 * 10e-2 * i;
        ss.push_back(s);
        V.push_back(V0*sin(phi0 + krf*s) - U0);
    }
    ring.cav.set_xy(ss,V);
    ring.en_lost_rad = U0;

    cout << ring.en_lost_rad << endl;

    const long num_part = 10000000;
    const long nturns   = 2000;
    Bunch_t bun (num_part,1e-3); //number of particles and current in A;
    generate_bunch(ring, bun);
    bun.sort();


    Wake_t wake;
    // wake.Wl.resonator = false;
    // wake.Wl.wr.push_back(30e9*TWOPI);
    // wake.Wl.Rs.push_back(1e4);
    // wake.Wl.Q.push_back(1);

    // my_Dvector x;
    // for (long i=0; i<=50000; i++){x.push_back(2e-5 * 5e-2 * i); }
    // my_Dvector y(wake.Wd.get_wake_at_points(x,1));
    // wake.Wl.WF.set_xy(x,y);
    // wake.Wl.resonator = false;
    // wake.Wd.wake_function  = true;
    Feedback_t fb;
    Results_t results (nturns, 10);

    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > s_;
    // my_Dvector&& dist = ring.get_distribution();
    std::chrono::time_point<clock_> beg_ = clock_::now();
    NumThreads::set_num_threads(32);
    single_bunch_tracking(ring,wake,fb,bun,results);
    // for (double i=1;i<=10;++i) {
    //     double espread (find_equilibrium_energy_spread(wake,ring, 1e-3 * i));
    //     fprintf(stdout,"%9.5g\n",espread);
    // }
    // convolution_same(dist,Wl);
    cout << "ET: " << chrono::duration_cast<s_> (clock_::now()-beg_).count() << " s" << endl;
    return 0;
}
