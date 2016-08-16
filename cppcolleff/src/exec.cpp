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
//     double s;
//     const long num_part = 10;
//
// }

int main()
{
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > s_;
    double s;
    const long num_part = 10000000;
    const long nturns   = 100;

    Ring_t ring;
    ring.energy      = 3e9;
    ring.tunex       = 0.13;
    ring.emitx       = 240e-12;
    ring.espread     = 9e-4;
    ring.bunlen      = 3e-3;
    ring.circum      = 518.396;
    ring.T0          = 518.396/light_speed;
    ring.mom_comp    = 1.7e-4;
    ring.harm_num    = 864;

    double phi0 (171.24/180*TWOPI), V0 (3e6);
    for (long i=-1000; i<=1000; i++){
        s = ((float )i ) * 1e-3 * 15e-2;
        ring.cav_s.push_back(s);
        ring.cav_V.push_back(3e6*sin(phi0 - ) * TWOPI*864/ring.circum*s);
    }

    Bunch_t bun (num_part,2e-3); //number of particles and current in A;
    generate_bunch(ring, bun);
    bun.sort();


    Wake_t wake;
    wake.Wl.resonator = true;
    wake.Wl.wr.push_back(10e9*TWOPI);
    wake.Wl.Rs.push_back(1e4);
    wake.Wl.Q.push_back(1);
    Feedback_t fb;
    Results_t results (nturns);

    std::chrono::time_point<clock_> beg_ = clock_::now();
    do_tracking(ring,wake,fb,bun,results);
    cout << chrono::duration_cast<s_> (clock_::now()-beg_).count() << endl;
    return 0;
}
