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
//     typedef std::chrono::high_resolution_clock clock_;
//     typedef std::chrono::duration<double, std::ratio<1> > s_;
//
//     std::chrono::time_point<clock_> beg1_ = clock_::now();
//     for (long i=0;i<10;++i) {double a = 3*4e-5*i;fprintf(stdout,"%f\n",a);}
//     cout << chrono::duration_cast<s_> (clock_::now()-beg1_).count() << endl;
//     std::chrono::time_point<clock_> beg2_ = clock_::now();
//     for (long i=0;i<10;++i) {double&& a = 3*4e-5*i;fprintf(stdout,"%f\n",a);}
//     cout << chrono::duration_cast<s_> (clock_::now()-beg2_).count() << endl;
//     // fprintf(stdout,"%f %f \n",b[2],b[3]);
//     // fprintf(stdout,"%f \n",a[2]);
// }

int main()
{
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > s_;
    const long num_part = 10000000;
    const long nturns   = 100;

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
    double phi0 (171.24/360*TWOPI), V0 (3e6), krf (TWOPI*ring.harm_num/ring.circum);
    double U0 (V0*sin(phi0));
    for (long i=-10000; i<=10000; i++){
        double&& s = 1e-4 * 10e-2 * i;
        ss.push_back(s);
        V.push_back(V0*sin(phi0 + krf*s) - U0);
    }
    ring.cav.set_xy(ss,V);

    Bunch_t bun (num_part,1e-3); //number of particles and current in A;
    generate_bunch(ring, bun);
    bun.sort();


    Wake_t wake;
    wake.Wl.resonator = true;
    wake.Wl.wr.push_back(30e9*TWOPI);
    wake.Wl.Rs.push_back(1e4);
    wake.Wl.Q.push_back(1);

    my_Dvector x;
    for (long i=0; i<=50000; i++){x.push_back(2e-5 * 5e-2 * i); }
    my_Dvector y(wake.Wl.get_wake_at_points(x,1));
    wake.Wl.W.set_xy(x,y);
    // wake.Wl.resonator = false;
    // wake.Wl.general  = true;
    Feedback_t fb;
    Results_t results (nturns);

    // my_Dvector&& dist = ring.get_distribution();
    std::chrono::time_point<clock_> beg_ = clock_::now();
    set_num_threads(2);
    do_tracking(ring,wake,fb,bun,results);
    // solve_Haissinski(wake,ring,5e-3);
    // convolution_same(dist,Wl);
    cout << "ET: " << chrono::duration_cast<s_> (clock_::now()-beg_).count() << " s" << endl;
    return 0;
}
