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
    const long num_part = 1000000;
    const long nturns   = 400;

    Ring_t ring;
    ring.energy      = 3e9;
    ring.tunex       = 0.13;
    ring.emitx       = 250e-12;
    ring.espread     = 8e-4;
    ring.bunlen      = 3e-3;
    ring.circum      = 518.396;
    ring.T0          = 518.396/light_speed;
    ring.mom_comp    = 1.7e-4;
    ring.harm_num    = 864;
    ring.betax       = 19;
    ring.gammax      = 1/ring.betax;

    double phi0 (171.24/180*TWOPI), V0 (3e6);
    my_Dvector ss ;
    my_Dvector V  ;
    for (long i=-1000; i<=1000; i++){
        double&& s = 1e-3 * 10e-2 * i;
        // fprintf(stdout,"%f ",s);
        ss.push_back(s);
        // fprintf(stdout,"%f ",ss.back());
        V.push_back(V0*sin(phi0 -TWOPI*ring.harm_num/ring.circum * s) - V0*sin(phi0));
    }
    ring.cav.set_xy(ss,V);

    Bunch_t bun (num_part,1e-3); //number of particles and current in A;
    generate_bunch(ring, bun);
    bun.sort();


    Wake_t wake;
    wake.Wl.resonator = true;
    wake.Wl.wr.push_back(10e9*TWOPI);
    wake.Wl.Rs.push_back(1e2);
    wake.Wl.Q.push_back(1);
    Feedback_t fb;
    Results_t results (nturns);

    std::chrono::time_point<clock_> beg_ = clock_::now();
    do_tracking(ring,wake,fb,bun,results);
    cout << chrono::duration_cast<s_> (clock_::now()-beg_).count() << endl;
    return 0;
}
