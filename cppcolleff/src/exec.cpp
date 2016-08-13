#include <string>
#include <iostream>
#include <cstdio>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/essentials.h>

// int main()
// {
//     double s;
//     const long num_part = 1000000;
//     const long nturns   = 100;
//
//     // Bunch_t bun (num_part);
//     // bun.generate_bunch();
//
//     for (int i=0;i<10;++i){
//         char test[14];
//         sprintf(test, "turn%04i.txt",i);
//         fprintf(stdout,"%s\n",test);
//         // bun.sort();
//     }
//     return 0;
// }

int main()
{
    double s;
    const long num_part = 1000000;
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

    Bunch_t bun (num_part);
    generate_bunch(ring, bun);
    bun.sort();


    Wake_t wake;
    Feedback_t fb;
    Results_t results (nturns);

    for (int i=0;i<num_part;i++){
        bun.xx[i] = ((1.0*i)/num_part -0.5)*5.0e-3;
        bun.xl[i] = 0.0;
        bun.de[i] = 0.0;
        bun.ss[i] = ((1.0*i)/num_part -0.5)*5.0e-3;
    }

    fb.track     = false;


    for (long i=-1000; i<=1000; i++){
        s = ((float )i ) * 1e-3 * 15e-2;
        ring.cav_s.push_back(s);
        ring.cav_V.push_back(3e6*cos(171.24/360*TWOPI) * TWOPI*864/ring.circum*s);
    }

    do_tracking(ring,wake,fb,bun,results);
    return 0;
}
