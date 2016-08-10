
#include <cppcolleff/Ring.h>

void Ring_t::track_one_turn(Bunch_t& bun)
{
    // the ring model is composed by one cavity,
    // then the magnetic elements and finally the wake-field.*/
     // if ((n % 1000000) == 0 || n == 0) {
     //    fprintf(stdout,"%09d ",n);
     //    for (int i=0;i<num_part;i++) {fprintf(stdout,"%12.6f %12.6f     ",ss[i]*1e3, de[i]*1e2);}
     //    fprintf(stdout,"\n");
     // }
    for (int p=0;p<bun.num_part;p++){
        // Longitudinal tracking:
        // bun.de[p] += 3e6*sin(171.24/360*TWOPI) * TWOPI*864/ring.circum*bun.ss[p]/ring.energy;
        bun.de[p] += interpola(cav_s, cav_V, bun.ss[p])/energy;
        bun.ss[p] += circum * mom_comp * bun.de[p];
        //Transverse tracking:
        bun.xx[p] -= etax  * bun.de[p]; // subtract the energy dependent fixed point;
        bun.xl[p] -= etaxl * bun.de[p];
        // calculate the invariant and the phase advance
        double phix  = TWOPI*(
         tunex +
         chromx*bun.de[p] +
         tunex_shift*(gammax*bun.xx[p]*bun.xx[p] + 2*alphax*bun.xx[p]*bun.xl[p] + betax*bun.xl[p]*bun.xl[p]) //Jx
        );
        double sinx  = sin(phix);
        double cosx  = cos(phix);
        // apply the one turn matrix
        double x_tmp  = bun.xx[p]*(cosx + alphax*sinx) + betax*bun.xl[p]*sinx;
        bun.xl[p]  =   -bun.xx[p]*gammax*sinx + bun.xl[p]*(cosx - alphax*sinx);
        bun.xx[p]  = x_tmp;

        bun.xx[p] += etax  * bun.de[p]; // add the energy dependent fixed point;
        bun.xl[p] += etaxl * bun.de[p];
    }
}
