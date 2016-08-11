
#include <cppcolleff/Ring.h>

void Ring_t::track_one_turn(Bunch_t& bun) const
{
    // the ring model is composed by one cavity,
    // then the magnetic elements and finally the wake-field.*/
     // if ((n % 1000000) == 0 || n == 0) {
     //    fprintf(stdout,"%09d ",n);
     //    for (int i=0;i<num_part;i++) {fprintf(stdout,"%12.6f %12.6f     ",ss[i]*1e3, de[i]*1e2);}
     //    fprintf(stdout,"\n");
     // }
    for (auto p=bun.particles.cbegin();p!=bun.particles.cend();++p){
        // Longitudinal tracking:
        // bun.de += 3e6*sin(171.24/360*TWOPI) * TWOPI*864/ring.circum*bun.ss/ring.energy;
        p->de += interpola(cav_s, cav_V, p->ss)/energy;
        p->ss += circum * mom_comp * p->de;
        //Transverse tracking:
        p->xx -= etax  * p->de; // subtract the energy dependent fixed point;
        p->xl -= etaxl * p->de;
        // calculate the invariant and the phase advance
        double phix  = TWOPI*(
         tunex +
         chromx*p->de +
         tunex_shift*(gammax*p->xx*p->xx + 2*alphax*p->xx*p->xl + betax*p->xl*p->xl) //Jx
        );
        double sinx  = sin(phix);
        double cosx  = cos(phix);
        // apply the one turn matrix
        double x_tmp  = p->xx*(cosx + alphax*sinx) + betax*p->xl*sinx;
        p->xl         =-p->xx*gammax*sinx + p->xl*(cosx - alphax*sinx);
        p->xx  = x_tmp;

        p->xx += etax  * p->de; // add the energy dependent fixed point;
        p->xl += etaxl * p->de;
    }
}
