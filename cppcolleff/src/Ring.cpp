
#include <cppcolleff/Ring.h>

void Ring_t::track_one_turn(Bunch_t& bun) const
{
  #ifdef OPENMP
    // the ring model is composed by one cavity,
    // then the magnetic elements and finally the wake-field.*/
    my_PartVector& par = bun.particles;

    #pragma omp parallel for schedule(guided,1)
    for (int i=0;i<par.size();++i){
        // Longitudinal tracking:
        par[i].de += cav.get_y(par[i].ss)/energy;
        par[i].ss += circum * mom_comp * par[i].de;
        //Transverse tracking:
        par[i].xx -= etax  * par[i].de; // subtract the energy dependent fixed point;
        par[i].xl -= etaxl * par[i].de;
        // calculate the invariant and the phase advance
        double&& phix  = TWOPI*(
            tunex +
            chromx*par[i].de +
            tunex_shift*(gammax*par[i].xx*par[i].xx + 2*alphax*par[i].xx*par[i].xl + betax*par[i].xl*par[i].xl) //Jx
        );
        double&& sinx  = sin(phix);
        double&& cosx  = cos(phix);
        // apply the one turn matrix
        double&& x_tmp  = par[i].xx*(cosx + alphax*sinx) + betax*par[i].xl*sinx;
        par[i].xl            =-par[i].xx*gammax*sinx + par[i].xl*(cosx - alphax*sinx);
        par[i].xx  = x_tmp;

        par[i].xx += etax  * par[i].de; // add the energy dependent fixed point;
        par[i].xl += etaxl * par[i].de;
    }
  #else
    // the ring model is composed by one cavity,
    // then the magnetic elements and finally the wake-field.*/
    for (auto& p:bun.particles){
        // Longitudinal tracking:
        p.de += cav.get_y(p.ss)/energy;
        p.ss += circum * mom_comp * p.de;
        //Transverse tracking:
        p.xx -= etax  * p.de; // subtract the energy dependent fixed point;
        p.xl -= etaxl * p.de;
        // calculate the invariant and the phase advance
        double&& phix  = TWOPI*(
         tunex +
         chromx*p.de +
         tunex_shift*(gammax*p.xx*p.xx + 2*alphax*p.xx*p.xl + betax*p.xl*p.xl) //Jx
        );
        double&& sinx  = sin(phix);
        double&& cosx  = cos(phix);
        // apply the one turn matrix
        double&& x_tmp  = p.xx*(cosx + alphax*sinx) + betax*p.xl*sinx;
        p.xl            =-p.xx*gammax*sinx + p.xl*(cosx - alphax*sinx);
        p.xx  = x_tmp;

        p.xx += etax  * p.de; // add the energy dependent fixed point;
        p.xl += etaxl * p.de;
    }
  #endif
}
