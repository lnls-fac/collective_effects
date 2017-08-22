
#include <cppcolleff/Ring.h>

my_Dvector Ring_t::_get_distribution(const my_Dvector& spos, const my_Dvector& V) const
{
    double scale (mom_comp * circum * espread*espread); //The potential is already normalized by the energy;
    // First I integrate the potential to get the potential well:
    my_Dvector pot_well (spos.size(),0.0);
    double min_ipot = pot_well.back();
    for (int i=spos.size()-2;i>=0;--i){
        pot_well[i] = pot_well[i+1]  +  (V[i] + V[i+1]) * (spos[i+1]-spos[i]) / 2;
        if (pot_well[i] < min_ipot) min_ipot = pot_well[i];
    }

    // Now I compute the distribution
    my_Dvector distri (spos.size(),0.0);
	double sumdist (0.0);
    distri[0] = exp(  -(pot_well[0]-min_ipot) / scale  );
  #ifdef OPENMP
    #pragma omp parallel for schedule(guided,1) reduction(+:sumdist)
  #endif
    for (int i=1;i<spos.size();++i) {
		distri[i] = exp(  -(pot_well[i]-min_ipot) / scale  );
		sumdist += (distri[i]+distri[i-1]) * (spos[i]-spos[i-1]) / 2;
	}
  #ifdef OPENMP
    #pragma omp parallel for schedule(guided,1)
  #endif
	for (int i=0;i<spos.size();++i){distri[i] /= sumdist;} //normalization of distribution;

    return distri;
}
my_Dvector Ring_t::get_distribution() const
{
    const my_Dvector& spos = cav.ref_to_xi();
    const my_Dvector& V = cav.ref_to_yi();
    return _get_distribution(spos,V);
}
my_Dvector Ring_t::get_distribution(const my_Dvector& V) const
{
    const my_Dvector& spos = cav.ref_to_xi();
    return _get_distribution(spos,V);
}
my_Dvector Ring_t::get_distribution(const my_Dvector& spos, const my_Dvector& V) const
{
    return _get_distribution(spos,V);
}


Interpola_t Ring_t::_get_integrated_distribution(const my_Dvector& spos, const my_Dvector& V) const
{
    my_Dvector distr (get_distribution(spos,V));
	//need to resize the vectors;
	my_Dvector idistr, s_idistr;
	idistr.push_back(0.0);
	s_idistr.push_back(spos[0]);
	for (int i=1;i<spos.size();++i){
		double&& idistri = (distr[i]+distr[i-1]) * (spos[i]-spos[i-1]) / 2;
		if (idistri >= 1e-15) { // much more than the number of particles; think about it!
			idistr.push_back(idistr.back() + idistri); // for interpolation to work properly there must
			s_idistr.push_back(spos[i]); // not exist repeated values in the integrated distribution.
		}
	}
	return Interpola_t (idistr,s_idistr);
}
Interpola_t Ring_t::get_integrated_distribution() const
{
    const my_Dvector& spos = cav.ref_to_xi();
    const my_Dvector& V    = cav.ref_to_yi();
	return _get_integrated_distribution(spos,V);
}
Interpola_t Ring_t::get_integrated_distribution(const my_Dvector& V) const
{
    const my_Dvector& spos = cav.ref_to_xi();
	return _get_integrated_distribution(spos,V);
}
Interpola_t Ring_t::get_integrated_distribution(const my_Dvector& spos, const my_Dvector& V) const
{
	return _get_integrated_distribution(spos,V);
}


void Ring_t::track_one_turn(Bunch_t& bun) const
{
    double gammax ((1+alphax*alphax)/betax);
    // the ring model is composed by one cavity,
    // then the magnetic elements and finally the wake-field.*/
    my_PartVector& p = bun.particles;

  #ifdef OPENMP
    #pragma omp parallel for schedule(guided,1)
  #endif
    for (int i=0;i<p.size();++i){
        // Longitudinal tracking:
        // The potential is already normilized by the energy!
        p[i].de += cav.get_y(p[i].ss);
        p[i].ss += circum * mom_comp * p[i].de;
        // Transverse tracking:

        // subtract the energy dependent fixed point;
        p[i].xx -= etax  * p[i].de;
        p[i].xl -= etaxl * p[i].de;

        // calculate the invariant and the phase advance
        double&& phix  = TWOPI*(
            tunex +
            chromx*p[i].de +
            tunex_shift*(gammax*p[i].xx*p[i].xx +
                         2*alphax*p[i].xx*p[i].xl +
                         betax*p[i].xl*p[i].xl)  //Jx
        );
        double&& sinx  = sin(phix);
        double&& cosx  = cos(phix);
        // double&& cosx  = sqrt(1-sinx*sinx);

        // apply the one turn matrix
        double&& x_tmp  = p[i].xx*(cosx + alphax*sinx) + betax*p[i].xl*sinx;
        p[i].xl       =-p[i].xx*gammax*sinx + p[i].xl*(cosx - alphax*sinx);
        p[i].xx  = x_tmp;

        p[i].xx += etax  * p[i].de; // add the energy dependent fixed point;
        p[i].xl += etaxl * p[i].de;
    }
}
