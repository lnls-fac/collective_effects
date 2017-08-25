
#include <cppcolleff/Ring.h>

static ThreadVars ThreadInfo(false);

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

int Ring_t::_track_one_turn(
    my_PartVector& p,
    const unsigned int th,
    const int init,
    const int final_) const
{
  // the ring model is composed by the magnetic elements,
  // then one cavity and finally the wake-field.*/
    double gammax ((1+alphax*alphax)/betax);

    // For damping and quantum excitation calculations
    double Fde (1 - damp_nume * en_lost_rad);
    double Fxl (1 - damp_numx * en_lost_rad);  //with rf contribution
    double Srde (sqrt(1 - Fde*Fde) * espread);
    double Srxl (sqrt((1 - Fxl*Fxl) * emitx / betax));
    normal_distribution<double> distribution(0.0,1.0);
    default_random_engine gen1(th);

    for (int i=init;i<final_;++i){
        // subtract the energy dependent fixed point;
        p[i].xx -= etax * p[i].de;
        p[i].xl -= etaxl * p[i].de;

        // Longitudinal tracking:
        p[i].ss += circum * mom_comp * p[i].de;
        p[i].de *= Fde;  // Damping
        p[i].de += Srde * distribution(gen1);  // Quantum excitation
        p[i].de += cav.get_y(p[i].ss);  // Already normalized by the energy!
        // Transverse tracking:

        // calculate the invariant and the phase advance
        double&& phix = TWOPI*(
            tunex + chromx*p[i].de + tunex_shift*(gammax * p[i].xx*p[i].xx +
                                                  2*alphax * p[i].xx*p[i].xl +
                                                  betax * p[i].xl*p[i].xl)  //Jx
            );
        double&& sinx = sin(phix);
        double&& cosx = cos(phix);
        // double&& cosx = sqrt(1-sinx*sinx);

        // apply the one turn matrix
        double&& x_tmp = p[i].xx*(cosx + alphax*sinx) + betax*p[i].xl*sinx;
        p[i].xl = -p[i].xx*gammax*sinx + p[i].xl*(cosx - alphax*sinx);
        p[i].xx = x_tmp;
        // Damping and Quantum excitation simulations:
        p[i].xl *= Fxl;
        p[i].xl += Srxl * distribution(gen1);

        p[i].xx += etax  * p[i].de; // add the energy dependent fixed point;
        p[i].xl += etaxl * p[i].de;
    }
    return 1;
}

void Ring_t::track_one_turn(Bunch_t& bun, ThreadPool& pool, int n) const
{
    my_PartVector& p = bun.particles;
    unsigned int nr_th = ThreadInfo.get_num_threads();
	my_Ivector lims (ThreadInfo.get_bounds(0,p.size()));

    std::vector< std::future<int> > results;

    for (unsigned int i=0;i<nr_th;++i){
        results.emplace_back(
            pool.enqueue(
                &Ring_t::_track_one_turn, this, ref(p), n*nr_th + i, lims[i], lims[i+1]
            )
        );
    }
    for(auto && result: results) result.get();
}
