
#include <cppcolleff/Ring.h>


// my_Dvector Ring_t::get_distribution(
//     ThreadPool& pool,
//     const my_Dvector& V_,
//     const my_Dvector& spos_) const
// {
//     my_Dvector V (V_);
//     if (V.size() == 0) V = cav.ref_to_yi();
//     my_Dvector spos (spos_);
//     if (spos.size() == 0) spos = cav.ref_to_xi();
//
//
//     double scale (mom_comp * circum * espread*espread); //The potential is already normalized by the energy;
//     // First I integrate the potential to get the potential well:
//     my_Dvector pot_well (spos.size(),0.0);
//     double min_ipot = pot_well.back();
//     for (int i=spos.size()-2;i>=0;--i){
//         pot_well[i] = pot_well[i+1]  +  (V[i] + V[i+1]) * (spos[i+1]-spos[i]) / 2;
//         if (pot_well[i] < min_ipot) min_ipot = pot_well[i];
//     }
//
//     // Now I compute the distribution
//     my_Dvector distri (spos.size(),0.0);
// 	double sumdist (0.0);
//     distri[0] = exp(  -(pot_well[0]-min_ipot) / scale  );
//
//     for (int i=1;i<spos.size();++i) {
// 		distri[i] = exp(  -(pot_well[i]-min_ipot) / scale  );
// 		sumdist += (distri[i]+distri[i-1]) * (spos[i]-spos[i-1]) / 2;
// 	}
//
// 	for (int i=0;i<spos.size();++i){distri[i] /= sumdist;} //normalization of distribution;
//
//     return distri;
// }


my_Dvector Ring_t::get_distribution(
    ThreadPool& pool,
    const my_Dvector& V_,
    const my_Dvector& spos_) const
{
    my_Dvector V (V_);
    if (V.size() == 0) V = cav.ref_to_yi();
    my_Dvector spos (spos_);
    if (spos.size() == 0) spos = cav.ref_to_xi();

    double scale (mom_comp * circum * espread*espread); //The potential is already normalized by the energy;

    unsigned int nr_th = get_num_threads();
	my_Ivector lims1 (get_bounds(1,spos.size()));
    my_Ivector lims2 (get_bounds(0,spos.size()));
    vector< std::future<int> > res1;
    vector< std::future<double> > res2;
    vector< std::future<int> > res3;

    // First I integrate the potential to get the potential well:
    my_Dvector pot_well (spos.size(),0.0);
    double min_ipot = pot_well.back();
    for (int i=spos.size()-2;i>=0;--i){
        pot_well[i] = pot_well[i+1]  +  (V[i] + V[i+1]) * (spos[i+1]-spos[i]) / 2;
        if (pot_well[i] < min_ipot) min_ipot = pot_well[i];
    }

    // Now I compute the distribution
    auto fun1 = [&](my_Dvector& dis, int ini, int fin)
    {
        for (int ii=ini;ii<fin;++ii)
            {dis[ii] = exp(-(pot_well[ii]-min_ipot) / scale  );}
        return 1;
    };
    my_Dvector distri (spos.size(),0.0);
    distri[0] = exp(  -(pot_well[0]-min_ipot) / scale  );
    for (unsigned int i=0;i<nr_th;++i){
        res1.emplace_back(pool.enqueue(fun1, ref(distri), lims1[i], lims1[i+1]));
    }
    for(auto&& result: res1) result.get();

    // Calculate its normalization:
    auto fun2 = [&](int ini, int fin)
    {
        double sums (0.0);
        for (int ii=ini;ii<fin;++ii)
            {sums += (distri[ii]+distri[ii-1]) * (spos[ii]-spos[ii-1]) / 2;}
        return sums;
    };
    double sumdist (0.0);
    for (unsigned int i=0;i<nr_th;++i){
        res2.emplace_back(pool.enqueue(fun2, lims1[i], lims1[i+1]));
    }
    for(auto&& result: res2) sumdist += result.get();


    // And normalize it:
    auto fun3 = [&] (my_Dvector& dis, int ini, int fin)
    {
        for (int ii=ini;ii<fin;++ii) {dis[ii] /= sumdist;}
        return 1;
    };
    for (unsigned int i=0;i<nr_th;++i){
        res3.emplace_back(pool.enqueue(fun3, ref(distri), lims2[i], lims2[i+1]));
    }
    for(auto&& result: res3) result.get();

    return distri;
}
my_Dvector Ring_t::get_distribution(const my_Dvector& V, const my_Dvector& spos) const
{
    ThreadPool pool (get_num_threads());
    return get_distribution(pool,V,spos);
}


Interpola_t Ring_t::_get_integrated_distribution(const my_Dvector& spos, const my_Dvector& V) const
{
    my_Dvector distr (get_distribution(V,spos));
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
    const unsigned int seed,
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
    normal_distribution<double> Gauss(0.0,1.0);
    default_random_engine gen1(seed);

    for (int i=init;i<final_;++i){
        // subtract the energy dependent fixed point;
        p[i].xx -= etax * p[i].de;
        p[i].xl -= etaxl * p[i].de;

        // Longitudinal tracking:
        p[i].ss += circum * mom_comp * p[i].de;
        p[i].de *= Fde;  // Damping
        p[i].de += Srde * Gauss(gen1);  // Quantum excitation
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
        p[i].xl += Srxl * Gauss(gen1);

        p[i].xx += etax  * p[i].de; // add the energy dependent fixed point;
        p[i].xl += etaxl * p[i].de;
    }
    return 1;
}

void Ring_t::track_one_turn(Bunch_t& bun, int n, ThreadPool& pool) const
{
    extern unsigned long seed;

    my_PartVector& p = bun.particles;
    unsigned int nr_th = get_num_threads();
	my_Ivector lims (get_bounds(0,p.size()));

    std::vector< std::future<int> > results;

    for (unsigned int i=0;i<nr_th;++i){
        seed += 1;
        results.emplace_back(
            pool.enqueue(
                &Ring_t::_track_one_turn, this, ref(p), seed, lims[i], lims[i+1]
            )
        );
    }
    for(auto && result: results) result.get();
}

void Ring_t::track_one_turn(Bunch_t& bun, int n) const
{
    ThreadPool pool (get_num_threads());
    return track_one_turn(bun, n, pool);
}
