
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

    // #ifdef OPENMP
    //   #pragma omp parallel for schedule(guided,1)
    // #endif
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
        results.emplace_back(pool.enqueue(
            &Ring_t::_track_one_turn, this, ref(p), seed, lims[i], lims[i+1]
        ));
    }
    for(auto && result: results) result.get();

    // When omp is preferred uncomment this option
    // _track_one_turn(p, seed, 0, p.size());

}

void Ring_t::track_one_turn(Bunch_t& bun, int n) const
{
    ThreadPool pool (get_num_threads());
    return track_one_turn(bun, n, pool);
}


void Ring_t::to_file(const char* filename) const
{
	ofstream fp(filename);
	if (fp.fail()) exit(1);
	fp.setf(fp.left | fp.scientific);
	fp.precision(15);
    fp << setw(30) << "% energy" << energy << " eV" << endl;
    fp << setw(30) << "% circumference" << circum << " m" << endl;
    fp << setw(30) << "% harmonic_number" << harm_num << endl;
    fp << setw(30) << "% momentum_compaction" << mom_comp << endl;
	fp << setw(30) << "% betax" << betax << " m" << endl;
    fp << setw(30) << "% etax" << etax << " m" << endl;
    fp << setw(30) << "% alphax" << alphax << endl;
    fp << setw(30) << "% etaxl" << etaxl << endl;
    fp << setw(30) << "% tunex" << tunex << endl;
    fp << setw(30) << "% chromx" << chromx << endl;
    fp << setw(30) << "% tune_shiftx" << tunex_shift << endl;
    fp << setw(30) << "% frac_energy_lost_per_turn" << en_lost_rad << endl;
    fp << setw(30) << "% energy_spread" << espread << endl;
    fp << setw(30) << "% emittancex" << emitx << " m.rad" << endl;
    fp << setw(30) << "% damping_part_num_long" << damp_nume << endl;
    fp << setw(30) << "% damping_part_num_x" << damp_numx << endl;
    auto& s = cav.ref_to_xi();
    auto& V = cav.ref_to_yi();
    fp << setw(30) << "% cav_num_points" << s.size() << endl;
    fp << setw(26) << "# ss [m]";
    fp << setw(26) << "(V - U0)/E" << endl;
    fp.setf(fp.left | fp.showpos | fp.scientific);
    for (auto i=0; i<s.size(); ++i){
        fp << setw(26) << s[i];
        fp << setw(26) << V[i] << endl;
    }
    fp.close();
}

void Ring_t::from_file(const char* filename)
{
	ifstream fp(filename);
	if (fp.fail()) return;

    my_Dvector pos, V;
	double s(0.0), v(0.0);
  	string line;
	unsigned long line_count = 0;
	while (getline(fp, line)) {
  		line_count++;
  		istringstream ss(line);
		char c = ss.get();
		while (c == ' ') c = ss.get();
  		if (c == '#' || c == '\n') continue;
  		if (c == '%') {
			string cmd;
	  		ss >> cmd;
            if (cmd.compare("energy") == 0) {ss >> energy; continue;}
            if (cmd.compare("circumference") == 0) {ss >> circum; continue;}
            if (cmd.compare("harmonic_number") == 0) {ss >> harm_num; continue;}
            if (cmd.compare("momentum_compaction") == 0) {ss >> mom_comp; continue;}
        	if (cmd.compare("betax") == 0) {ss >> betax; continue;}
            if (cmd.compare("etax") == 0) {ss >> etax; continue;}
            if (cmd.compare("alphax") == 0) {ss >> alphax; continue;}
            if (cmd.compare("etaxl") == 0) {ss >> etaxl; continue;}
            if (cmd.compare("tunex") == 0) {ss >> tunex; continue;}
            if (cmd.compare("chromx") == 0) {ss >> chromx; continue;}
            if (cmd.compare("tune_shiftx") == 0) {ss >> tunex_shift; continue;}
            if (cmd.compare("frac_energy_lost_per_turn") == 0) {ss >> en_lost_rad; continue;}
            if (cmd.compare("energy_spread") == 0) {ss >> espread; continue;}
            if (cmd.compare("emittancex") == 0) {ss >> emitx; continue;}
            if (cmd.compare("damping_part_num_long") == 0) {ss >> damp_nume; continue;}
            if (cmd.compare("damping_part_num_x") == 0) {ss >> damp_numx; continue;}
            if (cmd.compare("cav_num_points") == 0) {
                int np;
                ss >> np;
                pos.reserve(np);
                V.reserve(np);
                continue;
            }
  		}
		ss.unget();
  		ss >> s; ss >> v;
        pos.push_back(s);
        V.push_back(v);
	}
    cav.set_xy(pos, V);
	fp.close();
}
