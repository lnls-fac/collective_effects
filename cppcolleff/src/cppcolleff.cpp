
#include <cppcolleff/cppcolleff.h>

void generate_bunch(const Ring_t& ring, Bunch_t& bun)
{
	default_random_engine gen(19880419);
  	normal_distribution<double> espread_dist(0.0,ring.espread);
	exponential_distribution<double> emitx_dist(1/(2*ring.emitx));
	uniform_real_distribution<double> phix_dist(0.0,TWOPI);

    auto& cav_s = ring.cav.ref_to_xi();
    auto& cav_V = ring.cav.ref_to_yi();

    // First I integrate the potential to get the potential well:
    my_Dvector pot_well (cav_s.size(),0.0);
    double min_ipot = pot_well.back();
    for (int i=cav_s.size()-2;i>=0;--i){
        pot_well[i] = pot_well [i+1] + (cav_V[i] + cav_V[i+1])*(cav_s[i+1]-cav_s[i])/2;
        if (pot_well[i] < min_ipot) min_ipot = pot_well[i];
    }

    double&& exp_fact = 1/(ring.mom_comp*ring.circum*ring.energy)/(ring.espread*ring.espread);

    // Now I compute the equilibrium distribution and the integrated distribution to
    // be able to generate the particles
    my_Dvector idistr, s_distr;
    idistr.push_back(0.0);
    s_distr.push_back(cav_s[0]);
    for (int i=1;i<cav_s.size();++i){
        double&& distr = (exp(-exp_fact*(pot_well[i]  -min_ipot)) // I don't really need the distribution,
                        + exp(-exp_fact*(pot_well[i-1]-min_ipot)))*(cav_s[i]-cav_s[i-1])/2; // only the integrated.
        if (distr >= 1e-15) { // much more than the number of particles; think about it!
            idistr.push_back(idistr.back() + distr); // for interpolation to work properly the there must
            s_distr.push_back(cav_s[i]); // not exist repeated values in the integrated distribution.
        }
    }

    uniform_real_distribution<double> ss_dist(idistr.front(),idistr.back());
    Interpola_t idistr_interpol (idistr, s_distr);
	for (auto& p:bun.particles){
		double&& emitx = emitx_dist(gen);
		double&& phix  = phix_dist(gen);
		p.de = espread_dist(gen);
		p.ss = idistr_interpol.get_y(ss_dist(gen));
		p.xx =  sqrt(emitx*ring.betax)*cos(phix) + ring.etax*p.de;
		p.xl = -sqrt(emitx/ring.betax)*(ring.alphax*cos(phix) + sin(phix)) + ring.etaxl*p.de;
	}
	bun.is_sorted = false;
}


void do_tracking(
    const Ring_t& ring,
    const Wake_t& wake,
    Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results)
{
    //current dependent strength of the kick:
    const double kick_stren = ring.T0 / ring.energy * bun.Ib / bun.num_part;

    for (long n=0;n<results.get_nturns();n++){
        double&& xx_ave = results.calc_stats(bun,n);
        /* convention: positive ss, means particle behind the sinchronous particle;
         First do single particle tracking:*/
        ring.track_one_turn(bun);

        // After this sorting, the particles will be ordered from head to tail.
        // It means, from smaller ss to bigger ss.
        bun.general_sort();

        results.register_Wkicks(n, wake.apply_kicks(bun,kick_stren, ring.betax));
        results.register_FBkick(n,   fb.apply_kick(bun, xx_ave,     ring.betax));
    }
    results.calc_stats(bun,results.get_nturns());
}
