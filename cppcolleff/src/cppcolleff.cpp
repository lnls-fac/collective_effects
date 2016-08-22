
#include <cppcolleff/cppcolleff.h>

void generate_bunch(const Ring_t& ring, Bunch_t& bun, Results_t& results)
{

    auto& cav_s = ring.cav.ref_to_xi();
    auto& cav_V = ring.cav.ref_to_yi();

	double scale (ring.mom_comp * ring.circum * ring.energy);
    my_Dvector&& distr = get_distribution_from_potential(cav_s,cav_V,ring.espread,scale);
	Interpola_t&& idistr_interpol = get_integrated_distribution(cav_s, distr);

	default_random_engine gen(19880419);
	normal_distribution<double> espread_dist(0.0,ring.espread);
	exponential_distribution<double> emitx_dist(1/(2*ring.emitx));
	uniform_real_distribution<double> phix_dist(0.0,TWOPI);
	uniform_real_distribution<double> ss_dist(0.0,1.0);
  #ifdef OPENMP
    fprintf(stdout,"here\n");
    my_PartVector& par = bun.particles;
    #pragma omp parallel for schedule(guided,1)
    for (int i=0;i<par.size();++i){
        double&& emitx = emitx_dist(gen);
        double&& phix  = phix_dist(gen);
        par[i].de = espread_dist(gen);
        par[i].ss = idistr_interpol.get_y(ss_dist(gen));
        par[i].xx =  sqrt(emitx*ring.betax)*cos(phix) + ring.etax*par[i].de;
        par[i].xl = -sqrt(emitx/ring.betax)*(ring.alphax*cos(phix) + sin(phix)) + ring.etaxl*par[i].de;
    }
  #else
	for (auto& p:bun.particles){
		double&& emitx = emitx_dist(gen);
		double&& phix  = phix_dist(gen);
		p.de = espread_dist(gen);
		p.ss = idistr_interpol.get_y(ss_dist(gen));
		p.xx =  sqrt(emitx*ring.betax)*cos(phix) + ring.etax*p.de;
		p.xl = -sqrt(emitx/ring.betax)*(ring.alphax*cos(phix) + sin(phix)) + ring.etaxl*p.de;
	}
  #endif
	bun.is_sorted = false;
}

static double get_residue(const my_Dvector& spos, const my_Dvector& distr, const my_Dvector& distr_old)
{
	double res = 0.0;
	double ds = spos[1] - spos[0];
	for (int i=0;i<spos.size();++i){res += (distr[i] - distr_old[i]) * (distr[i] - distr_old[i]) * ds;}
	return res;
}
static double get_residue(const my_Dvector& spos, const my_Dvector& distr)
{
	double res = 0.0;
	double ds = spos[1] - spos[0];
	for (int i=0;i<spos.size();++i){res += distr[i] * distr[i] * ds;}
	return res;
}

void generate_bunch(const Ring_t& ring, const Wake_t& wake, Bunch_t& bun, Results_t& results)
{

    auto& cav_s = ring.cav.ref_to_xi();
    auto& cav_V = ring.cav.ref_to_yi();

	// get the wake function at the cavity longitudinal points
	my_Dvector&& wakeF = wake.Wl.get_wake_at_points(cav_s, bun.Ib * ring.rev_time);

	// Now we iterate to get the equilibrium distribution
	bool converged (false);
	double espread (ring.espread);
	do {
		double scale (ring.mom_comp * ring.circum * ring.energy);
		my_Dvector&& distr_old (get_distribution_from_potential(cav_s,cav_V,espread,scale));
		double&& residue_old (get_residue(cav_s,distr_old));
		int&& count (0);
		my_Dvector distr (distr_old);
		while (!converged){
			my_Dvector&& V (cav_V + convolution(distr, wakeF,'same'));

			distr = get_distribution_from_potential(ring,espread,cav_s,V);
			double&& residue (get_residue(cav_s,distr,distr_old));

			if (residue < residue_old){if (residue < 1e-20) {converged = true;}}
			else {if (++count > 10){espread *= 1.01; break;}}

			swap(distr_old, distr);
			swap(residue_old,residue);
		}
	} while (!converged);


	default_random_engine gen(19880419);
	normal_distribution<double> espread_dist(0.0,ring.espread);
	exponential_distribution<double> emitx_dist(1/(2*ring.emitx));
	uniform_real_distribution<double> phix_dist(0.0,TWOPI);
    uniform_real_distribution<double> ss_dist(idistr.front(),idistr.back());
    Interpola_t idistr_interpol (idistr, s_distr);
  #ifdef OPENMP
    fprintf(stdout,"here\n");
    my_PartVector& par = bun.particles;
    #pragma omp parallel for schedule(guided,1)
    for (int i=0;i<par.size();++i){
        double&& emitx = emitx_dist(gen);
        double&& phix  = phix_dist(gen);
        par[i].de = espread_dist(gen);
        par[i].ss = idistr_interpol.get_y(ss_dist(gen));
        par[i].xx =  sqrt(emitx*ring.betax)*cos(phix) + ring.etax*par[i].de;
        par[i].xl = -sqrt(emitx/ring.betax)*(ring.alphax*cos(phix) + sin(phix)) + ring.etaxl*par[i].de;
    }
  #else
	for (auto& p:bun.particles){
		double&& emitx = emitx_dist(gen);
		double&& phix  = phix_dist(gen);
		p.de = espread_dist(gen);
		p.ss = idistr_interpol.get_y(ss_dist(gen));
		p.xx =  sqrt(emitx*ring.betax)*cos(phix) + ring.etax*p.de;
		p.xl = -sqrt(emitx/ring.betax)*(ring.alphax*cos(phix) + sin(phix)) + ring.etaxl*p.de;
	}
  #endif
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
        // bun.general_sort();

        results.register_Wkicks(n, wake.apply_kicks(bun,kick_stren, ring.betax));
        results.register_FBkick(n,   fb.apply_kick(bun, xx_ave,     ring.betax));
    }
    results.calc_stats(bun,results.get_nturns());
}
