#include <cppcolleff/cppcolleff.h>

static void _generate_bunch(const Ring_t& ring, Bunch_t& bun, unsigned int seed)
{
	Interpola_t&& idistr_interpol = ring.get_integrated_distribution();
	const my_Dvector& idist = idistr_interpol.ref_to_xi();

	normal_distribution<double> espread_dist(0.0,ring.espread);
	exponential_distribution<double> emitx_dist(1/(2*ring.emitx));
	uniform_real_distribution<double> phix_dist(0.0,TWOPI);
	uniform_real_distribution<double> ss_dist(idist.front(),idist.back());
	default_random_engine gen(seed);
  #ifdef OPENMP
	my_PartVector& par = bun.particles;
	#pragma omp ṕarallel for schedule(guided,1)
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
void generate_bunch(const Ring_t& ring, Bunch_t& bun)
{
	unsigned int seed(19880419);
	_generate_bunch(ring,bun,seed);
}
void generate_bunch(const Ring_t& ring, Bunch_t& bun, unsigned int seed)
{
	_generate_bunch(ring,bun,seed);
}


static double get_residue(const my_Dvector& spos, const my_Dvector& distr, const my_Dvector& distr_old)
{
	double res = 0.0;
	double ds = spos[1] - spos[0];
	#ifdef OPENMP
	  #pragma omp parallel for schedule(guided,1) reduction(+:res)
	#endif
	for (int i=0;i<spos.size();++i){res += (distr[i] - distr_old[i]) * (distr[i] - distr_old[i]) * ds;}
	return res;
}
static double get_residue(const my_Dvector& spos, const my_Dvector& distr)
{
	double res = 0.0;
	double ds = spos[1] - spos[0];
	#ifdef OPENMP
	  #pragma omp parallel for schedule(guided,1) reduction(+:res)
	#endif
	for (int i=0;i<spos.size();++i){res += distr[i] * distr[i] * ds;}
	return res;
}
my_Dvector solve_Haissinski(const Wake_t& wake, Ring_t& ring, const double& Ib)
{
    auto& cav_s = ring.cav.ref_to_xi();
    auto& cav_V = ring.cav.ref_to_yi();

	my_Dvector V (cav_V.size(),0.0); // variable to be returned;

	// get the wake function at the cavity longitudinal points (actually it is the kick)
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.T0);

	// Now we iterate to get the equilibrium distribution
	bool converged (false);
	do {
		my_Dvector distr_old (ring.get_distribution());
		double residue_old (get_residue(cav_s,distr_old));
		int count (0);
		double&& ds = (cav_s[1]-cav_s[0]);
		fprintf(stdout,"%3i  %g\n",0,residue_old);
		while (!converged){
			V = move(convolution_same(KickF,distr_old));
		  #ifdef OPENMP
		  	#pragma omp parallel for schedule(guided,1)
		  #endif
			for (int i=0;i<V.size();++i){
				V[i] *= ds; // to fix the scale of the convolution.
				V[i] += cav_V[i];
			}

			my_Dvector&& distr = ring.get_distribution(V);
			double&& residue = get_residue(cav_s,distr,distr_old);
			fprintf(stdout,"%3i  %g\n",count,residue);
			if (residue < residue_old){
				if (residue < 1e-6) {converged = true;}
			}
			else {
				if (++count > 10){
					ring.espread *= 1.02;
					fprintf(stdout,"espread = %f\n",ring.espread);
					break;
				}
			}

			swap(distr_old, distr);
			swap(residue_old,residue);
		}
	} while (!converged);

	return V;
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
