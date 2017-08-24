#include <cppcolleff/cppcolleff.h>

static ThreadVars ThreadInfo(false);

static void _generate_bunch_thread(
	const Ring_t& ring,
	my_PartVector& p,
	const unsigned int th,
	const Interpola_t& idistr_interpol,
	const int init,
	const int final)
{
		const my_Dvector& idist = idistr_interpol.ref_to_xi();
		normal_distribution<double> espread_dist(0.0,ring.espread);
		exponential_distribution<double> emitx_dist(1/(2*ring.emitx));
		uniform_real_distribution<double> phix_dist(0.0,TWOPI);
		uniform_real_distribution<double> ss_dist(idist.front(),idist.back());

  	for (int i=init;i<final;++i){
	  		double&& emitx = emitx_dist(ThreadInfo.gens[th]);
	  		double&& phix  = phix_dist(ThreadInfo.gens[th]);
	  		double&& Ax    = sqrt(emitx/ring.betax);
	  		double&& Acosx = Ax*cos(phix);
	  		double&& Asinx = Ax*sin(phix);
	  		p[i].de = espread_dist(ThreadInfo.gens[th]);
	  		p[i].ss = idistr_interpol.get_y(ss_dist(ThreadInfo.gens[th]));
	  		p[i].xx =  Acosx*ring.betax          + ring.etax *p[i].de;
	  		p[i].xl = -Acosx*ring.alphax - Asinx + ring.etaxl*p[i].de;
  	}
}
void generate_bunch(const Ring_t& ring, Bunch_t& bun)
{
	Interpola_t idistr_interpol (ring.get_integrated_distribution());
	my_PartVector& p = bun.particles;

	int nr_th = ThreadInfo.get_num_threads();
	my_Ivector lims (ThreadInfo.get_bounds(0,p.size()));
	vector<thread> ths;

	for (unsigned int i=0;i<nr_th-1;++i){
		  ths.push_back(thread(_generate_bunch_thread, ref(ring),ref(p), i,
							   ref(idistr_interpol), lims[i], lims[i+1]));
	}
	_generate_bunch_thread(ring, p, nr_th-1, idistr_interpol,
									 lims[nr_th-1], lims[nr_th]);
	for(auto& th:ths) th.join();

	bun.is_sorted = false;
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
static my_Dvector _const_espread_iteration_Haissinski(
	const Ring_t ring,
	const my_Dvector KickF,
	const int niter,
	const my_Dvector distr_ini)
{
	my_Dvector distr_old (distr_ini);
	auto& cav_s = ring.cav.ref_to_xi();
	auto& cav_V = ring.cav.ref_to_yi();
	double residue_old (get_residue(cav_s,distr_old));
	int count (0);
	double&& ds = (cav_s[1]-cav_s[0]);
	while (true) {
		my_Dvector V (convolution_same(KickF,distr_old)); // variable to be returned;

	  #ifdef OPENMP
		#pragma omp parallel for schedule(guided,1)
	  #endif
		for (int i=0;i<V.size();++i){
			V[i] *= ds; // to fix the scale of the convolution.
			V[i] += cav_V[i];
		}

		my_Dvector&& distr = ring.get_distribution(V);
		double&& residue = get_residue(cav_s,distr,distr_old);

		// fprintf(stdout,"%3d: %8.3g\n",count,residue);
		if (residue < residue_old){if (residue < 1e-6) return V;}
		else {if (++count > niter) return my_Dvector ();};

		distr_old.swap(distr);
		swap(residue_old,residue);
	}
}

static my_Dvector _const_espread_iteration_Haissinski(
	const Ring_t ring,
	const my_Dvector KickF)
{
	int niter = 100;
	my_Dvector distr_ini (ring.get_distribution());
	return _const_espread_iteration_Haissinski(ring, KickF, niter, distr_ini);
}
static my_Dvector _const_espread_iteration_Haissinski(
	const Ring_t ring,
	const my_Dvector KickF,
	const int niter)
{
	my_Dvector distr_ini (ring.get_distribution());
	return _const_espread_iteration_Haissinski(ring, KickF, niter, distr_ini);
}
static my_Dvector _const_espread_iteration_Haissinski(
	const Ring_t ring,
	const my_Dvector KickF,
	const my_Dvector distr_ini)
{
	int niter = 100;
	return _const_espread_iteration_Haissinski(ring, KickF, niter, distr_ini);
}


my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib,
	const int niter,
	const my_Dvector distr_ini)
{
	// get the wake function at the cavity longitudinal points (actually it is the kick)
	auto& cav_s = ring.cav.ref_to_xi();
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.T0 / ring.energy);

	return _const_espread_iteration_Haissinski(ring, KickF, niter, distr_ini);
}
my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib)
{
	int niter (100);
	my_Dvector distr_ini (ring.get_distribution());
	return solve_Haissinski_get_potential(wake, ring, Ib, niter, distr_ini);
}
my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib,
	const int niter)
{
	my_Dvector distr_ini (ring.get_distribution());
	return solve_Haissinski_get_potential(wake, ring, Ib, niter, distr_ini);
}
my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib,
	const my_Dvector distr_ini)
{
	int niter (100);
	return solve_Haissinski_get_potential(wake, ring, Ib, niter, distr_ini);
}

double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib,
	const int niter,
	const my_Dvector distr_ini)
{
	// get the wake function at the cavity longitudinal points (actually it is the kick)
	auto& cav_s = ring.cav.ref_to_xi();
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.T0 / ring.energy);

	my_Dvector V (_const_espread_iteration_Haissinski(ring, KickF, niter, distr_ini));

	// Now use bissection to get the equilibrium distribution if needed
	double final_espread (ring.espread);
	if ( V.empty() ) {
		double init_espread (ring.espread);
		// double delta_spread (ring.espread); //begin with a delta equal to the natural energy spread
		double delta_spread (1e-4);
		bool conv_once (false);

		while (delta_spread > 1e-5) {
			if (V.empty()) {
				if (conv_once) delta_spread /= 2;
				ring.espread += delta_spread;
			}
			else {
				conv_once = true;
				final_espread = ring.espread;
				delta_spread /= 2;
				ring.espread -= delta_spread;
			}
			// fprintf(stdout,"ok\n");
			V = _const_espread_iteration_Haissinski(ring, KickF, niter, distr_ini);
			// fprintf(stdout,"ok2\n");
		}
		ring.espread = init_espread;
	}
	return final_espread;
}

double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib)
{
	int niter = 100;
	my_Dvector distr_ini (ring.get_distribution());
	return find_equilibrium_energy_spread(wake, ring, Ib, niter, distr_ini);
}
double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib,
	const int niter)
{
	my_Dvector distr_ini (ring.get_distribution());
	return find_equilibrium_energy_spread(wake, ring, Ib, niter, distr_ini);
}
double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib,
	const my_Dvector distr_ini)
{
	int niter = 100;
	return find_equilibrium_energy_spread(wake, ring, Ib, niter, distr_ini);
}

void single_bunch_tracking(
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

        results.register_Wkicks(n, wake.apply_kicks(bun,kick_stren, ring.betax));
        results.register_FBkick(n,   fb.apply_kick(bun, xx_ave,     ring.betax));
    }
    results.calc_stats(bun,results.get_nturns());
}


void multi_bunch_tracking(
	const Ring_t& ring,
	const Wake_t& long_range_wake,
	const Wake_t& short_range_wake,
	Feedback_t& fb,
	vector<Bunch_t>& buns,
	vector<Results_t>& results)
{

}
