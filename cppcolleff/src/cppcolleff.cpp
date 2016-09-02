#include <cppcolleff/cppcolleff.h>

static void _generate_bunch_thread(
	const Ring_t& ring,
	my_PartVector& p,
	const unsigned int seed,
	const Interpola_t& idistr_interpol,
	const int init,
	const int final)
{
	const my_Dvector& idist = idistr_interpol.ref_to_xi();
	normal_distribution<double> espread_dist(0.0,ring.espread);
	exponential_distribution<double> emitx_dist(1/(2*ring.emitx));
	uniform_real_distribution<double> phix_dist(0.0,TWOPI);
	uniform_real_distribution<double> ss_dist(idist.front(),idist.back());
	default_random_engine gen(seed);

	for (int i=init;i<final;++i){
		double&& emitx = emitx_dist(gen);
		double&& phix  = phix_dist(gen);
		double&& Ax    = sqrt(emitx/ring.betax);
		double&& Acosx = Ax*cos(phix);
		double&& Asinx = Ax*sin(phix);
		p[i].de = espread_dist(gen);
		p[i].ss = idistr_interpol.get_y(ss_dist(gen));
		p[i].xx =  Acosx*ring.betax          + ring.etax *p[i].de;
		p[i].xl = -Acosx*ring.alphax - Asinx + ring.etaxl*p[i].de;
	}
}
static void _generate_bunch(const Ring_t& ring, Bunch_t& bun, unsigned int seed)
{
	Interpola_t idistr_interpol (ring.get_integrated_distribution());
	my_PartVector& p = bun.particles;

	int nr_th = NumThreads::get_num_threads();
	my_Ivector lims (bounds_for_threads(nr_th,0,p.size()));
	vector<thread> ths;

	for (int i=0;i<nr_th-1;++i){
		ths.push_back(thread(
			_generate_bunch_thread, ref(ring),ref(p),seed+i,
									ref(idistr_interpol), lims[i], lims[i+1]
		));
	}
	_generate_bunch_thread(ring, p, seed+nr_th-1, idistr_interpol,
						   lims[nr_th-1], lims[nr_th]);
	for(auto& th:ths) th.join();

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
static my_Dvector _const_espread_iteration_Haissinski(const Ring_t ring, const my_Dvector KickF)
{
	my_Dvector distr_old (ring.get_distribution());
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

		if (residue < residue_old) if (residue < 1e-6) return V;
		else if (++count > 10) return my_Dvector ();

		distr_old.swap(distr);
		swap(residue_old,residue);
	}
}

my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib)
{
	// get the wake function at the cavity longitudinal points (actually it is the kick)
	auto& cav_s = ring.cav.ref_to_xi();
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.T0 / ring.energy);

	return my_Dvector (_const_espread_iteration_Haissinski(ring, KickF));
}

double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib)
{
	// get the wake function at the cavity longitudinal points (actually it is the kick)
	auto& cav_s = ring.cav.ref_to_xi();
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.T0 / ring.energy);

	my_Dvector V (_const_espread_iteration_Haissinski(ring, KickF));

	// Now use bissection to get the equilibrium distribution if needed
	double final_espread (ring.espread);
	if ( V.empty() ) {
		double init_espread (ring.espread);
		double delta_spread (ring.espread); //begin with a delta equal to the natural energy spread
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
			V = _const_espread_iteration_Haissinski(ring, KickF);
		}
		ring.espread = init_espread;
	}
	return final_espread;
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
