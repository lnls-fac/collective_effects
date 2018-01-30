#include <cppcolleff/cppcolleff.h>


static double get_residue(
	const my_Dvector& spos,
	const my_Dvector& distr,
	const my_Dvector& distr_old,
	ThreadPool& pool)
{
	double res (0.0);
	double ds (spos[1] - spos[0]);

	unsigned int nr_th = get_num_threads();
	my_Ivector lims (get_bounds(0, spos.size()));
	vector< std::future<double> > results;
	auto fun = [&](int ini, int fin){
		double r (0.0);
		for (int ii=ini;ii<fin;++ii){r += (distr[ii] - distr_old[ii]) * (distr[ii] - distr_old[ii]) * ds;}
		return r;};
	for (unsigned int i=0;i<nr_th;++i){
		results.emplace_back(pool.enqueue(fun, lims[i], lims[i+1]));
	}
	for(auto&& result: results) res += result.get();

	return res;
}
static my_Dvector _const_espread_iteration_Haissinski(
	const Ring_t ring,
	const my_Dvector KickF,
	const bool req_convergence,
	const int niter,
	const double weight,
	const my_Dvector distr_ini,
	ThreadPool& pool,
	Convolve_t& conv)
{
	my_Dvector distr_old (distr_ini);
	auto& cav_s = ring.cav.ref_to_xi();
	auto& cav_V = ring.cav.ref_to_yi();
	double residue_old (get_residue(cav_s, distr_old, my_Dvector(cav_s.size(),0.0), pool));
	int count (0);
	double&& ds = (cav_s[1]-cav_s[0]);

	unsigned int nr_th = get_num_threads();
	my_Ivector lims (get_bounds(0, cav_V.size()));
	vector< std::future<int> > res;
	while (true) {
		conv.prepare(KickF, distr_old);
		my_Dvector V (conv.execute_same()); // variable to be returned;
		// my_Dvector V (convolution_same(KickF, distr_old, pool)); // variable to be returned;
		// correct the scale of the convolution and add cavity to potential:
	    auto fun = [&](my_Dvector& V1, int ini, int fin)
	    {for (int i=ini;i<fin;++i){V1[i] *= ds; V1[i] += cav_V[i];}  return 1;};
	    for (unsigned int i=0;i<nr_th;++i){
	        res.emplace_back(pool.enqueue(fun, ref(V), lims[i], lims[i+1]));
	    }
	    for(auto&& result: res) result.get();
		res.clear();

		my_Dvector&& distr = ring.get_distribution(pool, V, cav_s);
		double&& residue = get_residue(cav_s, distr, distr_old, pool);

		if (residue < residue_old){if (residue < 1e-6) return V;}
		else {
			if (++count > niter) {
				if (req_convergence) return my_Dvector ();
				else return V;
			}
		}

		for (unsigned int i=0; i<distr.size(); ++i){
			distr[i] += weight*distr_old[i];
			distr[i] /= weight+1;
		}

		distr_old.swap(distr);
		swap(residue_old,residue);
	}
}


my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib,
	const bool req_convergence,
	const int niter,
	const double weight,
	const my_Dvector distr_ini)
{
	ThreadPool pool (get_num_threads());
	auto& cav_s = ring.cav.ref_to_xi();
	Convolve_t conv(cav_s.size(), cav_s.size());
	return solve_Haissinski_get_potential(wake, ring, Ib, pool, conv, req_convergence, niter, weight, distr_ini);
}

my_Dvector solve_Haissinski_get_potential(
	const Wake_t& wake,
	const Ring_t& ring,
	const double& Ib,
	ThreadPool& pool,
	Convolve_t& conv,
	const bool req_convergence,
	const int niter,
	const double weight,
	const my_Dvector distr_ini)
{
	// get the wake function at the cavity longitudinal points (actually it is the kick)
	my_Dvector ini_distr (distr_ini);
	if (ini_distr.size() == 0) {ini_distr = ring.get_distribution();}

	auto& cav_s = ring.cav.ref_to_xi();
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.circum / light_speed / ring.energy);
	return _const_espread_iteration_Haissinski(ring, KickF, req_convergence, niter, weight, ini_distr, pool, conv);
}


double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib,
	const int niter,
	const double weight,
	const my_Dvector distr_ini)
{
	ThreadPool pool (get_num_threads());
	auto& cav_s = ring.cav.ref_to_xi();
	Convolve_t conv(cav_s.size(), cav_s.size(), true);
	return find_equilibrium_energy_spread(wake, ring, Ib, pool, conv, niter, weight, distr_ini);
}

double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring,
	const double& Ib,
	ThreadPool& pool,
	Convolve_t& conv,
	const int niter,
	const double weight,
	const my_Dvector distr_ini)
{
	// get the wake function at the cavity longitudinal points (actually it is the kick)
	auto& cav_s = ring.cav.ref_to_xi();
	my_Dvector&& KickF = wake.Wl.get_wake_at_points(cav_s, -Ib * ring.circum / light_speed / ring.energy);

	my_Dvector V (_const_espread_iteration_Haissinski(ring, KickF, true, niter, weight, distr_ini, pool, conv));

	// Now use bissection to get the equilibrium distribution if needed
	double final_espread (ring.espread);
	if ( V.empty() ) {
		double init_espread (ring.espread);
		// double delta_spread (ring.espread); //begin with a delta equal to the natural energy spread
		double delta_spread (1e-4);
		bool conv_once (false);

		while (delta_spread > 1e-7) {
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
			V = _const_espread_iteration_Haissinski(ring, KickF, true, niter, weight, distr_ini, pool, conv);
			// fprintf(stdout,"ok2\n");
		}
		ring.espread = init_espread;
	}
	return final_espread;
}

my_Dvector long_simul_with_haissinki(
	const Wake_t& wake,
	Ring_t& ring,
	my_Dvector& currs,
    my_Dvector& init_espread,
	int niter,
	double weight)
{
	const my_Dvector& ss = ring.cav.ref_to_xi();
    my_Dvector espread (currs.size(), 0.0);
    double espread0 (ring.espread);
	ThreadPool pool (get_num_threads());
	Convolve_t conv(ss.size(), ss.size(), true);

	ofstream fp("results.txt");
	if (fp.fail()) exit(1);
	fp.setf(fp.left | fp.showpos | fp.scientific);
	fp << setw(30) << "% number_of_iterations" << niter << endl;
	fp << setw(22) << "Ib [A]";
	fp << setw(22) << "spread";
	fp << setw(22) << "bun len [m]";
	fp << setw(22) << "sync ph [m]" << endl;

	cout.setf(fp.left);
	cout << setw(22) << "Ib [mA]";
	cout << setw(22) << "spread x 1000";
	cout << setw(22) << "bun len [mm]";
	cout << setw(22) << "sync ph [mm]" << endl;

    my_Dvector&& dist = ring.get_distribution();
    for (int i=0; i<currs.size(); ++i){
        ring.espread = max(init_espread[i], ring.espread);
        espread[i] = find_equilibrium_energy_spread(wake, ring, currs[i],
                                            pool, conv, niter, weight, dist);
        ring.espread = espread[i];
		my_Dvector&& V = solve_Haissinski_get_potential(wake, ring, currs[i],
									pool, conv, true, niter*2, weight, dist);
        dist = ring.get_distribution(V);

		//trapezoidal integration to find synch pos and bunch length
		double&& s0 = dist[0]*ss[0];  double&& s02 = dist[0]*ss[0]*ss[0];
		s0 += dist.back()*ss.back();  s02 += dist.back()*ss.back()*ss.back();
		s0 /= 2;					  s02 /= 2;
		for (auto&& ii=1; ii<ss.size()-1; ++ii){
			double&& dummy = dist[ii]*ss[ii];
			// phase shift:
        	 s0 += dummy;
			 dummy *= ss[ii];
			 // bunch length
			 s02 += dummy;
		}
		double&& ds = ss[1]-ss[0];   s0 *= ds;   s02 *= ds;
		double&& bl = sqrt(s02 - s0*s0);

		fp << setw(22) << currs[i];
	 	fp << setw(22) << espread[i];
	 	fp << setw(22) << bl;
	 	fp << setw(22) << s0 << endl;

	 	cout << setw(22) << currs[i]*1e3;
	 	cout << setw(22) << espread[i]*1e3;
	 	cout << setw(22) << bl*1e3;
	 	cout << setw(22) << s0*1e3 << endl;

		char fname[50];
		sprintf(fname, "curr_%05.3fmA_distr_ss.txt", currs[i]*1e3);
        save_distribution_to_file(fname, dist, ss.front(), ss.back(), ss.size());
		sprintf(fname, "curr_%05.3fmA_potential.txt", currs[i]*1e3);
        save_distribution_to_file(fname, V, ss.front(), ss.back(), ss.size());
	}
    fp.close();
    ring.espread = espread0;
    return espread;
}


void single_bunch_tracking(
    const Ring_t& ring,
    Wake_t& wake,
    Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results)
{
    //current dependent strength of the kick:
    const double kick_stren = ring.circum / light_speed / ring.energy * bun.Ib / bun.num_part;

	ThreadPool pool (get_num_threads());

    for (long n=0;n<results.get_nturns();n++){
        double&& xx_ave = results.calc_stats(bun, n, pool);
        /* convention: positive ss means particle behind the sinchronous particle;
         First do single particle tracking:*/
        ring.track_one_turn(bun, n, pool);

        results.register_Wkicks(n, wake.apply_kicks(bun, kick_stren, ring.betax, pool));
        results.register_FBkick(n,   fb.apply_kick(bun, xx_ave,     ring.betax));
    }
    results.calc_stats(bun, results.get_nturns(), pool);
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
