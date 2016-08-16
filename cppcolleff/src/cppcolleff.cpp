
#include <cppcolleff/cppcolleff.h>

void generate_bunch(const Ring_t& ring, Bunch_t& bun)
{
	default_random_engine generator(19880419);
  	normal_distribution<double> espread_dist(0.0,ring.espread);
	exponential_distribution<double> emitx_dist(1/ring.emitx);
	uniform_real_distribution<double> phix_dist(0.0,TWOPI);

    auto& cav_s = ring.cav_s;
    auto& cav_V = ring.cav_V;
    my_Dvector ipot   (cav_s.size(),0.0);
    my_Dvector idistr, s_distr;

    double min_ipot = ipot[0];
    for (int i=1;i<cav_s.size();++i){
        ipot[i] = ipot [i-1] + (cav_V[i] + cav_V[i-1])*(cav_s[i]-cav_s[i-1])/2;
        if (ipot[i] < min_ipot) min_ipot = ipot[i];
    }

    double exp_fact = 1/(ring.mom_comp*ring.circum*ring.E)/(ring.espread*ring.espread);
    idistr.push_back(0.0);
    s_distr.push_back(cav_s[0]);
    for (int i=1; i<cav_s.size();++i){
        double distr = (exp(-exp_fact*(ipot[i]  -min_ipot))
                      - exp(-exp_fact*(ipot[i-1]-min_ipot)))*(cav_s[i]-cav_s[i-1])/2;
        if (distr !=0.0) {
            idistr.push_back(idistr.back() + distr);
            s_distr.push_back(cav_s[i]);
        }
    }

    uniform_real_distribution<double> ss_dist(idistr.front(),idistr.back());
	for (auto& p:bun.particles){
		double emitx = emitx_dist(generator);
		double phix  = phix_dist(generator);
		p.de = espread_dist(generator);
		p.ss = interpola2(ss_dist(gen),idistr,s_distr);
		p.xx = sqrt(emitx*ring.betax)*cos(phix) + ring.etax*p.de;
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

    fprintf(stdout,"%7s %12s  %12s %12s  %12s \n",
            "turn","<de> [%]","std(de) [%]","<xx> [mm]","std(xx) [mm]");
    for (long n=0;n<results.get_nturns();n++){
        double xx_ave = results.calc_stats(bun,n);
        /* convention: positive ss, means particle behind the sinchronous particle;
         First do single particle tracking:*/
        ring.track_one_turn(bun);

        // After this sorting, the particles will be ordered from head to tail.
        // It means, from smaller ss to bigger ss.
        bun.general_sort();

        if ((n % 10) == 0 || n == 0) {
            fprintf(stdout,"%07lu %12.6f  %12.6f %12.6f  %12.6f \n",n,
                    1e2*results.ave[n].de,1e2*results.std[n].de,
                    1e3*results.ave[n].xx,1e3*results.std[n].xx);
        }

        results.register_Wkicks(n, wake.apply_kicks(bun,kick_stren, ring.betax));
        results.register_FBkick(n,   fb.apply_kick(bun, xx_ave,     ring.betax));
    }
    results.calc_stats(bun,results.get_nturns());
}
