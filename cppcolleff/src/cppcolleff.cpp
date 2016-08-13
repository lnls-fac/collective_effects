
#include <cppcolleff/cppcolleff.h>

void generate_bunch(const Ring_t& ring, Bunch_t& bun)
{
	default_random_engine generator(19880419);
  	normal_distribution<double> norm(0.0,1.0);
	exponential_distribution<double> expo(ring.emitx);
	uniform_distribution<double> unif(0.0,TWOPI);

	for (auto& p:particles){
		double emitx = expo(generator);
		double phix  = unif(generator);
		p.de = norm(generator)*ring.espread;
		p.ss = norm(generator)*ring.bunlen;
		p.xx = sqrt(emitx*ring.betax)*cos(phix) + ring.etax*p.de;
		p.xl = -sqrt(emitx/ring.betax)*(ring.alphax*cos(phix) + sin(phix)) + ring.etaxl*p.de;
	}
	is_sorted = false;
}


void do_tracking(
    const Ring_t& ring,
    const Wake_t& wake,
    const Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results)
{
    //current dependent strength of the kick:
    const double kick_stren = ring.T0 / ring.energy * bun.Ib / bun.num_part;

    for (long n=0;n<results.get_nturns();n++){
        double xx_ave = results.calc_stats(n,bun);
        /* convention: positive ss, means particle behind the sinchronous particle;
         First do single particle tracking:*/
        ring.track_one_turn(bun);

        // After this sorting, the particles will be ordered from head to tail.
        // It means, from smaller ss to bigger ss.
        bun.insertion_sort();

        if ((n % 10) == 0 || n == 0) {
            fprintf(stdout,"%12.6f  %12.6f \n",1e3*results.de_ave[n],1e3*results.de_std[n]);
        }

        results.set_Wkicks(n, wake.apply_kicks(bun,kick_stren));
        results.set_FBkick(n, fb.apply_kick(xx_ave,bun));
    }
    results.calc_stats(results.get_nturns(),bun);
}
