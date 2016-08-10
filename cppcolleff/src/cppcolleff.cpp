
#include <cppcolleff/cppcolleff.h>

void do_tracking(
    Ring_t& ring,
    Wake_t& wake,
    Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results)
{

    //current dependent strength of the kick:
    const double kick_stren = ring.T0 / ring.energy * bun.Ib / bun.num_part;

    for (int n=0;n<results.nturns;n++){
        results.calc_stats(n,bun);
        /* convention: positive ss, means particle behind the sinchronous particle;
         First do single particle tracking:*/
        ring.track_one_turn(bun);

        // After this sorting, the particles will be ordered from head to tail.
        // It means, from smaller ss to bigger ss.
        bun.InsertionSort();

        if ((n % 10) == 0 || n == 0) {
            fprintf(stdout,"%12.6f  %12.6f \n",1e3*results.de_ave[n],1e3*results.de_std[n]);
        }

        results.set_Wkicks(n, wake.apply_kicks(bun,kick_stren));

        results.set_FBkick(n, fb.apply_kick(results.xx_ave,n,bun));
    }
}
