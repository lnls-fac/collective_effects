#include <cppcolleff/Results.h>

double calc_moments(
    const my_PartVector& p,
    Particle_t& ave,
    Particle_t& std,
    const int init,
    const int fin,
    const bool this_turn)
{
    if (!this_turn) {
        double avexx (0.0);
        for (int i=init;i<fin;++i) avexx += p[i].xx;
        return avexx;
    }
    for (int i=init;i<fin;++i){
        ave += p[i];
        std += p[i]*p[i];
    }
    return ave.xx;
}

double Results_t::calc_stats(
    const Bunch_t& bun,
    const long turn)
{
    ThreadPool pool (get_num_threads());
    return calc_stats(bun, turn, pool);
}

double Results_t::calc_stats(
    const Bunch_t& bun,
    const long turn,
    ThreadPool& pool)
{
    const bool this_turn = calc_this_turn(turn);
    ave.push_back(Particle_t (0.0));
    std.push_back(Particle_t (0.0));
    auto& rave = ave.back();
    auto& rstd = std.back();
    //
    const my_PartVector& p = bun.particles;

    unsigned int nr_th = get_num_threads();
    my_Ivector lims (get_bounds(0,p.size()));
    my_PartVector aves (nr_th,Particle_t());
    my_PartVector stds (nr_th,Particle_t());
    std::vector< std::future<double> > res;

    for (unsigned int i=0;i<nr_th;++i){
        res.emplace_back(pool.enqueue(calc_moments, ref(p), ref(aves[i]),
                         ref(stds[i]), lims[i], lims[i+1], this_turn));
    }
    if (!this_turn){
        double avexx (0.0);
        for(int i=0; i<nr_th; ++i) avexx += res[i].get();
        return avexx /= bun.num_part;
    }
    for(int i=0; i<nr_th; ++i){
        res[i].get();
        rave += aves[i];
        rstd += stds[i];
    }
    rave /= bun.num_part;
    rstd /= bun.num_part;
    rstd = sqrt(rstd - rave * rave);

    if (dump_bunch_to_file && dump_this_turn(turn)) {
        char filename[50];
        sprintf(filename,"turn%07lu.txt",turn);
        bun.write_bunch_to_file(filename);
    }
    if (print_in_screen) {
        if (turn == 0) {
            fprintf(stdout,"%7s %12s %12s   %12s %12s   %12s %12s   %12s %12s \n","turn",
                    "<xx> [um]","std(xx) [um]","<xl> [urad]","std(xl) [urad]",
                    "<de> [%]","std(de) [%]","<ss> [mm]","std(ss) [mm]"
                );
        }
        if (print_this_turn(turn)){
            fprintf(stdout,"%07lu %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f \n",turn,
                    1e6*ave.back().xx,1e6*std.back().xx,
                    1e6*ave.back().xl,1e6*std.back().xl,
                    1e2*ave.back().de,1e2*std.back().de,
                    1e3*ave.back().ss,1e3*std.back().ss
                );
        }
    }
    return rave.xx;
}

void Results_t::register_FBkick(const long turn, const double& kik)
{
    if (FB && calc_this_turn(turn)) {FBkick.push_back(kik);}
}

void Results_t::register_Wkicks(const long turn, const my_Dvector& kik)
{
    if (Wl && calc_this_turn(turn)) {Wlkick.push_back(kik[0]);}
    if (Wd && calc_this_turn(turn)) {Wdkick.push_back(kik[1]);}
    if (Wq && calc_this_turn(turn)) {Wqkick.push_back(kik[2]);}
}
