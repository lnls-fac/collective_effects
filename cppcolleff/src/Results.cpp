#include <cppcolleff/Results.h>

static ThreadVars ThreadInfo(false);

void Results_t::write_bunch_to_file(const Bunch_t& bun, const char* filename) const
{
    FILE* fp = fopen(filename,"w");
    fprintf(fp,"#%20s %20s %20s %20s\n","xx [m]", "xl [m]", "de", "ss [m]");

    for (auto& p:bun.particles){
        fprintf(fp,"%20.7g  %20.7g %20.7g %20.7g\n",p.xx,p.xl,p.de,p.ss);
    }
    fclose(fp);
}

int calc_moments(
    const my_PartVector& p,
    my_Dvector& moms,
    const int init,
    const int fin,
    const bool this_turn)
{
    for (int i=init;i<fin;++i){
        moms[0] += p[i].xx;
        if (!this_turn) continue;
        moms[1] += p[i].xl;
        moms[2] += p[i].de;
        moms[3] += p[i].ss;
        moms[4] += p[i].xx*p[i].xx;
        moms[5] += p[i].xl*p[i].xl;
        moms[6] += p[i].de*p[i].de;
        moms[7] += p[i].ss*p[i].ss;
    }
    return 1;
}


double Results_t::calc_stats(
    const Bunch_t& bun,
    ThreadPool& pool,
    const long turn)
{
    const bool this_turn = calc_this_turn(turn);
    ave.push_back(Particle_t (0.0));
    std.push_back(Particle_t (0.0));
    auto& rave = ave.back();
    auto& rstd = std.back();
    //
    const my_PartVector& p = bun.particles;

    unsigned int nr_th = ThreadInfo.get_num_threads();
    my_Ivector lims (ThreadInfo.get_bounds(0,p.size()));
    vector<vector<double>> moms (nr_th,vector<double>(8,0.0));
    std::vector< std::future<int> > res;

    for (unsigned int i=0;i<nr_th;++i){
        res.emplace_back(pool.enqueue(calc_moments, ref(p), ref(moms[i]),
                         lims[i], lims[i+1], this_turn));
    }

    for(int i=0; i<nr_th; ++i){
        res[i].get();
        rave.xx += moms[i][0];
        if (!this_turn) continue;
        rave.xl += moms[i][1];
        rave.de += moms[i][2];
        rave.ss += moms[i][3];
        rstd.xx += moms[i][4];
        rstd.xl += moms[i][5];
        rstd.de += moms[i][6];
        rstd.ss += moms[i][7];
    }
    rave.xx /= bun.num_part;
    if (!this_turn) return rave.xx;
    rave.xl /= bun.num_part;
    rave.de /= bun.num_part;
    rave.ss /= bun.num_part;
    rstd.xx /= bun.num_part;
    rstd.xl /= bun.num_part;
    rstd.de /= bun.num_part;
    rstd.ss /= bun.num_part;
    rstd.xx = sqrt(rstd.xx - rave.xx * rave.xx);
    rstd.xl = sqrt(rstd.xl - rave.xl * rave.xl);
    rstd.de = sqrt(rstd.de - rave.de * rave.de);
    rstd.ss = sqrt(rstd.ss - rave.ss * rave.ss);


    if (dump_bunch_to_file && dump_this_turn(turn)) {
        char filename[50];
        sprintf(filename,"turn%07lu.txt",turn);
        write_bunch_to_file(bun, filename);
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
