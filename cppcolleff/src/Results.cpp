#include <cppcolleff/Results.h>

void Results_t::dump_bunch_to_file(const Bunch_t& bun, const char* filename) const
{
    FILE* fp = fopen(filename,"w");
    fprintf(fp,"#%20s %20s %20s %20s\n","xx [m]", "xl [m]", "de", "ss [m]");

    for (auto& p:bun.particles){
        fprintf(fp,"%20.7g  %20.7g %20.7g %20.7g\n",p.xx,p.xl,p.de,p.ss);
    }
    fclose(fp);
}

double Results_t::calc_stats(const Bunch_t& bun, const long turn)
{
    if (this_turn(turn)) { // this_turn is a private method of class Results_t
        ave.push_back(Particle_t (0.0));     auto& rave = ave.back();
        std.push_back(Particle_t (0.0));     auto& rstd = std.back();
      #ifdef OPENMP
        const my_PartVector& par = bun.particles;
        double xx(0.0), xl(0.0), de(0.0), ss(0.0), xx2(0.0), xl2(0.0), de2(0.0), ss2(0.0);
        #pragma omp parallel for reduction(+:xx,xl,de,ss,xx2,xl2,de2,ss2) schedule(guided,1)
        for (int i=0;i<par.size();++i){
            xx  += par[i].xx;
            xl  += par[i].xl;
            de  += par[i].de;
            ss  += par[i].ss;
            xx2 += par[i].xx*par[i].xx;
            xl2 += par[i].xl*par[i].xl;
            de2 += par[i].de*par[i].de;
            ss2 += par[i].ss*par[i].ss;
        }
        rave.xx = xx;        rave.xl = xl;        rave.de = de;        rave.ss = ss;
        rstd.xx = xx2;       rstd.xl = xl2;       rstd.de = de2;       rstd.ss = ss2;
      #else
        for (const auto& p:bun.particles){
            rave.xx += p.xx;
            rave.xl += p.xl;
            rave.de += p.de;
            rave.ss += p.ss;
            rstd.xx += p.xx*p.xx;
            rstd.xl += p.xl*p.xl;
            rstd.de += p.de*p.de;
            rstd.ss += p.ss*p.ss;
        }
      #endif
        rave.xx /= bun.num_part;
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


        if (to_file) {
            char filename[50];
            sprintf(filename,"turn%07lu.txt",turn);
            dump_bunch_to_file(bun, filename);
        }

        if (print_screen) {
            if (turn == 0) {
                fprintf(stdout,"%7s %12s %12s   %12s %12s   %12s %12s   %12s %12s \n","turn",
                        "<xx> [um]","std(xx) [um]","<xl> [urad]","std(xl) [urad]",
                        "<de> [%]","std(de) [%]","<ss> [mm]","std(ss) [mm]"
                    );
            }
            if (this_turn_print(turn)){
                fprintf(stdout,"%07lu %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f \n",turn,
                        1e6*ave.back().xx,1e6*std.back().xx,
                        1e6*ave.back().xl,1e6*std.back().xl,
                        1e2*ave.back().de,1e2*std.back().de,
                        1e3*ave.back().ss,1e3*std.back().ss
                    );
            }
        }

        return rave.xx;
    } else {
        double ave_xx = 0;
      #ifdef OPENMP
        const my_PartVector& par = bun.particles;
        #pragma omp parallel for
        for (int i=0;i<par.size();++i){ ave_xx += par[i].xx;}
      #else
        for (const auto& p:bun.particles){ ave_xx += p.xx;}
      #endif
        return ave_xx/bun.num_part;
    }
}

void Results_t::register_FBkick(const long turn, const double& kik)
{
    if (FB && this_turn(turn)) {FBkick.push_back(kik);}
}

void Results_t::register_Wkicks(const long turn, const my_Dvector& kik)
{
    if (Wl && this_turn(turn)) {Wlkick.push_back(kik[0]);}
    if (Wd && this_turn(turn)) {Wdkick.push_back(kik[1]);}
}
