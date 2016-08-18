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
                        "<xx> [mm]","std(xx) [mm]","<xl> [mm]","std(xl) [mm]",
                        "<de> [%]","std(de) [%]","<ss> [mm]","std(ss) [mm]"
                    );
            }
            if (this_turn_print(turn)){
                fprintf(stdout,"%07lu %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f \n",turn,
                        1e3*ave.back().xx,1e3*std.back().xx,
                        1e3*ave.back().xl,1e3*std.back().xl,
                        1e2*ave.back().de,1e2*std.back().de,
                        1e3*ave.back().ss,1e3*std.back().ss
                    );
            }
        }

        return rave.xx;
    } else {
        double ave_xx = 0;
        for (const auto& p:bun.particles){ ave_xx += p.xx;}
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
