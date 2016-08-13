#include <cppcolleff/Results.h>

void dump_bunch_to_file(const Bunch_t& bun, char filename) const
{
    FILE* fp = fopen(filename,"w");
    fprintf(fp,"#%20s %20s %20s %20s\n","xx [m]", "xl [m]", "de", "ss [m]");

    for (auto& p:bun.particles){
        fprintf(fp,"%20.7g  %20.7g %20.7g %20.7g\n",p.xx,p.xl,p.de,p.ss);
    }
    fclose(fp);
}

double Results_t::calc_stats(const long turn, const Bunch_t& bun)
{
    if this_turn(turn) { // this_turn is a private method of class Results_t
        ave.push_back(Particle_t (0.0));     rave = ave.back();
        std.push_back(Particle_t (0.0));     rstd = std.back();
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


        if to_file {
            sprintf(char filename[50],"turn%07l.txt",turn);
            dump_bunch_to_file(bun, filename);
        }

        return rave.xx;
    } else {
        double ave_xx =0;
        for (const auto& p:bun.particles){ ave_xx += p.xx;}
        return ave_xx/bun.num_part;
    }
}

void Results_t::set_FBkick(const long turn, const double& kik)
{
    if (FB && this_turn(turn)) {FBkick[get_p(turn)] = kik;}
}

void Results_t::set_Wkicks(const long turn, const my_Dvector& kik)
{
    if (Wl && this_turn(turn)) {Wlkick[get_p(turn)] = kik[0];}
    if (Wd && this_turn(turn)) {Wdkick[get_p(turn)] = kik[1];}
}
