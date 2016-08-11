#include <cppcolleff/Results.h>


void Results_t::calc_stats(const long turn, const Bunch_t& bun)
{
    double n = get_p(turn);
    //calculate the average and std values;
    for (auto p=bun.particles.cbegin();p!=bun.particles.cend();++p){
        ave[n].xx += p->xx;
        ave[n].xl += p->xl;
        ave[n].de += p->de;
        ave[n].ss += p->ss;
        std[n].xx += p->xx * p->xx;
        std[n].xl += p->xl * p->xl;
        std[n].de += p->de * p->de;
        std[n].ss += p->ss * p->ss;
    }
    ave[n].xx /= bun.num_part;
    ave[n].xl /= bun.num_part;
    ave[n].de /= bun.num_part;
    ave[n].ss /= bun.num_part;
    std[n].xx /= bun.num_part;
    std[n].xl /= bun.num_part;
    std[n].de /= bun.num_part;
    std[n].ss /= bun.num_part;
    std[n].xx -= ave[n].xx * ave[n].xx;
    std[n].xl -= ave[n].xl * ave[n].xl;
    std[n].de -= ave[n].de * ave[n].de;
    std[n].ss -= ave[n].ss * ave[n].ss;
    std[n].xx = sqrt(std[n].xx);
    std[n].xl = sqrt(std[n].xl);
    std[n].de = sqrt(std[n].de);
    std[n].ss = sqrt(std[n].ss);
}

void Results_t::set_FBkick(const long turn, const double& kik)
{
    if (this->FB && (turn % every)==0) {FBkick[get_p(turn)] = kik;}
}

void Results_t::set_Wkicks(const long turn, const my_Dvector& kik)
{
    if (this->Wl && (turn % every)==0) {Wlkick[get_p(turn)] = kik[0];}
    if (this->Wd && (turn % every)==0) {Wdkick[get_p(turn)] = kik[1];}
}
