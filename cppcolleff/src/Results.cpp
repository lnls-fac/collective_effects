#include <cppcolleff/Results.h>


void Results_t::calc_stats(const long turn, const Bunch_t& bun)
{
    double n = get_p(turn);
    //calculate the average and std values;
    for (int p=0;p<=bun.num_part;p++){
        xx_ave[n] += bun.xx[p];
        xl_ave[n] += bun.xl[p];
        de_ave[n] += bun.de[p];
        ss_ave[n] += bun.ss[p];
        xx_std[n] += bun.xx[p] * bun.xx[p];
        xl_std[n] += bun.xl[p] * bun.xl[p];
        de_std[n] += bun.de[p] * bun.de[p];
        ss_std[n] += bun.ss[p] * bun.ss[p];
    }
    xx_ave[n] /= bun.num_part;
    xl_ave[n] /= bun.num_part;
    de_ave[n] /= bun.num_part;
    ss_ave[n] /= bun.num_part;
    xx_std[n] /= bun.num_part;
    xl_std[n] /= bun.num_part;
    de_std[n] /= bun.num_part;
    ss_std[n] /= bun.num_part;
    xx_std[n] -= xx_ave[n] * xx_ave[n];
    xl_std[n] -= xl_ave[n] * xl_ave[n];
    de_std[n] -= de_ave[n] * de_ave[n];
    ss_std[n] -= ss_ave[n] * ss_ave[n];
    xx_std[n] = sqrt(xx_std[n]);
    xl_std[n] = sqrt(xl_std[n]);
    de_std[n] = sqrt(de_std[n]);
    ss_std[n] = sqrt(ss_std[n]);
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
