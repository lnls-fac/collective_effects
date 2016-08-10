#ifndef _RESULTS_H
#define _RESULTS_H

#include <cmath>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

class Results_t {
    bool to_file, FB, Wd, Wl;
    long every;
    long get_p(const long n) {
        return n/every;
    }
  public:
    long nturns;
    my_Dvector xx_ave;
    my_Dvector xx_std;
    my_Dvector xl_ave;
    my_Dvector xl_std;
    my_Dvector de_ave;
    my_Dvector de_std;
    my_Dvector ss_ave;
    my_Dvector ss_std;
    my_Dvector Wlkick;
    my_Dvector Wdkick;
    my_Dvector FBkick;
    Results_t (const long nt):
        nturns(nt), every(1L), to_file(false), FB(false), Wd(false),Wl(false),
        xx_ave(nt,0.0),  xl_ave(nt,0.0),  de_ave(nt,0.0),  ss_ave(nt,0.0),
        xx_std(nt,0.0),  xl_std(nt,0.0),  de_std(nt,0.0),  ss_std(nt,0.0) {}
    Results_t (const long nt, const long eve):
        nturns(nt), every(eve), to_file(false), FB(false), Wd(false),Wl(false),
        xx_ave(nt,0.0),  xl_ave(nt,0.0),  de_ave(nt,0.0),  ss_ave(nt,0.0),
        xx_std(nt,0.0),  xl_std(nt,0.0),  de_std(nt,0.0),  ss_std(nt,0.0){}
    // ~Results_t() = default;

    void calc_stats(const long turn, const Bunch_t& bun);
    void set_Wkicks(const long turn, const my_Dvector& kik);
    void set_FBkick(const long turn, const double& kik);
};


#endif
