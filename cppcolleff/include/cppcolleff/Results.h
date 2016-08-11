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
    my_PartVector ave;
    my_PartVector std;
    my_Dvector Wlkick;
    my_Dvector Wdkick;
    my_Dvector FBkick;
    Results_t (const long nt):
        nturns(nt), every(1L), to_file(false), FB(false), Wd(false),Wl(false),
        ave(nt,0.0), std(nt,0.0) {};
    Results_t (const long nt, const long eve):
        nturns(nt), every(eve), to_file(false), FB(false), Wd(false),Wl(false),
        ave(nt,0.0), std(nt,0.0) {};
    ~Results_t() = default;

    void calc_stats(const long turn, const Bunch_t& bun);
    void set_Wkicks(const long turn, const my_Dvector& kik);
    void set_FBkick(const long turn, const double& kik);
};


#endif
