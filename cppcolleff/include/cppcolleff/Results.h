#ifndef _RESULTS_H
#define _RESULTS_H

#include <cmath> //std::sqrt std::sin std::cos
#include <cstdio> //std::sprintf std::fprintf
#include <string> //std::string.c_str()
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

class Results_t {
    long every, nturns;
    bool FB, Wd, Wl;
    bool this_turn(const long n) const
    {
        if ((n % every)==0 || n==nturns) return true;
        else return false;
    }
    void reserve_memory()
    {
        long np = nt/every + 2;
        ave.reserve(np);
        std.reserve(np);
        if (FB) FBKick.reserve(np); else FBKick.reserve(0);
        if (Wd) WdKick.reserve(np); else WdKick.reserve(0);
        if (Wl) WlKick.reserve(np); else WlKick.reserve(0);
    }
  public:
    bool to_file;
    my_PartVector ave;
    my_PartVector std;
    my_Dvector Wlkick;
    my_Dvector Wdkick;
    my_Dvector FBkick;
    Results_t (const long nt):
        nturns(nt), every(1L), to_file(false),
        FB(false), Wd(false),Wl(false)
    {
        reserve_memory();
    }
    Results_t (const long nt, const long eve):
        nturns(nt), every(eve), to_file(false),
        FB(false), Wd(false),Wl(false)
    {
        reserve_memory();
    }
    Results_t (const long nt, const long eve, const bool kicks):
        nturns(nt), every(eve), to_file(false),
        FB(kicks), Wd(kicks), Wl(kicks)
    {
        reserve_memory();
    }
    ~Results_t() = default;
    void set_keepFB(bool keep)    {FB     = keep; reserve_memory();}
    void set_keepWd(bool keep)    {Wd     = keep; reserve_memory();}
    void set_keepWl(bool keep)    {Wl     = keep; reserve_memory();}
    void set_every(const long eve){every  = eve;  reserve_memory();}
    void set_nturns(const long nt){nturns = nt;   reserve_memory();}

    long get_nturns() const {return nturns;}
    double calc_stats(const long turn, const Bunch_t& bun);
    void register_Wkicks(const long turn, const my_Dvector& kik);
    void register_FBkick(const long turn, const double& kik);
    void dump_bunch_to_file(const Bunch_t& bun, string& filename) const;
};


#endif
