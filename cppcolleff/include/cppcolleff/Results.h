#ifndef _RESULTS_H
#define _RESULTS_H

#include <cmath> //std::sqrt std::sin std::cos
#include <cstdio> //std::sprintf std::fprintf
#include <string> //std::string.c_str()
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/ThreadPool/ThreadPool.h>

class Results_t {
    private:
        unsigned long calc_every, print_every, dump_every, nturns;
        bool FB, Wd, Wl, Wq;

        bool calc_this_turn(const long n) const
        {
            if (calc_every == 0) return false;
            if ((n % calc_every)==0 || n==nturns) return true;
            else return false;
        }
        bool print_this_turn(const long n) const
        {
            if (!calc_this_turn(n) || print_every == 0) return false;
            if ((n % (calc_every*print_every))==0 || n==nturns) return true;
            else return false;
        }
        bool dump_this_turn(const long n) const
        {
            if (!calc_this_turn(n) || dump_every == 0) return false;
            if ((n % (calc_every*dump_every))==0 || n==nturns) return true;
            else return false;
        }
        void reserve_memory()
        {
            long np = nturns/calc_every + 2;
            ave.reserve(np);
            std.reserve(np);
            if (FB) FBkick.reserve(np); else FBkick.reserve(0);
            if (Wd) Wdkick.reserve(np); else Wdkick.reserve(0);
            if (Wq) Wqkick.reserve(np); else Wqkick.reserve(0);
            if (Wl) Wlkick.reserve(np); else Wlkick.reserve(0);
        }
    public:
        bool dump_bunch_to_file, print_in_screen;
        my_PartVector ave;
        my_PartVector std;
        my_Dvector Wlkick;
        my_Dvector Wdkick;
        my_Dvector Wqkick;
        my_Dvector FBkick;
        Results_t (const unsigned long nt):
            calc_every(1L), print_every(10L), dump_every(0L),
            nturns(nt), dump_bunch_to_file(false), print_in_screen(true),
            FB(false), Wd(false), Wq(false), Wl(false)      {reserve_memory();}
        Results_t (const unsigned long nt, const unsigned long eve):
            calc_every(eve), print_every(10L), dump_every(0L),
            nturns(nt), dump_bunch_to_file(false), print_in_screen(true),
            FB(false), Wd(false), Wq(false), Wl(false)      {reserve_memory();}
        Results_t (const unsigned long nt, const unsigned long eve, const bool kicks):
            calc_every(eve), print_every(10L), dump_every(0L),
            nturns(nt), dump_bunch_to_file(false), print_in_screen(true),
            FB(kicks), Wd(kicks), Wq(false), Wl(kicks)      {reserve_memory();}
        ~Results_t() = default;

        void set_keepFB(const bool keep)    {FB = keep; reserve_memory();}
        void set_keepWd(const bool keep)    {Wd = keep; reserve_memory();}
        void set_keepWq(const bool keep)    {Wq = keep; reserve_memory();}
        void set_keepWl(const bool keep)    {Wl = keep; reserve_memory();}
        void set_nturns(const long nt){nturns = nt; reserve_memory();}
        void set_calc_every(const unsigned long eve){calc_every = eve; reserve_memory();}
        void set_print_every(const unsigned long eve){print_every = eve;}
        void set_dump_to_file_every(const unsigned long eve){dump_every = eve;}

        unsigned long get_nturns() const {return nturns;}
        unsigned long get_calc_every() const {return calc_every;}
        unsigned long get_print_every() const {return print_every;}
        unsigned long get_dump_every() const {return dump_every;}
        double calc_stats(const Bunch_t& bun, const long turn);
        void register_Wkicks(const long turn, const my_Dvector& kik);
        void register_FBkick(const long turn, const double& kik);
        void write_bunch_to_file(const Bunch_t& bun, const char* filename) const;
};


#endif
