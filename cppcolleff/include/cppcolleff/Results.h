#ifndef _RESULTS_H
#define _RESULTS_H

#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

class Results_t {
    private:
        unsigned long calc_every, print_every, save_bunch_every, save_distributions_every, nturns;
        bool FB, Wd, Wl, Wq;

        bool calc_this_turn(const long n) const
        {
            if (calc_every == 0) return false;
            if ((n % calc_every)==0 || n==nturns) return true;
            return false;
        }
        bool print_this_turn(const long n) const
        {
            if (!calc_this_turn(n) || !print_in_screen || print_every == 0) return false;
            if ((n % (calc_every*print_every))==0 || n==nturns) return true;
            return false;
        }
        bool save_bunch_this_turn(const long n) const
        {
            if (!save_bunch || save_bunch_every == 0) return false;
            if ((n % save_bunch_every)==0 || n==nturns) return true;
            return false;
        }
        bool save_distribution_this_turn(const long n) const
        {
            if (save_distributions_every == 0) return false;
            if (!save_distribution_xx && !save_distribution_xx &&
                !save_distribution_de && !save_distribution_ss) return false;
            if ((n % save_distributions_every)==0 || n==nturns) return true;
            return false;
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
        bool save_bunch, print_in_screen, save_distribution_xx;
        bool save_distribution_xl, save_distribution_de, save_distribution_ss;
        my_PartVector ave;
        my_PartVector std;
        my_Dvector Wlkick;
        my_Dvector Wdkick;
        my_Dvector Wqkick;
        my_Dvector FBkick;
        Results_t (const unsigned long nt):
            calc_every(1L), print_every(10L), save_bunch_every(0L), save_distributions_every(0L),
            nturns(nt), save_bunch(false), print_in_screen(true), save_distribution_xx(false),
            save_distribution_xl(false), save_distribution_de(false), save_distribution_ss(false),
            FB(false), Wd(false), Wq(false), Wl(false)      {reserve_memory();}
        Results_t (const unsigned long nt, const unsigned long eve):
            calc_every(eve), print_every(10L), save_bunch_every(0L), save_distributions_every(0L),
            nturns(nt), save_bunch(false), print_in_screen(true), save_distribution_xx(false),
            save_distribution_xl(false), save_distribution_de(false), save_distribution_ss(false),
            FB(false), Wd(false), Wq(false), Wl(false)      {reserve_memory();}
        Results_t (const unsigned long nt, const unsigned long eve, const bool kicks):
            calc_every(eve), print_every(10L), save_bunch_every(0L), save_distributions_every(0L),
            nturns(nt), save_bunch(false), print_in_screen(true), save_distribution_xx(false),
            save_distribution_xl(false), save_distribution_de(false), save_distribution_ss(false),
            FB(kicks), Wd(kicks), Wq(kicks), Wl(kicks)      {reserve_memory();}
        ~Results_t() = default;

        void set_keepFB(const bool keep)    {FB = keep; reserve_memory();}
        void set_keepWd(const bool keep)    {Wd = keep; reserve_memory();}
        void set_keepWq(const bool keep)    {Wq = keep; reserve_memory();}
        void set_keepWl(const bool keep)    {Wl = keep; reserve_memory();}
        void set_nturns(const long nt){nturns = nt; reserve_memory();}
        void set_calc_every(const unsigned long eve){calc_every = eve; reserve_memory();}
        void set_print_every(const unsigned long eve){print_every = eve;}
        void set_save_bunch_every(const unsigned long eve){save_bunch_every = eve;}
        void set_save_distributions_every(const unsigned long eve){save_distributions_every = eve;}

        unsigned long get_nturns() const {return nturns;}
        unsigned long get_calc_every() const {return calc_every;}
        unsigned long get_print_every() const {return print_every;}
        unsigned long get_save_bunch_every() const {return save_bunch_every;}
        unsigned long get_save_distributions_every() const {return save_distributions_every;}
        double calc_stats(const Bunch_t& bun, const long turn, ThreadPool& pool);
        double calc_stats(const Bunch_t& bun, const long turn);
        void register_Wkicks(const long turn, const my_Dvector& kik);
        void register_FBkick(const long turn, const double& kik);

        void to_file(const char* filename) const;
        void from_file(const char* filename);
};

#endif
