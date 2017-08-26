#ifndef _WAKE_H
#define _WAKE_H

#include <thread>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/ThreadPool/ThreadPool.h>

struct WakePl
{
    bool wake_function, wake_potential, resonator;
    Interpola_t  WF, WP;
    my_Dvector wr, Rs, Q;
    WakePl(): wake_function(false), wake_potential(false), resonator(false) {};
    ~WakePl() = default;
    my_Dvector get_wake_at_points(const my_Dvector& spos, const double& stren) const;
};

struct Wake_t
{
    WakePl Wd, Wq, Wl;
    Wake_t() {};
    ~Wake_t() = default;
    my_Dvector apply_kicks(Bunch_t& bun, const double stren, const double betax) const;
};

#endif
