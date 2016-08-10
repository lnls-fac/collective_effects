#ifndef _WAKE_H
#define _WAKE_H

#include <complex>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

struct WakePl
{
    bool general, resonator;
    my_Dvector  s, W;
    my_Dvector wr, Rs, Q;
    WakePl(): general(false), resonator(false) {};
    // ~WakePl() = default;
};

struct Wake_t
{
    WakePl Wd, Wq, Wl;
    Wake_t() {};
    // ~Wake_t() = default;
    my_Dvector apply_kicks(Bunch_t& bun, const double stren);
};

#endif
