#ifndef _RING_H
#define _RING_H

#include <cmath>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

class Ring_t{
  public:
    double betax, alphax, etax, etaxl, gammax, // optical functions
           tunex, chromx, tunex_shift,         // tune related parameters
           circum, mom_comp, T0, energy;       // general ring parameters
    my_Dvector cav_s, cav_V;                   // cavity parameters

    //const double gammax = (1 + alphax*alphax) / betax;
    Ring_t():
        betax(1.0), alphax(0.0), etax(0.0), etaxl(0.0), gammax(1.0),
        tunex(0.0), chromx(0.0), tunex_shift(0.0),
        circum(0.0), mom_comp(0.0), T0(0.0), energy(0.0) {};
    ~Ring_t() = default;
    void track_one_turn(Bunch_t& bun) const;
};

#endif
