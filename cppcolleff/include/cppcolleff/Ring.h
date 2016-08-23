#ifndef _RING_H
#define _RING_H

#include <cmath> // std::sin std::cos
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
//#include <omp.h>

class Ring_t{
  private:
    my_Dvector _get_distribution(const my_Dvector& spos, const my_Dvector& V) const;
  public:
    int harm_num;                              // harmonic number
    double betax, alphax, etax, etaxl,         // optical functions
           tunex, chromx, tunex_shift,         // tune related parameters
           circum, mom_comp, T0, energy,       // general ring parameters
           emitx, espread;                     // equilibrium parameters
    Interpola_t cav;                           // cavity parameters

    //const double gammax = (1 + alphax*alphax) / betax;
    Ring_t():
        betax(1.0), alphax(0.0), etax(0.0), etaxl(0.0),
        tunex(0.0), chromx(0.0), tunex_shift(0.0),
        circum(0.0), mom_comp(0.0), T0(0.0), energy(0.0),
        emitx(0.0), espread(0.0) {};
    ~Ring_t() = default;
    my_Dvector get_distribution(const my_Dvector& spos, const my_Dvector& V) const;
    my_Dvector get_distribution(const my_Dvector& V) const;
    my_Dvector get_distribution() const;
    Interpola_t get_integrated_distribution() const;
    void track_one_turn(Bunch_t& bun) const;
};

#endif
