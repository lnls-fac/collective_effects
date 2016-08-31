#ifndef _CPPCOLLEFF_H
#define _CPPCOLLEFF_H

#include <random>  //std::generator and some distributions
#include <utility> //std::swap
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/Results.h>
#include <cppcolleff/Feedback.h>
#include <cppcolleff/Wake.h>
#include <cppcolleff/Ring.h>

void generate_bunch(const Ring_t& ring, Bunch_t& bun);
void generate_bunch(const Ring_t& ring, Bunch_t& bun, unsigned int seed);

my_Dvector solve_Haissinski(const Wake_t& wake,	const Ring_t& ring,	const double& Ib);
double find_equilibrium_energy_spread(const Wake_t& wake, Ring_t& ring,	const double& Ib);

void do_tracking(
    const Ring_t& ring,
    const Wake_t& wake,
    Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results);



#endif
