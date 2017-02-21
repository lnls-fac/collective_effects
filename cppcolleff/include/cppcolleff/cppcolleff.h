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

my_Dvector solve_Haissinski_get_potential(
    const Wake_t& wake,
    const Ring_t& ring,
    const double& Ib,
    const int niter,
    const my_Dvector distr_ini);
my_Dvector solve_Haissinski_get_potential(
    const Wake_t& wake,
    const Ring_t& ring,
    const double& Ib);
my_Dvector solve_Haissinski_get_potential(
    const Wake_t& wake,
    const Ring_t& ring,
    const double& Ib,
    const int niter);
my_Dvector solve_Haissinski_get_potential(
    const Wake_t& wake,
    const Ring_t& ring,
    const double& Ib,
    const my_Dvector distr_ini);

double find_equilibrium_energy_spread(
    const Wake_t& wake,
    Ring_t& ring, // Only changes the energy spread internally;
    const double& Ib,
    const int niter,
    const my_Dvector distr_ini);
double find_equilibrium_energy_spread(
	const Wake_t& wake,
	Ring_t& ring, // Only changes the energy spread internally;
	const double& Ib);
double find_equilibrium_energy_spread(
    const Wake_t& wake,
    Ring_t& ring, // Only changes the energy spread internally;
    const double& Ib,
    const int niter);
double find_equilibrium_energy_spread(
    const Wake_t& wake,
    Ring_t& ring, // Only changes the energy spread internally;
    const double& Ib,
    const my_Dvector distr_ini);

void single_bunch_tracking(
    const Ring_t& ring,
    const Wake_t& wake,
    Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results);

void multi_bunch_tracking(
	const Ring_t& ring,
	const Wake_t& long_range_wake,
	const Wake_t& short_range_wake,
	Feedback_t& fb,
	vector<Bunch_t>& buns,
	vector<Results_t>& results);


#endif
