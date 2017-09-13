#ifndef _CPPCOLLEFF_H
#define _CPPCOLLEFF_H

#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/Results.h>
#include <cppcolleff/Feedback.h>
#include <cppcolleff/Wake.h>
#include <cppcolleff/Ring.h>

my_Dvector solve_Haissinski_get_potential(
    const Wake_t& wake,
    const Ring_t& ring,
    const double& Ib,
    const bool req_convergence = true,
    const int niter = 100,
    const double weight = 3,
    const my_Dvector distr_ini = my_Dvector ());

my_Dvector solve_Haissinski_get_potential(
    const Wake_t& wake,
    const Ring_t& ring,
    const double& Ib,
    ThreadPool& pool,
    Convolve_t& conv,
    const bool req_convergence = true,
    const int niter = 100,
    const double weight = 3,
    const my_Dvector distr_ini = my_Dvector ());

double find_equilibrium_energy_spread(
    const Wake_t& wake,
    Ring_t& ring, // Only changes the energy spread internally;
    const double& Ib,
    const int niter = 100,
    const double weight = 3,
    const my_Dvector distr_ini = my_Dvector ());

double find_equilibrium_energy_spread(
    const Wake_t& wake,
    Ring_t& ring, // Only changes the energy spread internally;
    const double& Ib,
    ThreadPool& pool,
    Convolve_t& conv,
    const int niter = 100,
    const double weight = 3,
    const my_Dvector distr_ini = my_Dvector ());

my_Dvector long_simul_with_haissinki(
	const Wake_t& wake,
	Ring_t& ring,  // Only changes the energy spread internally;
	my_Dvector& currs,
    my_Dvector& init_espread,
	int niter);

void single_bunch_tracking(
    const Ring_t& ring,
    Wake_t& wake,
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
