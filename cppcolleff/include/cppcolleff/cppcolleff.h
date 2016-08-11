#ifndef _CPPCOLLEFF_H
#define _CPPCOLLEFF_H

#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/Results.h>
#include <cppcolleff/Feedback.h>
#include <cppcolleff/Wake.h>
#include <cppcolleff/Ring.h>

void do_tracking(
    const Ring_t& ring,
    const Wake_t& wake,
    const Feedback_t& fb,
    Bunch_t& bun,
    Results_t& results);

#endif
