#ifndef _FEEDBACK_H
#define _FEEDBACK_H

#include <cmath>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

struct Feedback_t {
    bool track;
    unsigned int npoints, delay;
    double phase, freq, gain, satur, bpmbeta, kikbeta;
    Feedback_t(): track(false) {};
    ~Feedback_t() = default;
    double apply_kick(const my_Dvector& xx_ave,const long nturns, Bunch_t& bun) const;
};

#endif
