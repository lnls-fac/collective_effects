#ifndef _FEEDBACK_H
#define _FEEDBACK_H

#include <cmath>  // std::sin std::cos std::sqrt
#include <deque>  // class std::deque
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

class Feedback_t {
private:
    deque<double> xx_ave;
    bool track;
    unsigned int npoints, delay;
    double phase, freq, gain, satur, bpmbeta, kikbeta;
    Feedback_t(): track(false) {};
    ~Feedback_t() = default;
    double apply_kick(const double xx_mean, Bunch_t& bun) const;
};

#endif
