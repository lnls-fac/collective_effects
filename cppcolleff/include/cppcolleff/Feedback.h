#ifndef _FEEDBACK_H
#define _FEEDBACK_H

#include <cmath>  // std::sin std::cos std::sqrt
#include <deque>  // class std::deque
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

class Feedback_t {
  private:
    deque<double> xx_ave;
  public:
    bool track;
    unsigned int npoints, delay;
    double phase, freq, gain, satur, bpmbeta, kikbeta;
    Feedback_t(): track(false) {};
    ~Feedback_t() = default;
    double apply_kick(Bunch_t& bun, const double xx_mean, const double betax);
};

#endif
