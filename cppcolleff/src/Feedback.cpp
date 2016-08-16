#include <cppcolleff/Feedback.h>

double Feedback_t::apply_kick(Bunch_t& bun, const double xx_mean, const double betax)
{
    xx_ave.push_back(xx_mean);
    if (xx_ave.size() > npoints+delay) xx_ave.pop_front();
    double kick = 0.0;
    if (track && (xx_ave.size() == (npoints + delay))) {
        for (unsigned int i=npoints;  i>0; --i){
            kick += cos(TWOPI*freq*i + phase) * sin(TWOPI*i/(2*npoints)) / (i*TWOPI/2) * //filter
                    gain * xx_ave.at(npoints-i) * sqrt(bpmbeta/betax);
        }
        if      (kick > satur){kick = satur;}
        else if (kick <-satur){kick =-satur;}

        kick *= sqrt(kikbeta/betax);
        for (auto& p:bun.particles) {p.xl += kick;}
    }
    return kick;
}
