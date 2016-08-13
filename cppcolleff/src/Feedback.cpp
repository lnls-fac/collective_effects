#include <cppcolleff/Feedback.h>

double Feedback_t::apply_kick(const double xx_mean, Bunch_t& bun) const
{
    xx_ave.push_back(xx_mean);
    if (xx_ave.size() > npoints+delay) xx_ave.pop_front();
    double kick = 0.0;
    if (track && (xx_ave.size() == (npoints + delay))) {
        for (int i=npoints, const auto x= xx_ave.cbegin;  i>=1;  --i,++x){
            kick += cos(TWOPI*freq*i + phase) * sin(TWOPI*i/(2*npoints)) / (i*TWOPI/2) * //filter
                    gain * (*xx_ave) * sqrt(bpmbeta);
        }
        if      (kick > satur){kick = satur;}
        else if (kick <-satur){kick =-satur;}

        kick *= sqrt(kikbeta);
        for (auto p:bun.particles)){p.xl += kick;}
    }
    return kick;
}
