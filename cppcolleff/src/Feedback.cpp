#include <cppcolleff/Feedback.h>

double Feedback_t::apply_kick(const my_Dvector& xx_ave, const long turn, Bunch_t& bun)
{
    double kick = 0.0;
    if (track && (turn >= (npoints + delay))) {
        for (int i=1;i<=npoints;i++){
            kick += cos(TWOPI*freq*i + phase) * sin(TWOPI*i/(2*npoints)) / (i*TWOPI/2) * //filter
                    gain * xx_ave[turn-(npoints-i)-delay+1] * sqrt(bpmbeta);
        }
        if      (kick > satur){kick = satur;}
        else if (kick <-satur){kick =-satur;}

        kick *= sqrt(kikbeta);
        for (int p=0;p<bun.num_part;p++){bun.xl[p] += kick;}
    }
    return kick;
}
