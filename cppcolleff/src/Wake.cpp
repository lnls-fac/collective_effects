
#include <cppcolleff/Wake.h>

my_Dvector WakePl::get_wake_at_points(const my_Dvector& spos, const double& stren) const
{
    my_Dvector wakeF(spos.size(),0.0);
    if (general) {
        for (int i=0;i<spos.size();++i){
            wakeF[i] = W.get_y(spos[i]) * stren;
        }
    }
    if (resonator){
        for (int r=0; r<wr.size(); r++){
            double&& kr  = wr[r] / light_speed;
            double&& Ql  = sqrt(Q[r] * Q[r] - 0.25);
            double&& Amp = wr[r] * Rs[r] / Q[r] * stren;
            double&& krl (kr * Ql / Q[r]);
            complex<double> cpl_kr (kr/(2*Q[r]), krl);
            complex<double> W_pot (0.0,0.0);
            #ifdef OPENMP
              #pragma omp parallel for schedule(guided,1)
            #endif
            for (int i=0;i<spos.size();++i){
                if (spos[i] < 0.0) {continue;}
                if ((spos[i] < 1e-10) && (spos[i] > -1e-10)) {
                    complex<double>&& kik = exp( -spos[i]*cpl_kr);
                    wakeF[i] += 0.5 * Amp * (1.0*kik.real() + 1.0*kik.imag()/(2*Ql));
                }
                else {
                    complex<double>&& kik = exp( -spos[i]*cpl_kr);
                    wakeF[i] += 1.0 * Amp * (1.0*kik.real() + 1.0*kik.imag()/(2*Ql));
                }
            }
        }
    }
    return wakeF;
}


my_Dvector Wake_t::apply_kicks(Bunch_t& bun, const double stren, const double betax) const
{
    my_Dvector Wkick (2,0.0);

    //pw --> witness particle   ps --> source particle
    if (Wd.general || Wq.general || Wl.general){
        // #ifdef OPENMP
        //   #pragma omp parallel for schedule(guided,1)
        // #endif
        for (auto pw=bun.particles.begin();pw!=bun.particles.end();++pw){ // Begin from the particle ahead
            // for (auto ps=bun.particles.end()-1;ps!=bun.particles.begin()-1;--ps){ // loop over all particles ahead of it.
            for (auto ps=pw;ps!=bun.particles.begin()-1;--ps){ // loop over all particles ahead of it.
                double&& ds = pw->ss - ps->ss;
                if (Wl.general) {
                    double&& kick = - Wl.W.get_y(ds) * stren;
                    // if (kick > 1e-30 || kick <-1e-30){
                        // kick *= stren; // The kick is the negative of the wake;
                        Wkick[0] += kick;
                        pw->de   += kick;
                    // }
                }
                if (Wd.general) {
                    double&& kick = -ps->xx * Wd.W.get_y(ds) * stren / betax; // The kick is the negative of the wake;
                    Wkick[1]   += kick;
                    pw->xl += kick;
                }
                if (Wq.general) {
                    double&& kick = -pw->xx * Wq.W.get_y(ds) * stren / betax; // The kick is the negative of the wake;
                    Wkick[1]   += kick;
                    pw->xl += kick;
                }
            }
        }
    }
    if (Wl.resonator){
        for (int r=0; r<Wl.wr.size(); r++){
            double&& kr  = Wl.wr[r] / light_speed;
            double&& Ql  = sqrt(Wl.Q[r] * Wl.Q[r] - 0.25);
            double&& Amp = Wl.wr[r] * Wl.Rs[r] / Wl.Q[r] * stren;
            double&& krl = kr * Ql / Wl.Q[r];
            complex<double> cpl_kr (kr/(2*Wl.Q[r]), krl);
            complex<double> W_pot (0.0,0.0);
            for (auto& p:bun.particles){
                complex<double>&& kik = W_pot * exp(-p.ss*cpl_kr);
                W_pot +=                        exp( p.ss*cpl_kr);

                double&& kick = - Amp * (0.5 + 1.0*kik.real() + 1.0*kik.imag()/(2*Ql));
                Wkick[0]  += kick;
                p.de += kick;
            }
        }
    }
    if (Wd.resonator){
        for (int r=0; r<Wd.wr.size(); r++){
            double&& kr  = Wd.wr[r] / light_speed;
            double&& Ql  = sqrt(Wd.Q[r]*Wd.Q[r] - 0.25);
            double&& Amp = Wd.wr[r] * Wd.Rs[r] / Ql  * stren / betax;
            double&& krl = kr * Ql / Wd.Q[r];
            complex<double> cpl_kr (kr/(2*Wd.Q[r]), krl);
            complex<double> W_pot (0.0,0.0);
            for (auto& p:bun.particles){
                complex<double>&& kik = W_pot * p.xx * exp(-p.ss*cpl_kr);
                W_pot +=                               exp( p.ss*cpl_kr);

                double&& kick = -Amp * kik.imag();
                Wkick[1] += kick;
                p.xl     += kick;
            }
        }
    }
    if (Wq.resonator){
        for (int r=0; r<Wq.wr.size(); r++){
            double&& kr  = Wq.wr[r] / light_speed;
            double&& Ql  = sqrt(Wq.Q[r]*Wq.Q[r] - 0.25);
            double&& Amp = Wq.wr[r] * Wq.Rs[r] / Ql  * stren / betax;
            double&& krl = kr * Ql / Wq.Q[r];
            complex<double> cpl_kr (kr/(2*Wq.Q[r]), krl);
            complex<double> W_pot (0.0,0.0);
            for (auto& p:bun.particles){
                complex<double>&& kik = W_pot * exp(-p.ss*cpl_kr);
                W_pot +=                        exp( p.ss*cpl_kr);

                double&& kick = p.xx * (-Amp * kik.imag());
                Wkick[1] += kick;
                p.xl     += kick;
            }
        }
    }
    return Wkick;
}
