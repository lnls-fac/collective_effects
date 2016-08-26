
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

void Wl_res_kick(
    my_PartVector& p,
    wakePl& Wl,
    double stren, double Wkl)
{
    double&& kr  = Wl.wr[r] / light_speed;
    double&& Ql  = sqrt(Wl.Q[r] * Wl.Q[r] - 0.25);
    double&& Amp = Wl.wr[r] * Wl.Rs[r] / Wl.Q[r] * stren;
    double&& krl = kr * Ql / Wl.Q[r];
    complex<double> cpl_kr (kr/(2*Wl.Q[r]), krl);
    complex<double> W_pot (0.0,0.0);
    for (auto w=0;w<p.size();++w){
        complex<double>&& kik = W_pot * exp(-p[w].ss*cpl_kr);
        W_pot +=                        exp( p[w].ss*cpl_kr);

        double&& kick = - Amp * (0.5 + 1.0*kik.real() + 1.0*kik.imag()/(2*Ql));
        Wkl     += kick;
        p[w].de += kick;
    }
}
void Wd_res_kick(
    my_PartVector& p,
    wakePl& Wd,
    double stren, double betax, double Wkd)
{
    double&& kr  = Wd.wr[r] / light_speed;
    double&& Ql  = sqrt(Wd.Q[r] * Wd.Q[r]  -  0.25);
    double&& Amp = Wd.wr[r] * Wd.Rs[r] / Ql  * stren / betax;
    double&& krl = kr * Ql / Wd.Q[r];
    complex<double> cpl_kr (kr/(2*Wd.Q[r]), krl);
    complex<double> W_pot (0.0,0.0);
    for (auto w=0;w<p.size();++w){
        complex<double>&& kik = W_pot * p[w].xx * exp(-p[w].ss*cpl_kr);
        W_pot +=                                  exp( p[w].ss*cpl_kr);

        double&& kick = -Amp * kik.imag();
        Wkd     += kick;
        p[w].xl += kick;
    }

}

void Wq_res_kick(
    my_PartVector& p,
    wakePl& Wq,
    double stren, double betax, double Wkd)
{
    double&& kr  = Wq.wr[r] / light_speed;
    double&& Ql  = sqrt(Wq.Q[r] * Wq.Q[r]  -  0.25);
    double&& Amp = Wq.wr[r] * Wq.Rs[r] / Ql  * stren / betax;
    double&& krl = kr * Ql / Wq.Q[r];
    complex<double> cpl_kr (kr/(2*Wq.Q[r]), krl);
    complex<double> W_pot (0.0,0.0);
    for (auto w=0;w<p.size();++w){
        complex<double>&& kik = W_pot * exp(-p[w].ss*cpl_kr);
        W_pot +=                        exp( p[w].ss*cpl_kr);

        double&& kick = p[w].xx * (-Amp * kik.imag());
        Wkd     += kick;
        p[w].xl += kick;
    }
}
my_Dvector Wake_t::apply_kicks(Bunch_t& bun, const double stren, const double betax) const
{
    my_Dvector Wkick (2,0.0);
    double Wgl(0.0), Wgd(0.0);
    auto& p = bun.particles;

    //pw --> witness particle   ps --> source particle
    if (Wd.general || Wq.general || Wl.general){
      #ifdef OPENMP
        #pragma omp parallel for schedule(guided,1) reduction(+:Wkl,Wkd)
      #endif
        for (auto w=0;w<p.size();++w){ // Begin from the particle ahead
            // for (auto ps=bun.particles.end()-1;ps!=bun.particles.begin()-1;--ps){ // loop over all particles ahead of it.
            for (auto s=w;s>=0;--s){ // loop over all particles ahead of it.
                double&& ds = p[w].ss - p[s].ss;
                if (Wl.general) {
                    double&& kick = - Wl.W.get_y(ds) * stren;
                    Wgl     += kick;
                    p[w].de += kick;
                }
                if (Wd.general) {
                    double&& kick = -p[s].xx * Wd.W.get_y(ds) * stren / betax; // The kick is the negative of the wake;
                    Wgd     += kick;
                    p[w].xl += kick;
                }
                if (Wq.general) {
                    double&& kick = -p[w].xx * Wq.W.get_y(ds) * stren / betax; // The kick is the negative of the wake;
                    Wgd     += kick;
                    p[w].xl += kick;
                }
            }
        }
    }

    if (Wd.general || Wq.general || Wl.general){
        my_Ivector ready(global_num_threads,0); // the idea is to use this vector to check if the thread is finished
        my_Ivector res(global_num_threads,0); // the idea is to use this vector to check if the thread is finished
    }

    for (int r=0; r<Wl.wr.size(); r++){
        Wl_res_kick(ref(p), ref(Wl), ref(Wrl), stren, r);
    }
    for (int r=0; r<Wd.wr.size(); r++){
        Wd_res_kick(ref(p), ref(Wd), ref(Wkd), stren, betax, r);
    }
    for (int r=0; r<Wq.wr.size(); r++){
        Wq_res_kick(ref(p), ref(Wq), ref(Wkd), stren, betax, r);
    }

    Wkick[0] += Wkl;
    Wkick[1] += Wkd;
    return Wkick;
}
