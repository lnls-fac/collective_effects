
#include <cppcolleff/Wake.h>

my_Dvector Wake_t::apply_kicks(Bunch_t& bun, const double stren)
{
    my_Dvector Wkick (2,0.0);

    //pw --> witness particle   ps --> source particle
    if (Wd.general || Wq.general || Wl.general){
        for (int pw=0;pw<bun.num_part;pw++){ // Begin from the particle ahead
            for (int ps=pw;ps>=0;ps--){ // loop over all particles ahead of it.
                double ds = bun.ss[pw] - bun.ss[ps];
                if (Wl.general) {
                    double kick = interpola(Wl.s,Wl.W,ds) * stren;
                    Wkick[0]   += kick;
                    bun.de[pw] += kick;
                }
                if (Wd.general) {
                    double kick = bun.xx[ps]*interpola(Wd.s,Wd.W,ds) * stren;
                    Wkick[1]   += kick;
                    bun.xl[pw] += kick;
                }
                if (Wq.general) {
                    double kick = bun.xx[pw]*interpola(Wq.s,Wq.W,ds) * stren;
                    Wkick[1]   += kick;
                    bun.xl[pw] += kick;
                }
            }
        }
    }
    if (Wl.resonator){
        for (int r=0; r<Wl.wr.size(); r++){
            double Ql  = sqrt(Wl.Q[r]*Wl.Q[r] - 1/4);
            double wrl = Wl.wr[r] * Ql / Wl.Q[r];
            complex<double> cpl_wr = complex<double> (Wl.wr[r]/(2*Wl.Q[r]), wrl);
            complex<double> W_pot = (0.0,0.0);
            for (int p=0;p<bun.num_part;p++){
                complex<double> kik = W_pot * exp( bun.ss[p]*cpl_wr) * stren;
                W_pot +=                      exp(-bun.ss[p]*cpl_wr);

                double kick = - Wl.wr[r]*Wl.Rs[r]/Wl.Q[r] * (1/2 + kik.real() + 1*kik.imag()/(2*Ql));
                Wkick[0]  += kick;
                bun.de[p] += kick;
            }
        }
    }
    if (Wd.resonator){
        for (int r=0; r<Wd.wr.size(); r++){
            double Ql  = sqrt(Wd.Q[r]*Wd.Q[r] - 1/4);
            double wrl = Wd.wr[r] * Ql / Wd.Q[r];
            complex<double> cpl_wr (Wd.wr[r]/(2*Wd.Q[r]), wrl);
            complex<double> W_pot (0.0,0.0);
            for (int p=0;p<bun.num_part;p++){
                complex<double> kik = W_pot * bun.xx[p] * exp( bun.ss[p]*cpl_wr) * stren;
                W_pot +=                                  exp(-bun.ss[p]*cpl_wr);

                double kick = - Wd.wr[r] * Wd.Rs[r] / Ql * kik.imag();
                Wkick[1]  += kick;
                bun.xl[p] += kick;
            }
        }
    }
    if (Wq.resonator){
        for (int r=0; r<Wq.wr.size(); r++){
            double Ql  = sqrt(Wq.Q[r]*Wq.Q[r] - 1/4);
            double wrl = Wq.wr[r] * Ql / Wq.Q[r];
            complex<double> cpl_wr (Wq.wr[r]/(2*Wq.Q[r]), wrl);
            complex<double> W_pot (0.0,0.0);
            for (int p=0;p<bun.num_part;p++){
                complex<double> kik = W_pot * exp( bun.ss[p]*cpl_wr) * stren;
                W_pot +=                      exp(-bun.ss[p]*cpl_wr);

                double kick = bun.xx[p] * (-Wq.wr[r] * Wq.Rs[r] / Ql * kik.imag());
                Wkick[1]  += kick;
                bun.xl[p] += kick;
            }
        }
    }
    return Wkick;
}
