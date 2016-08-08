
#include <vector>
#include <complex>

const double light_speed              = 299792458;         // [m/s]   - definition

static void InsertionSort(
    std::vector<double> &xx, std::vector<double> &xl,
    std::vector<double> &de, std::vector<double> &ss){

	unsigned int i, j;
    double atualxx, atualxl, atualde, atualss;

	for (i=1; i < xx.size(); i++){
		atualxx = xx[i];
        atualxl = xl[i];
        atualde = de[i];
        atualss = ss[i];
		j = i - 1;

		while ((j>=0) && (atualss < s[j])){
			ss[j+1] = ss[j];
            xx[j+1] = xx[j];
            xl[j+1] = xl[j];
            de[j+1] = de[j];
            j--;
		}

		ss[j+1] = atualss;
        xx[j+1] = atualxx;
        xl[j+1] = atualxl;
        de[j+1] = atualde;
	}
}

static double interpola(const std::vector<double> &si, const std::vector<double> &Wi, const double &s){
    unsigned long i

    i = (unsigned long) ((s - s[0])/(si[1]-si[0]));
    return (Wi[i+1] - Wi[i]) / (si[i+1] - si[i]);
}

struct wake_t {
    bool general, resonator;
    std::vector<double>  s;
    std::vector<double>  W;
    std::vector<double> wr;
    std::vector<double> Rs;
    std::vector<double>  Q;
}

struct ring_t {
    double betax, alphax, etax, etaxl,
           tunex, chromx, tunex_shift, circum, mom_comp;
}

void do_tracking(
    const param_t ring,
    const unsigned int nturns,
    const wake_t &Wll, const wake_t &Wdx, const wake_t &Wqx){

    std::vector<double> opt_cav_s (N,0.0);
    std::vector<double> opt_cav_V (N,0.0);
    std::vector<double> de (num_part,0.0);
    std::vector<double> xx (num_part,0.0);
    std::vector<double> xl (num_part,0.0);
    std::vector<double> ss (num_part,0.0);
    std::vector<double> Wkickd (nturns,0.0);
    std::vector<double> Wkickx (nturns,0.0);
    std::vector<double> FBkickx (nturns,0.0);
    std::vector<double> x_ave  (nturns,0.0);
    std::vector<double> x_std  (nturns,0.0);
    std::vector<double> xl_ave (nturns,0.0);
    std::vector<double> xl_std (nturns,0.0);
    std::vector<double> d_ave  (nturns,0.0);
    std::vector<double> d_std  (nturns,0.0);
    std::vector<double> s_ave  (nturns,0.0);
    std::vector<double> s_std  (nturns,0.0);

    unsigned int FB_npoints, FB_delay;
    double FB_phase, FB_freq, FB_gain, FB_satur, FB_bpmbeta, FB_kikbeta;

    std::complex<double> W_pot, kik, cpl_wrl;
    double kick, ds, x_tmp, phix, fil, Ql, wrl;
    double kick_stren = opt_T0 / opt_energy * bunch_Ib / bunch_num_part;

    InsertionSort(xx, xl, de, ss)



    for (int n=0;i<nturns;i++){
        // convention: positive ss, means particle behind the sinchronous particle;
        // First do single particle tracking:
        // the ring model is composed by one cavity, then the magnetic elements and finally the wake-field.
        for (int p=0;p<num_part;p++){
            de[p] += interpola(opt_cav_s, opt_cav_V, ss[p])/opt_energy;
            ss[p] += opt_circum * opt_mom_comp * de[p];
            phix  = TWOPI*(opt_tunex +
                           opt_chromx*de[p] +
                           opt_tunex_shift*((xx[p]-opt_etax* de[p])*(xx[p]-opt_etax* de[p])/opt_betax +
                                            (xl[p]-opt_etaxl*de[p])*(xl[p]-opt_etaxl*de[p])*opt_betax));
            x_tmp = opt_etax*de[p] + (xx[p]-opt_etax* de[p])*std::cos(phix) +
                         opt_betax * (xl[p]-opt_etaxl*de[p])*std::sin(phix);
            xl[p] = -1.0/opt_betax * (xx[p]-opt_etax* de[p])*std::sin(phix) +
                   opt_etaxl*de[p] + (xl[p]-opt_etaxl*de[p])*std::cos(phix);
            xx[p] = x_tmp;
        }

        // After this sorting, the particles will be ordered from head to tail.
        // It means, from smaller ss to bigger ss.
        InsertionSort(xx, xl, de, ss)

        //pw --> witness particle   ps --> source particle
        if (Wdx.general || Wqx.general || Wll.general){
            for (int pw=0;pw<num_part;pw++){ // Begin from the particle ahead
                for (int ps=pw;ps>=0;ps--){ // loop over all particles ahead of it.
                    ds = ss[pw] - ss[ps];
                    if (Wdx.general) {
                        kick       = xx[ps]*interpola(Wdx.s,Wdx.W,ds) * kick_stren;
                        Wkickx[n] += kick;
                        xl[pw]    += kick;
                    }
                    if (Wqx.general) {
                        kick       = xx[pw]*interpola(Wqx.s,Wqx.W,ds) * kick_stren;
                        Wkickx[n] += kick;
                        xl[pw]    += kick;
                    }
                    if (Wll.general) {
                        kick       = interpola(Wll.s,Wll.W,ds) * kick_stren;
                        Wkickd[n] += kick;
                        de[pw]    += kick;
                    }
                }
            }
        }
        if (Wll.resonator){
            for int(r=0; r<Wll.wr.size(); r++){
                Ql = std::sqrt(Wll.Q[r]*Wll.Q[r] - 1/4);
                wrl = Wll.wr[r] * Ql / Wll.Q[r];
                cpl_wr = std::complex<double> (Wll.wr[r]/(2*Wll.Q[r]), wrl);
                W_pot = (0.0,0.0);
                for (int p=0;p<num_part;p++){
                    kik    = W_pot * std::complex::exp( ss[p]*cpl_wr) * kick_stren;
                    W_pot +=         std::complex::exp(-ss[p]*cpl_wr);

                    kick = - Wll.wr[r]*Wll.Rs[r]/Wll.Q[r] * (1/2 + kik.real() + 1*kik.imag()/(2*Ql));
                    Wkickd[n] += kick;
                    de[p]     += kick;
                }
            }
        }
        if (Wdx.resonator){
            for int(r=0; r<Wdx.wr.size(); r++){
                Ql = std::sqrt(Wdx.Q[r]*Wdx.Q[r] - 1/4);
                wrl = Wdx.wr[r] * Ql / Wdx.Q[r];
                cpl_wr = std::complex<double> (Wdx.wr[r]/(2*Wdx.Q[r]), wrl);
                W_pot = (0.0,0.0);
                for (int p=0;p<num_part;pw++){
                    kik    = W_pot * xx[p] * std::complex::exp( ss[p]*cpl_wr) * kick_stren;
                    W_pot +=                 std::complex::exp(-ss[p]*cpl_wr);

                    kick = - Wdx.wr[r]*Wdx.Rs[r]/Ql * kik.imag();
                    Wkickd[n] += kick;
                    de[p]     += kick;
                }
            }
        }
        if (Wqx.resonator){
            for int(r=0; r<num_res; r++){
                Ql = sqrt(Wqx.Q[r]*Wqx.Q[r] - 1/4);
                wrl = Wqx.wr[r] * Ql / Wqx.Q[r];
                cpl_wr = std::complex<double> (Wqx.wr[r]/(2*Wqx.Q[r]), wrl);
                W_pot = (0.0,0.0);
                for (int p=0;p<num_part;pw++){
                    kik    = W_pot * std::complex::exp( ss[p]*cpl_wr) * kick_stren;
                    W_pot +=         std::complex::exp(-ss[p]*cpl_wr);

                    kick = xx[p] * (-Wqx.wr[r]*Wqx.Rs[r]/Ql * kik.imag());
                    Wkickd[n] += kick;
                    de[p]    += kick;
                }
            }
        }

        if (FB_track && (n >= (FB_npoints+FB_delay))) {
            kick = 0;
            for (int i=1;i<=FB_npoints;i++){
                fil = std::cos(TWOPI*FB_freq*i + FB_phase) * std::sin(TWOPI*i/(2*FB_npoints)) / (i*PI);
                kick += fil * FB_gain * x_m[n-(FB_npoints-i)-FB_delay+1] * std::sqrt(FB_bpmbeta);
            }
            if      (kick > FB_satur){kick = satur;}
            else if (kick <-FB_satur){kick =-satur;}

            kick *= std::sqrt(FB_kikbeta);
        }
    }
}
