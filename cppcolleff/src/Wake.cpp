#include <cppcolleff/Wake.h>

my_Dvector WakePl::get_wake_at_points(const my_Dvector& spos, const double& stren) const
{
    my_Dvector wakeF(spos.size(),0.0);
    if (wake_function) {
        for (int i=0;i<spos.size();++i){
            wakeF[i] = WF.get_y(spos[i]) * stren;
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
            //   #pragma omp parallel for schedule(guided,1)
            #endif
            for (int i=0;i<spos.size();++i){
                if (spos[i] < 0.0) {continue;}
                complex<double>&& kik = exp( -spos[i]*cpl_kr);
                wakeF[i] += Amp * (1.0*kik.real() + 1.0*kik.imag()/(2*Ql));
                if (abs(spos[i])<1e-13) { wakeF[i] /= 2; } //# a particle feels half of its wake
            }
        }
    }
    return wakeF;
}

my_Dvector Wake_t::apply_wake_function_kick(
    my_PartVector& par,
    double stren,
    double strenT,
    ThreadPool& pool) const
{
    auto fun = [&](my_PartVector& p, int com, int ter)
    {
        my_Dvector vec (3,0.0);
        for (auto w=com;w<ter;++w){ // Begin from the particle ahead
            if (Wl.wake_function) {
                double&& kick = Wl.WF.get_y(0.0);
                kick *= -stren;
                vec[0] += kick;
                p[w].de += kick;
            }
            for (auto s=w-1;s>=0;--s){ // loop over all particles ahead of it.
                double&& ds = p[w].ss - p[s].ss;
                if (Wl.wake_function) {
                    double&& kick = Wl.WF.get_y(ds);
                    kick *= -stren;
                    vec[0] += kick;
                    p[w].de += kick;
                }
                if (Wd.wake_function) {
                    double&& kick = Wd.WF.get_y(ds);
                    kick *= -p[s].xx * strenT; // The kick is the negative of the wake;
                    vec[1] += kick;
                    p[w].xl += kick;
                }
                if (Wq.wake_function) {
                    double&& kick = Wq.WF.get_y(ds);
                    kick *= -p[w].xx * strenT; // The kick is the negative of the wake;
                    vec[2] += kick;
                    p[w].xl += kick;
                }
            }
        }
        return vec;
    };

    // auto&& off = ((i-1) % 2 == 0) ? (nr_th-(i/2)+1) : (i/2);
    // res1.emplace_back(pool.enqueue(fun, ref(par), lims2[off-(i/2)], lims2[off-(i/2)+1]));

    unsigned int nr_th = 64*get_num_threads();
    // unsigned int&& nr_job = max(int(par.size()/100), min(int(nr_th), int(par.size())));
    vector< std::future<my_Dvector> > res1;

    my_Ivector lims2 (get_bounds(par.size()/2, par.size(), nr_th));
    for (unsigned int i=1;i<=nr_th;++i){
        // auto&& off = (i%2 == 0) ? (nr_th-(i/2)-1) : (i/2);
        // res1.emplace_back(pool.enqueue(fun, ref(par), lims2[off], lims2[off+1]));

        res1.emplace_back(pool.enqueue(fun, ref(par), lims2[nr_th-i], lims2[nr_th-i+1]));
    }

    my_Ivector lims1 (get_bounds(0, par.size()/2, nr_th));
    for (unsigned int i=0;i<nr_th;++i){
        res1.emplace_back(pool.enqueue(fun, ref(par), lims1[i], lims1[i+1]));
    }

    my_Dvector Wgs (3,0.0);
    my_Dvector v (3,0.0);
    for(auto&& result: res1)
        {v = result.get(); Wgs[0] += v[0]; Wgs[1] += v[1]; Wgs[2] += v[2];}
    return Wgs;
}


double W_res_kick_threads(
    my_PartVector& p,
    double& Amp,
    complex<double>& cpl_kr,
    double& Ql,
    int Ktype, // 0 for longituinal, 1 dipolar, 2 for quadrupolar
    my_Cvector& W_pot,
    my_Ivector& lims,
    int i,
    bool ch_W_pot)
{
    double tot_kick (0.0);

    for (auto w=lims[i];w<lims[i+1];++w){
        complex<double>&& ex = exp(-p[w].ss*cpl_kr);
        complex<double>&& kik = W_pot[i] * ex;
        // in the first round of threads I have to calculate the potential
        // in the cavity, but in the second, I haven't.
        if (ch_W_pot) {
            if (Ktype==1){W_pot[i] += p[w].xx / ex;} //dip
            else         {W_pot[i] +=   1.0   / ex;}
        }

        if (Ktype==0){ // longitudinal
            double&& kick = - Amp * ( 0.5 + kik.real() + kik.imag()/(2.0*Ql) );
            tot_kick += kick;
            p[w].de += kick;
        }
        else if (Ktype==1){ // dipolar
            double&& kick = -Amp * kik.imag();
            tot_kick += kick;
            p[w].xl += kick;
        }
        else {  // quadrupolar
            double&& kick = -Amp * kik.imag() * p[w].xx;
            tot_kick += kick;
            p[w].xl += kick;
        }
    }
    return tot_kick;
}

double Wake_t::apply_wake_resonator_kick(
    my_PartVector& p,
    int Ktype, // 0 for longituinal, 1 dipolar, 2 for quadrupolar
    double stren,
    my_Ivector& lims,
    ThreadPool& pool) const
{
    double Total_Kick;

    // First I get the size of the pool of threads:
    int nr_th = lims.size()-1;

    const WakePl& W = (Ktype == LL) ? Wl: ((Ktype == XD) ? Wd : Wq);

    for (int r=0; r<W.wr.size(); r++){
        // Now I calculate some important variables'
        double&& kr  = W.wr[r] / light_speed;
        double&& Ql  = sqrt(W.Q[r] * W.Q[r] - 0.25);
        double&& krl = kr * Ql / W.Q[r];
        // these two are the only ones which will actually be used:
        double Amp;
        if (Ktype==0) Amp = W.wr[r] * W.Rs[r] / W.Q[r] * stren;
        else          Amp = W.wr[r] * W.Rs[r] / Ql  * stren;
        complex<double> cpl_kr (kr/(2*W.Q[r]), krl);

        // I need to define some variables to be used as temporaries in the parallelization:
        my_Cvector W_pot (nr_th,(0.0,0.0)); // This one will keep the wake phasors
        bool ch_W_pot (true); // flow control to be passed to the threads
        my_Dvector Kick (nr_th,0.0); // This one will keep the kicks received by the particles;

        vector< std::future<double> > res1, res2;
        // submit the first round of threads:
        for (unsigned int i=0;i<=nr_th;++i){
            res1.emplace_back(pool.enqueue(
                W_res_kick_threads, ref(p), ref(Amp), ref(cpl_kr), ref(Ql),
                Ktype, ref(W_pot), ref(lims), i, ch_W_pot));
        }

        // Now I have to prepare the second round of threads.
        Total_Kick += res1[0].get(); // get results of job 0 first;
        ch_W_pot = false; // I don't want the potential W_pot to be updated;
        complex<double> W_pot_sum (W_pot[0]); // temporary variables
        for(int i=1;i<nr_th;++i){
            Total_Kick += res1[i].get(); // get results of job i;
            complex<double> temp_var (W_pot[i]); // create the new wake potential
            W_pot[i]   = W_pot_sum;              // to be applied in the particles
            W_pot_sum += temp_var;
            res2.emplace_back(pool.enqueue(
                W_res_kick_threads, ref(p), ref(Amp), ref(cpl_kr), ref(Ql),
                Ktype, ref(W_pot), ref(lims), i, ch_W_pot));
        }
        // join the second round of threads and update the total kick received by the particles
        for (int i=0;i<res2.size();++i){
            Total_Kick += res2[i].get();
        }
    }
    return Total_Kick;
}


my_Dvector Wake_t::apply_kicks(
    Bunch_t& bun,
    const double stren,
    const double betax,
    ThreadPool& pool) const
{
    my_Dvector Wkick (3,0.0);
    double Wgl(0.0), Wgd(0.0), Wgq(0.0);
    auto& p = bun.particles;
    double&& strenT = stren / betax;
    //Determine the bounds of for loops in each thread
    my_Ivector lims (get_bounds(0,p.size()));

    // After this sorting, the particles will be ordered from head to tail.
    // It means, from smaller ss to bigger ss.
    if (Wd.wake_function || Wq.wake_function || Wl.wake_function ||
        Wd.resonator     || Wq.resonator     || Wl.resonator ){
        bun.general_sort();
    }

    if (Wl.resonator)
        Wgl += apply_wake_resonator_kick(p, LL, stren, lims, pool);
    if (Wd.resonator)
        Wgd += apply_wake_resonator_kick(p, XD, strenT, lims, pool);
    if (Wq.resonator)
        Wgq += apply_wake_resonator_kick(p, XQ, strenT, lims, pool);

    //pw --> witness particle   ps --> source particle
    if (Wd.wake_function || Wq.wake_function || Wl.wake_function){
        #ifdef OPENMP
          #pragma omp parallel for schedule(guided,1) reduction(+:Wgl,Wgd,Wgq)
        #endif
        for (auto w=0;w<p.size();++w){ // Begin from the particle ahead
            if (Wl.wake_function) {
                double&& kick = Wl.WF.get_y(0.0);
                kick *= -stren;
                Wgl     += kick;
                p[w].de += kick;
            }
            for (auto s=w-1;s>=0;--s){ // loop over all particles ahead of it.
                double&& ds = p[w].ss - p[s].ss;
                if (Wl.wake_function) {
                    double&& kick = Wl.WF.get_y(ds);
                    kick *= -stren;
                    Wgl     += kick;
                    p[w].de += kick;
                }
                if (Wd.wake_function) {
                    double&& kick = Wd.WF.get_y(ds);
                    kick *= -p[s].xx * strenT; // The kick is the negative of the wake;
                    Wgd     += kick;
                    p[w].xl += kick;
                }
                if (Wq.wake_function) {
                    double&& kick = Wq.WF.get_y(ds);
                    kick *= -p[w].xx * strenT; // The kick is the negative of the wake;
                    Wgq     += kick;
                    p[w].xl += kick;
                }
            }
        }

        //my parallelization (for now openmp is better)
        // my_Dvector&& Wgs (apply_wake_function_kick(p, stren, strenT, pool));
        // Wgl = Wgs[0]; Wgd = Wgs[1]; Wgq = Wgs[2];
    }

    //pw --> witness particle   ps --> source particle
    if (Wd.wake_potential || Wq.wake_potential || Wl.wake_potential){
      #ifdef OPENMP
        #pragma omp parallel for schedule(guided,1) reduction(+:Wgl,Wgd,Wgq)
      #endif
        for (auto w=0;w<p.size();++w){ // Begin from the particle ahead
            for (auto s=p.size();s>=0;--s){ // loop over all particles.
                double&& ds = p[w].ss - p[s].ss;
                if (Wl.wake_potential) {
                    double&& kick = Wl.WP.get_y(ds);
                    if (abs(kick)> 1e-30){
                        kick *= -stren;
                        Wgl     += kick;
                        p[w].de += kick;
                    }
                }
                if (Wd.wake_potential) {
                    double&& kick = Wd.WP.get_y(ds);
                    if (abs(kick)> 1e-30){
                        kick *= -p[s].xx * strenT; // The kick is the negative of the wake;
                        Wgd     += kick;
                        p[w].xl += kick;
                    }
                }
                if (Wq.wake_potential) {
                    double&& kick = Wq.WP.get_y(ds);
                    if (abs(kick)> 1e-30){
                        kick *= -p[w].xx * strenT; // The kick is the negative of the wake;
                        Wgq     += kick;
                        p[w].xl += kick;
                    }
                }
            }
        }
    }

    Wkick[0] = Wgl;
    Wkick[1] = Wgd;
    Wkick[2] = Wgq;
    return Wkick;
}
