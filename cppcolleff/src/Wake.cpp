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


void WakePl::to_stream(ostream& fp, const bool isFile) const
{
    auto g_bool = [](bool cond){return (cond) ? "true":"false";};

	fp.setf(fp.left | fp.scientific);
	fp.precision(15);
    fp << setw(30) << "% use_resonator" << g_bool(resonator) << endl;
    fp << setw(30) << "% use_wake_function" << g_bool(wake_function) << endl;
    if (!wr.empty()) {
        fp << "% start_resonators" << endl;
        fp << setw(30) << "% num_resonators" << wr.size() << endl;
        fp << setw(26) << "# Q";
        fp << setw(26) << "wr [rad/s]";
        fp << setw(26) << "Rs [Ohm or Ohm/m]" << endl;
        fp.setf(fp.left | fp.showpos | fp.scientific);
        for (auto i=0; i<wr.size(); ++i){
            fp << setw(26) << Q[i];
            fp << setw(26) << wr[i];
            fp << setw(26) << Rs[i] << endl;
        }
        fp << "% end_resonators" << endl;
    }
    if (!WF.empty()){
        auto& pos = WF.ref_to_xi();
        auto& w = WF.ref_to_yi();
        fp << "% start_wake_function" << endl;
        fp << setw(30) << "% size_wake_function" << pos.size() << endl;
        if (isFile){
            fp << setw(26) << "# pos [m]";
            fp << setw(26) << "W [V/C or V/C/m]" << endl;
            fp.setf(fp.left | fp.showpos | fp.scientific);
            for (auto i=0; i<pos.size(); ++i){
                fp << setw(26) << pos[i];
                fp << setw(26) << w[i] << endl;
            }
        }
        fp << "% end_wake_function" << endl;
    }
}

void WakePl::show_properties() const
{
    if (wr.empty() && WF.empty()) return;

	ostringstream fp;
	if (fp.fail()) exit(1);
	to_stream(fp, false);
    cout << fp.str();
}

void WakePl::to_file(const char* filename) const
{
    if (wr.empty() && WF.empty()) return;

	ofstream fp(filename);
	if (fp.fail()) exit(1);
	to_stream(fp, true);
    fp.close();
}

void WakePl::from_file(const char* filename)
{
	ifstream fp(filename);
	if (fp.fail()) return;

    auto g_bool = [](string& cmd){return (cmd.compare("true")==0) ? true:false;};

    my_Dvector vf1, vf2, vp1, vp2;
	double d1(0.0), d2(0.0), d3(0.0);
    bool fill_re(false), fill_wf(false);
  	string line;
	unsigned long line_count = 0;
	while (getline(fp, line)) {
  		line_count++;
  		istringstream ss(line);
		char c = ss.get();
		while (c == ' ') c = ss.get();
  		if (c == '#' || c == '\n') continue;
  		else if (c == '%') {
			string cmd;
	  		ss >> cmd;
            if (cmd.compare("use_resonator") == 0)
		  		{ss >> cmd; resonator = g_bool(cmd);}
            else if (cmd.compare("use_wake_function") == 0)
		  		{ss >> cmd; wake_function = g_bool(cmd);}
            else if (cmd.compare("start_resonators") == 0)
                {fill_re = true; fill_wf = false;}
            else if (cmd.compare("start_wake_function") == 0)
                {fill_re = false; fill_wf = true;}
            else if (cmd.compare("end_resonators") == 0)
                {fill_re = false;}
            else if (cmd.compare("end_wake_function") == 0)
                {fill_wf = false; WF.set_xy(vf1,vf2);}
            else if (cmd.compare("num_resonators") == 0)
                {int np; ss >> np; Q.reserve(np); wr.reserve(np); Rs.reserve(np);}
            else if (cmd.compare("size_wake_function") == 0)
                {int np; ss >> np; vf1.reserve(np); vf2.reserve(np);}
            continue;
  		}
		ss.unget();
        if (fill_re){
            ss >> d1; ss >> d2; ss >> d3;
            Q.push_back(d1); wr.push_back(d2); Rs.push_back(d3);
        }else if (fill_wf){
            ss >> d1; ss >> d2;
            vf1.push_back(d1); vf2.push_back(d2);
        }
	}
	fp.close();
}


double Wake_t::apply_wake_function_kick(
    Bunch_t& bun,
    double stren,
    int Ktype)
{
    my_PartVector& par = bun.particles;

    WakePl& W = (Ktype == LL) ? Wl: ((Ktype == XD) ? Wd : Wq);

    const my_Dvector& spos = W.WF.ref_to_xi();
    const my_Dvector& Wa = W.WF.ref_to_yi();

    my_Dvector&& distr = (Ktype == XD) ? bun.calc_first_moment(spos, bun.XX):
                                         bun.calc_distribution(spos);
    W.WFC.prepare(distr, Wa, true);
    Interpola_t Kick (spos, W.WFC.execute_same());

    stren *= (spos[1]-spos[0]); // to normalize convolution;
    stren *= bun.num_part; // the value comes divided by the number of particles.
    double Wg (0.0);
    if (Ktype == LL){
        #ifdef OPENMP
          #pragma omp parallel for schedule(guided,1) reduction(+:Wg)
        #endif
        for (auto w=0;w<par.size();++w){
            double&& kick = Kick.get_y(par[w].ss);
            kick *= -stren;
            Wg += kick;
            par[w].de += kick;
        }
    }else if (Ktype == XD){
        #ifdef OPENMP
          #pragma omp parallel for schedule(guided,1) reduction(+:Wg)
        #endif
        for (auto w=0;w<par.size();++w){
            double&& kick = Kick.get_y(par[w].ss);
            kick *= -stren; // The kick is the negative of the wake;
            Wg += kick;
            par[w].xl += kick;
        }
    }else {
        #ifdef OPENMP
          #pragma omp parallel for schedule(guided,1) reduction(+:Wg)
        #endif
        for (auto w=0; w<par.size(); ++w){
            double&& kick = Kick.get_y(par[w].ss);
            kick *= -par[w].xx * stren; // The kick is the negative of the wake;
            Wg += kick;
            par[w].xl += kick;
        }
    }
    return Wg;
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
        for (unsigned int i=0;i<nr_th;++i){
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
    ThreadPool& pool)
{
    //necessary for keeping track of kick received by some particles:
    const my_Ivector& track_indcs = bun.get_track_indcs();
    my_Dvector xl, de;
    for (auto& i: track_indcs){
        xl.push_back(bun.particles[i].xl);
        de.push_back(bun.particles[i].de);
    }
    //

    my_Dvector Wkick (3, 0.0);
    auto& p = bun.particles;
    double&& strenT = stren / betax;
    //Determine the bounds of for loops in each thread
    my_Ivector lims (get_bounds(0,p.size()));

    // After this sorting, the particles will be ordered from head to tail.
    // It means, from smaller ss to bigger ss.
    if (Wd.resonator || Wq.resonator || Wl.resonator) bun.sort();

    if (Wl.resonator)
        Wkick[0] += apply_wake_resonator_kick(p, LL, stren, lims, pool);
    if (Wd.resonator)
        Wkick[1] += apply_wake_resonator_kick(p, XD, strenT, lims, pool);
    if (Wq.resonator)
        Wkick[2] += apply_wake_resonator_kick(p, XQ, strenT, lims, pool);

    if (Wl.wake_function)
        Wkick[0] += apply_wake_function_kick(bun, stren, LL);
    if (Wd.wake_function)
        Wkick[1] += apply_wake_function_kick(bun, strenT, XD);
    if (Wq.wake_function)
        Wkick[2] += apply_wake_function_kick(bun, strenT, XQ);

    Wkick[0] /= bun.num_part;
    Wkick[1] /= bun.num_part;
    Wkick[2] /= bun.num_part;
    for (auto&& i=0; i<track_indcs.size(); ++i){
        Wkick.push_back(bun.particles[track_indcs[i]].de - de[i]);
        Wkick.push_back(bun.particles[track_indcs[i]].xl - xl[i]);
    }
    return Wkick;
}


void Wake_t::show_properties() const
{
	cout << "#########   Longitudinal Wake   #########" << endl;
    Wl.show_properties();
    cout << "#########     Dipolar Wake      #########" << endl;
    Wd.show_properties();
    cout << "#########   Quadrupolar Wake    #########" << endl;
    Wq.show_properties();
}

void Wake_t::to_file(const char* filename) const
{
    string fnameD (filename);
    auto p = fnameD.rfind('.');
    fnameD.insert(p, "_wakeD");
    Wd.to_file(fnameD.c_str());

    string fnameQ (filename);
    fnameQ.insert(p, "_wakeQ");
    Wd.to_file(fnameQ.c_str());

    string fnameL (filename);
    fnameL.insert(p, "_wakeL");
    Wl.to_file(fnameL.c_str());
}

void Wake_t::from_file(const char* filename)
{
    string fnameD (filename);
    auto p = fnameD.rfind('.');
    fnameD.insert(p, "_wakeD");
    Wd.from_file(fnameD.c_str());

    string fnameQ (filename);
    fnameQ.insert(p, "_wakeQ");
    Wd.from_file(fnameQ.c_str());

    string fnameL (filename);
    fnameL.insert(p, "_wakeL");
    Wl.from_file(fnameL.c_str());
}
