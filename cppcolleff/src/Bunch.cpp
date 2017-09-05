#include <cppcolleff/Bunch.h>

my_Dvector Bunch_t::get_xx() const
{
	my_Dvector x (particles.size(),0.0);
	for(int i=0;i<x.size();++i) x[i] = particles[i].xx;
	return x;
}
my_Dvector Bunch_t::get_xl() const
{
	my_Dvector x (particles.size(),0.0);
	for(int i=0;i<x.size();++i) x[i] = particles[i].xl;
	return x;
}
my_Dvector Bunch_t::get_de() const
{
	my_Dvector x (particles.size(),0.0);
	for(int i=0;i<x.size();++i) x[i] = particles[i].de;
	return x;
}
my_Dvector Bunch_t::get_ss() const
{
	my_Dvector x (particles.size(),0.0);
	for(int i=0;i<x.size();++i) x[i] = particles[i].ss;
	return x;
}

const my_Ivector& Bunch_t::get_track_indcs() const
{
	return track_indcs;
}

void Bunch_t::set_track_indcs(my_Ivector indcs)
{
	track_indcs.clear();
	for (auto&& i: indcs)
		if (i>=0 && i<particles.size()) track_indcs.push_back(i);
}

inline bool sort_increasing_ss(const Particle_t& first, const Particle_t& second)
{
	if (first.ss < second.ss) return true;
	else return false;
}
void Bunch_t::general_sort()
{
	__gnu_parallel::sort(particles.begin(),particles.end(),sort_increasing_ss);
}
// void Bunch_t::general_sort() {std::sort(particles.begin(),particles.end(),sort_increasing_ss);}

void Bunch_t::insertion_sort()
{
	for (auto p=1;p<particles.size();++p)
		for(auto p2 = p-1; p2>=0;--p2)
			 if (particles[p2+1].ss < particles[p2].ss) {swap(particles[p2+1],particles[p2]);}else break;
}

void Bunch_t::selection_sort()
{
	for (auto p=0;p<particles.size()-1;++p)
		for(auto p2 = p+1; p2<particles.size();++p2)
			 if (particles[p2].ss < particles[p].ss) {swap(particles[p2],particles[p]); break;}
}

void Bunch_t::sort(){
	my_PartVector trv;
	auto& p = particles;
	for (int& i: track_indcs){
		if (i >= 0 && i < p.size()) trv.emplace_back(p[i]);
	}

	if (is_sorted) {insertion_sort();}
	else {general_sort();}
	is_sorted = true;

	for (int i=0; i<track_indcs.size(); ++i){
		int old (track_indcs[i]);
		if (old < 0 || old >= p.size()) continue;
		// dividing the loop in two variables I hope it will be faster:
		for (auto ip=old, in=old-1; (ip<p.size() || in >= 0); ++ip, --in){
			if (ip<p.size() && trv[i] == p[ip]){
				track_indcs[i] = ip;
				break;
			}
			if (in >= 0 && trv[i] == p[in]){
				track_indcs[i] = in;
				break;
			}
		}
	}
}

void Bunch_t::add_particles(const my_PartVector& parts)
{
	num_part += parts.size();
	particles.reserve(num_part);
	for (auto& p:parts) particles.push_back(Particle_t(p));
}

my_PartVector Bunch_t::pick_particles(const int n) const
{
	my_PartVector part;
	default_random_engine gen(seed);
	++seed;
	uniform_int_distribution <int> dist(0,n-1);
	for (int&& i=0; i<n; ++i) part.push_back(Particle_t(particles[dist(gen)]));
	return part;
}

void Bunch_t::add_offsets(const double xx, const double de,
						  const double xl, const double ss)
{
	my_PartVector& p = particles;
	#ifdef OPENMP
	#pragma omp parallel for schedule(guided,1)
	#endif
	for (int i=0;i<p.size();++i){
		p[i].xx += xx;
		p[i].xl += xl;
		p[i].de += de;
		p[i].ss += ss;
	}
}

void Bunch_t::scale_longitudinal(const double scale)
{
	my_PartVector& p = particles;
	#ifdef OPENMP
	#pragma omp parallel for schedule(guided,1)
	#endif
	for (int i=0;i<p.size();++i){
		p[i].de *= scale;
		p[i].ss *= scale;
	}
}

void Bunch_t::scale_transverse(const double scale)
{
	my_PartVector& p = particles;
	#ifdef OPENMP
	#pragma omp parallel for schedule(guided,1)
	#endif
	for (int i=0;i<p.size();++i){
		p[i].xx *= scale;
		p[i].xl *= scale;
	}
}

my_Dvector Bunch_t::calc_dipole_moment(
	const my_Dvector& spos,
	const int plane) const
{
	double ini = spos.front();
	double fin = spos.back();
	const int nbin = spos.size();
	return calc_dipole_moment(ini, fin, nbin, plane);
}

my_Dvector Bunch_t::calc_dipole_moment(
	double ini,
	double fin,
	const int nbin,
	const int plane) const
{
	my_Dvector distr (nbin, 0.0);

	double&& delta = (fin-ini)/(nbin-1);
	double&& offset = -ini/delta + 0.5;

	const my_PartVector& p = particles;
	for (long&& i=0;i<p.size();++i){
		int k;
		if (plane==XX) k = p[i].xx/delta + offset;
		else if (plane==XL) k = p[i].xl/delta + offset;
		else if (plane==DE) k = p[i].de/delta + offset;
		else if (plane==SS) k = p[i].ss/delta + offset;
		if (k < 0) k=0;
		if (k>=nbin) k = nbin-1;
		distr[k]++;
	}
	for (long&& i=0;i<distr.size();++i){distr[i] /= delta*p.size();}
	return distr;
}

my_Dvector Bunch_t::calc_particles_distribution(
	const my_Dvector& spos,
	const int plane) const
{
	double ini = spos.front();
	double fin = spos.back();
	const int nbin = spos.size();
	return calc_particles_distribution(ini, fin, nbin, plane);
}

my_Dvector Bunch_t::calc_particles_distribution(
	double ini,
	double fin,
	const int nbin,
	const int plane) const
{
	my_Dvector distr (nbin, 0.0);

	double&& delta = (fin-ini)/(nbin-1);
	double&& offset = -ini/delta + 0.5;

	const my_PartVector& p = particles;
	for (long&& i=0;i<p.size();++i){
		int k;
		if (plane==XX) k = p[i].xx/delta + offset;
		else if (plane==XL) k = p[i].xl/delta + offset;
		else if (plane==DE) k = p[i].de/delta + offset;
		else if (plane==SS) k = p[i].ss/delta + offset;
		if (k < 0) k=0;
		if (k>=nbin) k = nbin-1;
		distr[k]++;
	}
	for (long&& i=0;i<distr.size();++i){distr[i] /= delta*p.size();}
	return distr;
}

void Bunch_t::to_stream(ostream& fp, const bool isFile) const
{
	fp.setf(fp.left | fp.scientific);
	fp.precision(15);
	fp << setw(30) << "% current" << Ib << " A" << endl;
	fp << setw(30) << "% is_sorted" << (is_sorted ? "true": "false") << endl;
	fp << setw(30) << "% n_particles" << num_part << endl;
	if (!track_indcs.empty()){
		fp << setw(30) << "% track_particles";
		for (auto& i:track_indcs) fp << setw(10) << i;
		fp << endl;
	}
	if (!isFile) return;
	fp << setw(26) << "# xx [m]";
	fp << setw(26) << "xl";
	fp << setw(26) << "de";
	fp << setw(26) << "ss [m]" << endl;
	fp.setf(fp.left | fp.showpos | fp.scientific);
    for (auto& p:particles){
		fp << setw(26) << p.xx;
		fp << setw(26) << p.xl;
		fp << setw(26) << p.de;
		fp << setw(26) << p.ss << endl;
	}
}

void Bunch_t::show_properties() const
{
	ostringstream fp;
	if (fp.fail()) exit(1);
	to_stream(fp, false);
    cout << fp.str();
}

void Bunch_t::to_file(const char* filename) const
{
	ofstream fp(filename);
	if (fp.fail()) exit(1);
	to_stream(fp, true);
    fp.close();
}

void Bunch_t::from_file(const char* filename)
{
	ifstream fp(filename);
	if (fp.fail()) return;

	particles.clear();
	track_indcs.clear();
	double x(0.0), l(0.0), e(0.0), s(0.0);
  	string line;
	unsigned long line_count = 0;
	while (getline(fp, line)) {
  		line_count++;
  		istringstream ss(line);
		char c = ss.get();
		while (c == ' ') c = ss.get();
  		if (c == '#' || c == '\n') continue;
  		if (c == '%') {
			string cmd;
	  		ss >> cmd;
	  		if (cmd.compare("current") == 0) { ss >> Ib; continue; }
	  		if (cmd.compare("n_particles") == 0){
		  		ss >> num_part;
		  		particles.reserve(num_part);
		  		continue;
	  		}
	  		if (cmd.compare("is_sorted") == 0) {
		  		ss >> cmd;
		  		is_sorted = (cmd.compare("true")==0) ? true:false;
		  		continue;
	  		}
			if (cmd.compare("track_particles") == 0){
				int ind;
				while (ss >> ind) track_indcs.push_back(ind);
				continue;
			}
  		}
		ss.unget();
  		ss >> x; ss >> l; ss >> e; ss >> s;
  		particles.push_back(Particle_t(x, l, e, s));
		num_part = particles.size();
	}
	fp.close();
}

void Bunch_t::distribution_to_file(
	const char* filename,
	const double ini,
	const double fin,
	const int nbin,
	const int plane) const
{
	string unit, pl;
	if (plane == XX){unit = "[m]"; pl = "xx";}
	else if (plane == XL){unit = ""; pl = "xl";}
	else if (plane == DE){unit = ""; pl = "de";}
	else if (plane == SS){unit = "[m]"; pl = "ss";}

	my_Dvector&& distr = calc_particles_distribution(ini, fin, nbin, plane);

	save_distribution_to_file(filename, distr, ini, fin, nbin, unit.c_str(), pl.c_str());
}
