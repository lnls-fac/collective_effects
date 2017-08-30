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

inline bool sort_increasing_ss(const Particle_t& first, const Particle_t& second)
{
	if (first.ss < second.ss) return true;
	else return false;
}
void Bunch_t::general_sort() {__gnu_parallel::sort(particles.begin(),particles.end(),sort_increasing_ss);}
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
	if (is_sorted) {insertion_sort();}
	else {general_sort();}
	is_sorted = true;
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

my_Dvector Bunch_t::calc_particles_distribution(
	const my_Dvector& spos,
	const int plane) const
{
	double delta (spos[1]-spos[0]);
	my_Dvector distr (spos.size(), 0.0);

	const my_PartVector& p = particles;
	for (long&& i=0;i<p.size();++i){
		int k;
		if (plane==XX) k = (p[i].xx - spos[0])/delta + 0.5;
		else if (plane==XL) k = (p[i].xl - spos[0])/delta + 0.5;
		else if (plane==DE) k = (p[i].de - spos[0])/delta + 0.5;
		else if (plane==SS) k = (p[i].ss - spos[0])/delta + 0.5;
		if (k < 0) k=0;
		if (k>=spos.size()) k = spos.size()-1;
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
	if (!isFile) return;
	fp << setw(26) << "# xx [m]";
	fp << setw(26) << "xl [m]";
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
  		}
		ss.unget();
  		ss >> x; ss >> l; ss >> e; ss >> s;
  		particles.push_back(Particle_t(x, l, e, s));
		num_part = particles.size();
	}
	fp.close();
}
