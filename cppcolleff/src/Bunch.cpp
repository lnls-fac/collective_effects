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
void Bunch_t::add_offsets(const double xx, const double de)
{
	double xl(0.0), ss(0.0);
	add_offsets(xx, de, xl, ss);
}
void Bunch_t::add_offsets(const double xx)
{
	double xl(0.0), de(0.0), ss(0.0);
	add_offsets(xx, de, xl, ss);
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
		if (plane==XX) k = (p[i].xx - spos[0]-delta) / delta;
		else if (plane==XL) k = (p[i].xl - spos[0]-delta) / delta;
		else if (plane==DE) k = (p[i].de - spos[0]-delta) / delta;
		else if (plane==SS) k = (p[i].ss - spos[0]-delta) / delta;
		if (k < 0) k=0;
		if (k>=spos.size()) k = spos.size()-1;
		distr[k]++;
	}
	for (long&& i=0;i<distr.size();++i){distr[i] /= delta*p.size();}
	return distr;
}
