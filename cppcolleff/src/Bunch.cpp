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

void Bunch_t::generate_bunch()
{
	default_random_engine generator(19880419);
  	normal_distribution<double> distribution(0.0,1.0e-3);

	for (auto& p:particles){
		p.xx = distribution(generator);
		p.xl = distribution(generator);
		p.de = distribution(generator);
		p.ss = distribution(generator);
	}
	is_sorted = false;
}
