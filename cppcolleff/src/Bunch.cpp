
#include <cppcolleff/Bunch.h>


static bool sort_increasing_ss(const Particle_t& first, const Particle_t& second)
{
	if (first.ss < second.ss) return true;
	else return false;
}
void Bunch_t::sort() {std::sort(particles.begin(),particles.end(),sort_increasing_ss);}

// void Bunch_t::sort()
// {
// 	for (long i=1; i < xx.size(); i++){
// 		double atualxx=xx[i]; double atualxl=xl[i]; double atualde=de[i]; double atualss=ss[i];
//         long j = i - 1;
// 		while ((j>=0) && (atualss < ss[j])){
// 			ss[j+1]=ss[j];   xx[j+1]=xx[j];    xl[j+1]=xl[j];    de[j+1]=de[j];
//             --j;
// 		}
// 		ss[j+1] = atualss;    xx[j+1] = atualxx;    xl[j+1] = atualxl;    de[j+1] = atualde;
// 	}
// }

void Bunch_t::generate_bunch()
{
	default_random_engine generator(19880419);
  	normal_distribution<double> distribution(0.0,1.0);

	for (auto p=particles.begin();p!=particles.end();++p){
		(*p).xx = distribution(generator);
		(*p).xl = distribution(generator);
		(*p).de = distribution(generator);
		(*p).ss = distribution(generator);
	}
}
