
#include <cppcolleff/essentials.h>

void Interpola_t::check_consistency()
{
    if (xi.size() != yi.size()){
        fprintf(stdout,"yi must be the same size as xi.\n");
        exit(1);
    }

    equally_spaced = true;
    double ds0 = xi[1]-xi[0];
    for (auto i=1;i<xi.size(); ++i){
        double&& ds = (xi[i] - xi[i-1]);
        if (ds <= 0.0) {
            fprintf(stdout,"xi must be strictly increasing.\n");
            exit(1);
        }
        if (abs(ds - ds0) > abs(ds0)*1e-10) {equally_spaced = false;}
    }
}

my_Dvector& get_distribution_from_potential(
	const my_Dvector& spos,
	const my_Dvector& V,
	const double& espread,
	const double& scale)
{
    // First I integrate the potential to get the potential well:
    my_Dvector pot_well (spos.size(),0.0);
    double min_ipot = pot_well.back();
    for (int i=spos.size()-2;i>=0;--i){
        pot_well[i] = pot_well[i+1]  +  (V[i] + V[i+1]) * (spos[i+1]-spos[i]) / 2;
        if (pot_well[i] < min_ipot) min_ipot = pot_well[i];
    }

    // Now I compute the distribution
	double&& exp_fact = 1.0  /  (scale * espread*espread);
	double sumdist (0.0);
    for (int i=0;i<spos.size();++i) {
		distr.push_back(exp(  -exp_fact * (pot_well[i]-min_ipot)  ));
		sumdist += distr.back();
	}
	for (int i=0;i<spos.size();++i){distr[i] /= sumdist;} //normalization of distribution;
	return distr;
}

Interpola_t& get_integrated_distribution(const my_Dvector& spos, const my_Dvector& distr) {
	//need to resize the vectors;
	my_Dvector idistr, s_distr;
	idistr.push_back(0.0);
	s_idistr.push_back(spos[0]);
	for (int i=1;i<spos.size();++i){
		double&& idistri = (distr[i]+distr[i-1]) * (spos[i]-spos[i-1]) / 2;
		if (idistri >= 1e-15) { // much more than the number of particles; think about it!
			idistr.push_back(idistr.back() + idistri); // for interpolation to work properly there must
			s_idistr.push_back(spos[i]); // not exist repeated values in the integrated distribution.
		}
	}
	return Interpola_t(idistr,s_idistr);
}

my_Dvector& convolution_same(const my_Dvector& vec1, const my_Dvector vec2)
{
    my_Dvector conv (vec1.size()+vec2.size()-1,0.0);
    for (int k=0;k<conv.size();++k){
        int i1 = k;
        for (int l=0; l<conv.size()+2;  ++l) {
            conv[k] += vec1[l]*vec2[l+vec2.size()-1-k];
        }
    for (i=0; i<conv.size(); i++){
		for (int j=0, k=i; j<vec2.size(), k>=0; ++j, --k)
			if i1<vec1.size())
                conv[i] += vec1[i1]*vec2[j];
	}
    }
}
