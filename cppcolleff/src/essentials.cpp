
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

my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector vec2)
{
    my_Dvector conv (vec1.size()+vec2.size()-1,0.0);
    for (int i=0;i<conv.size();++i){
        for (int j=0,k=i;j<vec2.size(),k>=0;++j,--k){
            if (k<vec1.size()){ conv[i] += vec1[k]*vec2[j];}
        }
    }
    return conv;
}

// this function follows matlab's convention of same, not numpy's.
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector vec2)
{
    my_Dvector conv (vec1.size(),0.0);
    int init = vec2.size()/2;
    for (int i=0;i<conv.size();++i){
        for (int j=0,k=(i+init);j<vec2.size(),k>=0;++j,--k){
            if (k<vec1.size()){ conv[i] += vec1[k]*vec2[j];}
        }
    }
    return conv;
}
