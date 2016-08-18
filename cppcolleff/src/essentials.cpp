
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
