#include <cppcolleff/essentials.h>

static void ThreadVars::_set_seed_num_threads(
    const unsigned long s,
    const int nr)
{
    num_threads = nr;
    omp_set_num_threads(nr);
    seed = s;
    gens.clear();
    for (int i=0; i<nr; ++i){
        default_random_engine gen (s + i);
        gens.push_back(gen);
    }
}

my_Ivector ThreadVars::get_bounds(const int ini, const int final)
{
    my_Ivector bounds (num_threads+1,0);
    int incr = (final-ini) / num_threads;
    int rem  = (final-ini) % num_threads;
    int ad   = 0;
    int N = ini;
    for (auto& bnd:bounds) {
        bnd = N;
        if (rem-- > 0) {ad = 1;} else {ad = 0;}
        N += incr + ad;
    }
    return bounds;
}

static unsigned int _tmp = omp_get_num_threads();
ThreadVars ThreadInfo (_tmp);

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

// this function follows matlab's convention of same, not numpy's.
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2)
{
    my_Dvector conv (vec1.size(),0.0);
    int init = vec2.size()/2;
  #ifdef OPENMP
    #pragma omp parallel for schedule(guided,1)
  #endif
    for (int i=0;i<conv.size();++i){
        for (int j=0,k=(i+init);j<vec2.size() && k>=0;++j,--k){
            if (k<vec1.size()){ conv[i] += vec1[k]*vec2[j];}
        }
    }
    return conv;
}
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2)
{
    my_Dvector conv (vec1.size()+vec2.size()-1,0.0);
  #ifdef OPENMP
    #pragma omp parallel for schedule(guided,1)
  #endif
    for (int i=0;i<conv.size();++i){
        for (int j=0,k=i;j<vec2.size() && k>=0;++j,--k){
            if (k<vec1.size()){ conv[i] += vec1[k]*vec2[j];}
        }
    }
    return conv;
}


/* The functions below implement the convolutions just like above, but with an
optimization for the case when the vectors to be convoluted have zeros at the
beginning and end. The speedup is of the order o 16 times for the convolutions
between the wake function and the bunch distribution.*/
static void _trim_vec(const my_Dvector& vec, int& mini, int& maxi)
{
    bool minb (true), maxb (true);
    for (int i=0,j=vec.size(); i<vec.size(); ++i,--j){
        if (minb && (vec[i]  >1e-30 || vec[i]  <-1e-30)) {minb = false; mini = i;}
        if (maxb && (vec[j-1]>1e-30 || vec[j-1]<-1e-30)) {maxb = false; maxi = j;}
        if ((!minb) && (!maxb)) {break;}
    }
}
static void _convolution(const my_Dvector& vec1, const my_Dvector& vec2, my_Dvector& conv, int init)
{
    // To do the convolution I trim the beginning and ending points if they are zero (1e-30);
    int min1 (0), max1 (vec1.size());
    _trim_vec(vec1, min1, max1);

    int min2 (0), max2 (vec2.size());
    _trim_vec(vec2, min2, max2);

    int ini = max(min1+min2-init,0);
    int fin = min(int(conv.size()), max1 + max2 - 1);

  #ifdef OPENMP
    #pragma omp parallel for schedule(guided,1)
  #endif
    for (int i=ini;i<fin;++i){
        const int delta = max(i+init-min2-max1,0);
        for (int j=min2+delta,k=(i+init-min2-delta);j<max2 && k>=min1;++j,--k){
            conv[i] += vec1[k]*vec2[j];
        }
    }
}
// this function follows matlab's convention of same, not numpy's.
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2)
{
    my_Dvector conv (vec1.size(),0.0);
    int init = vec2.size()/2;

    _convolution(vec1, vec2, conv, init);
    return conv;
}
my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2)
{
    my_Dvector conv (vec1.size()+vec2.size()-1,0.0);
    int init (0);

    _convolution(vec1, vec2, conv, init);
    return conv;
}
