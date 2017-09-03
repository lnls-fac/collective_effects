
#ifndef _ESSENTIALS_H
#define _ESSENTIALS_H

// IO
#include <iostream> // std::cout, std:cerr, std::cin
#include <ostream> // ostream class
#include <iomanip> //std::setw
#include <fstream> // ifstream and ofstream classes
#include <sstream> // istringstream and ostringstream classes
#include <cstdio> //std::sprintf std::fprintf
#include <string> //std::string.c_str()
//Utilities:
#include <vector>
#include <utility> //std::swap
#include <random>  //std::generator and some distributions
#include <cmath> //std::sqrt std::sin std::cos
#include <complex>
//parallelization
#include <thread>
#include <parallel/algorithm> //std::sort
#include <omp.h>
#include <cppcolleff/ThreadPool/ThreadPool.h>
//extern libraries for special calculations
#include <fftw3.h> // for fft and convolution

using namespace std;

#define TWOPI  6.28318530717959
#define light_speed 299792458.0         // [m/s]   - definition

typedef vector<int> my_Ivector;
typedef vector<double> my_Dvector;
typedef vector<my_Dvector> my_Dmatrix;
typedef vector<complex<double>> my_Cvector;

extern unsigned long seed;
extern int num_threads;

void set_num_threads(int nr);
int get_num_threads();
void set_seed_num(int nr);
unsigned long get_seed();
my_Ivector get_bounds(const int ini, const int fin);
my_Ivector get_bounds(const int ini, const int fin, const int nr);


class Particle_t {
public:
    double xx, xl, de, ss;
    Particle_t(const Particle_t& p):xx(p.xx), xl(p.xl), de(p.de), ss(p.ss){};
    Particle_t(
        const double x = 0.0,
        const double l = 0.0,
        const double e = 0.0,
        const double s = 0.0): xx(x),xl(l),de(e),ss(s) {};
    Particle_t(const long ini): xx(ini),xl(ini),de(ini),ss(ini) {};
    ~Particle_t() = default;
    bool operator ==(const Particle_t& b)
    {
        return (fabs(xx - b.xx) < 1e-13 && fabs(xl - b.xl) < 1e-13 &&
                fabs(de - b.de) < 1e-10 && fabs(ss - b.ss) < 1e-10);
    }
    bool operator != (const Particle_t& b) {return !(*this==b);}
    Particle_t& operator += (const Particle_t& b)
                {xx += b.xx; xl += b.xl; de += b.de; ss += b.ss; return *this;}
    Particle_t& operator -= (const Particle_t& b)
                {xx -= b.xx; xl -= b.xl; de -= b.de; ss -= b.ss; return *this;}
    Particle_t& operator *= (const Particle_t& b)
                {xx *= b.xx; xl *= b.xl; de *= b.de; ss *= b.ss; return *this;}
    Particle_t& operator /= (const Particle_t& b)
                {xx /= b.xx; xl /= b.xl; de /= b.de; ss /= b.ss; return *this;}
    template <class T> Particle_t& operator += (const T& b)
                            {xx += b; xl += b; de += b; ss += b; return *this;}
    template <class T> Particle_t& operator -= (const T& b)
                            {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
    template <class T> Particle_t& operator *= (const T& b)
                            {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
    template <class T> Particle_t& operator /= (const T& b)
                            {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
    template <class T> Particle_t operator + (const T& b) const
                            {Particle_t c = *this; c += b; return c;}
    template <class T> Particle_t operator - (const T& b) const
                            {Particle_t c = *this; c -= b; return c;}
    template <class T> Particle_t operator * (const T& b) const
                            {Particle_t c = *this; c *= b; return c;}
    template <class T> Particle_t operator / (const T& b) const
                            {Particle_t c = *this; c /= b; return c;}
};

inline Particle_t sqrt(const Particle_t& p)
{
    Particle_t c;
    c.xx = sqrt(p.xx);
    c.xl = sqrt(p.xl);
    c.de = sqrt(p.de);
    c.ss = sqrt(p.ss);
    return c;
}

typedef vector<Particle_t> my_PartVector;
typedef vector<my_PartVector> my_PartMatrix;

class Interpola_t {
private:
    my_Dvector xi, yi;
    my_Dvector dy, b; // pre-calculation for interpolation;
    double dx, offset;
    bool equally_spaced;
    void check_consistency();
    void initialize_interp() {
        dx = xi[1] - xi[0];
        offset = xi[0]/dx;
        for (int i=0;i<xi.size()-1;++i){
            dy.push_back(  (yi[i+1]-yi[i]) / (xi[i+1]-xi[i])  );
            b.push_back(  yi[i]   -   dy.back() * xi[i]   );
        }
    }
public:
    //Interpola_t(my_Dvector&& xref, my_Dvector&& yref) {xi = move(xref); yi = move(yref);check_consistency();}
    Interpola_t(my_Dvector& xref, my_Dvector& yref):
    xi(xref), yi(yref) {check_consistency();initialize_interp();}
    Interpola_t() = default;
    ~Interpola_t() = default;
    void set_x(my_Dvector& xref){xi=xref; check_consistency(); initialize_interp();}
    void set_y(my_Dvector& yref){yi=yref; check_consistency(); initialize_interp();}
    void set_xy(my_Dvector& xref, my_Dvector& yref){xi=xref; yi=yref; check_consistency(); initialize_interp();}
    const my_Dvector& ref_to_xi() const {return xi;}
    const my_Dvector& ref_to_yi() const {return yi;}
    bool empty() const {return xi.empty();}

    inline double get_y(const double& x) const
    {
        if      (x > xi.back())  {return 0.0;}
        else if (x < xi.front()) {return 0.0;}

        unsigned int i;
        if (equally_spaced){i = (unsigned int) (x/dx - offset);}
        else {for (i=0;i<xi.size();++i){if (x<xi[i]){break;}}}

        return dy[i] * x + b[i];
    }
};

my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool);
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool);

// this function follows matlab's convention of same, not numpy's:
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool);
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool);


void save_distribution_to_file(
    const char* filename,
    const my_Dvector& distr,
    const double ini,
    const double fin,
    const int nbin,
    const char* unit = "[m]",
    const char* pl = "ss");

my_Dvector load_distribution_from_file(const char* filename, my_Dvector& lims);

#endif
