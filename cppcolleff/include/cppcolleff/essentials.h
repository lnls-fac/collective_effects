
#ifndef _ESSENTIALS_H
#define _ESSENTIALS_H

#include <vector>
#include <iostream>
#include <random>
#include <thread>
#include <cmath>
#include <complex>
#include <omp.h>
using namespace std;

#define TWOPI  6.28318530717959
#define light_speed 299792458.0         // [m/s]   - definition

typedef vector<int> my_Ivector;
typedef vector<double> my_Dvector;
typedef vector<complex<double>> my_Cvector;
typedef vector<default_random_engine> my_RndVector;

class ThreadVars
{
    private:
        static int num_threads;
        static unsigned long seed;
        static void _set_seed_num_threads(const unsigned long s, const int nr)
        {
            num_threads = nr;
            omp_set_num_threads(nr);
            seed = s;
            gens.clear();
            for (int i=0;i<nr;++i) gens.push_back(default_random_engine(s + i));
        }
    public:
        static my_RndVector gens;
        ThreadVars(const bool init)
                            {if (init) set_num_threads(omp_get_num_threads());}
        ~ThreadVars() = default;

        static void set_seed(const unsigned long s)
                                        {_set_seed_num_threads(s, num_threads);}
        static unsigned long get_seed(){return seed;}
        static void set_num_threads(const int nr)
                                            {_set_seed_num_threads(seed, nr);}
        static int get_num_threads(){return num_threads;}

        static my_Ivector get_bounds(const int ini, const int fin);
};

class Particle_t {
public:
    double xx, xl, de, ss;
    Particle_t():xx(0.0),xl(0.0),de(0.0),ss(0.0) {};
    Particle_t(double ini): xx(ini),xl(ini),de(ini),ss(ini) {};
    ~Particle_t() = default;
    // Particle_t& operator += (const Particle_t& b);
    // Particle_t& operator -= (const Particle_t& b);
    // Particle_t& operator *= (const Particle_t& b);
    // Particle_t& operator /= (const Particle_t& b);
    Particle_t& operator += (const Particle_t& b)
                            {xx += b.xx; xl += b.xl; de += b.de; ss += b.ss; return *this;}
    Particle_t& operator -= (const Particle_t& b)
                            {xx -= b.xx; xl -= b.xl; de -= b.de; ss -= b.ss; return *this;}
    Particle_t& operator *= (const Particle_t& b)
                            {xx *= b.xx; xl *= b.xl; de *= b.de; ss *= b.ss; return *this;}
    Particle_t& operator /= (const Particle_t& b)
                            {xx /= b.xx; xl /= b.xl; de /= b.de; ss /= b.ss; return *this;}
    // Particle_t& operator += (const double& b);
    // Particle_t& operator -= (const double& b);
    // Particle_t& operator *= (const double& b);
    // Particle_t& operator /= (const double& b);
    // Particle_t& operator += (const long& b);
    // Particle_t& operator -= (const long& b);
    // Particle_t& operator *= (const long& b);
    // Particle_t& operator /= (const long& b);
    // Particle_t& operator += (const int& b);
    // Particle_t& operator -= (const int& b);
    // Particle_t& operator *= (const int& b);
    // Particle_t& operator /= (const int& b);
    Particle_t& operator += (const double& b)
                            {xx += b; xl += b; de += b; ss += b; return *this;}
    Particle_t& operator -= (const double& b)
                            {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
    Particle_t& operator *= (const double& b)
                            {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
    Particle_t& operator /= (const double& b)
                            {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
    Particle_t& operator += (const long& b)
                            {xx += b; xl += b; de += b; ss += b; return *this;}
    Particle_t& operator -= (const long& b)
                            {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
    Particle_t& operator *= (const long& b)
                            {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
    Particle_t& operator /= (const long& b)
                            {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
    Particle_t& operator += (const int& b)
                            {xx += b; xl += b; de += b; ss += b; return *this;}
    Particle_t& operator -= (const int& b)
                            {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
    Particle_t& operator *= (const int& b)
                            {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
    Particle_t& operator /= (const int& b)
                            {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
    // Particle_t operator + (const Particle_t& b) const;
    // Particle_t operator - (const Particle_t& b) const;
    // Particle_t operator * (const Particle_t& b) const;
    // Particle_t operator / (const Particle_t& b) const;
    // Particle_t operator + (const long& b) const;
    // Particle_t operator - (const long& b) const;
    // Particle_t operator * (const long& b) const;
    // Particle_t operator / (const long& b) const;
    // Particle_t operator + (const int& b) const;
    // Particle_t operator - (const int& b) const;
    // Particle_t operator * (const int& b) const;
    // Particle_t operator / (const int& b) const;
    Particle_t operator + (const Particle_t& b) const
                            {Particle_t c; c += *this; c += b; return c;}
    Particle_t operator - (const Particle_t& b) const
                            {Particle_t c; c -= *this; c -= b; return c;}
    Particle_t operator * (const Particle_t& b) const
                            {Particle_t c; c *= *this; c *= b; return c;}
    Particle_t operator / (const Particle_t& b) const
                            {Particle_t c; c /= *this; c /= b; return c;}
    Particle_t operator + (const long& b) const
                            {Particle_t c; c += *this; c += b; return c;}
    Particle_t operator - (const long& b) const
                            {Particle_t c; c -= *this; c -= b; return c;}
    Particle_t operator * (const long& b) const
                            {Particle_t c; c *= *this; c *= b; return c;}
    Particle_t operator / (const long& b) const
                            {Particle_t c; c /= *this; c /= b; return c;}
    Particle_t operator + (const int& b) const
                            {Particle_t c; c += *this; c += b; return c;}
    Particle_t operator - (const int& b) const
                            {Particle_t c; c -= *this; c -= b; return c;}
    Particle_t operator * (const int& b) const
                            {Particle_t c; c *= *this; c *= b; return c;}
    Particle_t operator / (const int& b) const
                            {Particle_t c; c /= *this; c /= b; return c;}
};

// Particle_t& Particle_t::operator += (const Particle_t& b)
//                         {xx += b.xx; xl += b.xl; de += b.de; ss += b.ss; return *this;}
// Particle_t& Particle_t::operator -= (const Particle_t& b)
//                         {xx -= b.xx; xl -= b.xl; de -= b.de; ss -= b.ss; return *this;}
// Particle_t& Particle_t::operator *= (const Particle_t& b)
//                         {xx *= b.xx; xl *= b.xl; de *= b.de; ss *= b.ss; return *this;}
// Particle_t& Particle_t::operator /= (const Particle_t& b)
//                         {xx /= b.xx; xl /= b.xl; de /= b.de; ss /= b.ss; return *this;}
// Particle_t& Particle_t::operator += (const double& b)
//                         {xx += b; xl += b; de += b; ss += b; return *this;}
// Particle_t& Particle_t::operator -= (const double& b)
//                         {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
// Particle_t& Particle_t::operator *= (const double& b)
//                         {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
// Particle_t& Particle_t::operator /= (const double& b)
//                         {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
// Particle_t& Particle_t::operator += (const long& b)
//                         {xx += b; xl += b; de += b; ss += b; return *this;}
// Particle_t& Particle_t::operator -= (const long& b)
//                         {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
// Particle_t& Particle_t::operator *= (const long& b)
//                         {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
// Particle_t& Particle_t::operator /= (const long& b)
//                         {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
// Particle_t& Particle_t::operator += (const int& b)
//                         {xx += b; xl += b; de += b; ss += b; return *this;}
// Particle_t& Particle_t::operator -= (const int& b)
//                         {xx -= b; xl -= b; de -= b; ss -= b; return *this;}
// Particle_t& Particle_t::operator *= (const int& b)
//                         {xx *= b; xl *= b; de *= b; ss *= b; return *this;}
// Particle_t& Particle_t::operator /= (const int& b)
//                         {xx /= b; xl /= b; de /= b; ss /= b; return *this;}
// Particle_t Particle_t::operator + (const Particle_t& b) const
//                         {Particle_t c; c += *this; c += b; return c;}
// Particle_t Particle_t::operator - (const Particle_t& b) const
//                         {Particle_t c; c -= *this; c -= b; return c;}
// Particle_t Particle_t::operator * (const Particle_t& b) const
//                         {Particle_t c; c *= *this; c *= b; return c;}
// Particle_t Particle_t::operator / (const Particle_t& b) const
//                         {Particle_t c; c /= *this; c /= b; return c;}
// Particle_t Particle_t::operator + (const long& b) const
//                         {Particle_t c; c += *this; c += b; return c;}
// Particle_t Particle_t::operator - (const long& b) const
//                         {Particle_t c; c -= *this; c -= b; return c;}
// Particle_t Particle_t::operator * (const long& b) const
//                         {Particle_t c; c *= *this; c *= b; return c;}
// Particle_t Particle_t::operator / (const long& b) const
//                         {Particle_t c; c /= *this; c /= b; return c;}
// Particle_t Particle_t::operator + (const int& b) const
//                         {Particle_t c; c += *this; c += b; return c;}
// Particle_t Particle_t::operator - (const int& b) const
//                         {Particle_t c; c -= *this; c -= b; return c;}
// Particle_t Particle_t::operator * (const int& b) const
//                         {Particle_t c; c *= *this; c *= b; return c;}
// Particle_t Particle_t::operator / (const int& b) const
//                         {Particle_t c; c /= *this; c /= b; return c;}

Particle_t sqrt(const Particle_t& p);
// {
//     Particle_t c;
//     c.xx = sqrt(p.xx);
//     c.xl = sqrt(p.xl);
//     c.de = sqrt(p.de);
//     c.ss = sqrt(p.ss);
//     return c;
// }

typedef vector<Particle_t> my_PartVector;

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

    inline double get_y(const double& x) const
    {
        if      (x > xi.back())  {return 0.0;}
        else if (x < xi.front()) {return 0.0;}

        unsigned int i;
        if (equally_spaced){i = (unsigned int) (x/dx - offset);}
        else {for (i=0;i<xi.size();++i){if (x<xi[i]){break;}}}

        return dy[i] * x + b[i];
    }
    // inline double get_y(const double& x)const {
    //     if (equally_spaced){
    //         int&& i = (int) (x/dx - offset);
    //
    //         if (i > xi.size()-2)  {return 0.0;}
    //         else if (i < 0) {return 0.0;}
    //
    //         return dy[i] * x + b[i];
    //     } else{
    //         if      (x > xi.back())  {return 0.0;}
    //         else if (x < xi.front()) {return 0.0;}
    //
    //         unsigned int i;
    //         for (i=0;i<xi.size();++i){if (x<xi[i]){break;}}
    //         return dy[i] * x + b[i];
    //     }
    // }
};

my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2);

// this function follows matlab's convention of same, not numpy's:
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2);

#endif
