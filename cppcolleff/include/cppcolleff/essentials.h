
#ifndef _ESSENTIALS_H
#define _ESSENTIALS_H

#include <vector>
#include <iostream>
#include <thread>
#include <cmath>
#include <complex>
#include <omp.h>
using namespace std;

#define TWOPI  6.28318530717959
#define light_speed 299792458.0         // [m/s]   - definition

class NumThreads {
  private:
    static int num_threads;
  public:
    static void set_num_threads(int nr) {num_threads = nr; omp_set_num_threads(nr);}
    static int get_num_threads(){return num_threads;}
};

struct Particle_t {
    double xx, xl, de, ss;
    Particle_t():xx(0.0),xl(0.0),de(0.0),ss(0.0) {};
    Particle_t(double ini): xx(ini),xl(ini),de(ini),ss(ini) {};
    ~Particle_t() = default;
};

typedef vector<int> my_Ivector;
typedef vector<double> my_Dvector;
typedef vector<complex<double>> my_Cvector;
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

my_Ivector bounds_for_threads(const int parts, const int ini, const int final);

my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2);

// this function follows matlab's convention of same, not numpy's:
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2);

#endif
