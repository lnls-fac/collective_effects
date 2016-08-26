
#ifndef _ESSENTIALS_H
#define _ESSENTIALS_H

#include <vector>
#include <iostream>
#include <thread>
#include <cmath>
using namespace std;

#define TWOPI  6.28318530717959
#define light_speed 299792458.0         // [m/s]   - definition

extern int global_num_threads;

struct Particle_t {
    double xx, xl, de, ss;
    Particle_t():xx(0.0),xl(0.0),de(0.0),ss(0.0) {};
    Particle_t(double ini): xx(ini),xl(ini),de(ini),ss(ini) {};
    ~Particle_t() = default;
};

typedef vector<double> my_Dvector;
typedef vector<Particle_t> my_PartVector;
typedef vector<int> my_Ivector;


class Interpola_t {
private:
    my_Dvector xi, yi;
    bool equally_spaced;
    void check_consistency();
public:
    //Interpola_t(my_Dvector&& xref, my_Dvector&& yref) {xi = move(xref); yi = move(yref);check_consistency();}
    Interpola_t(my_Dvector& xref, my_Dvector& yref): xi(move(xref)), yi(move(yref)) {check_consistency();}
    Interpola_t() = default;
    ~Interpola_t() = default;
    void set_x(my_Dvector& xref){xi = move(xref);check_consistency();}
    void set_y(my_Dvector& yref){yi = move(yref);check_consistency();}
    void set_xy(my_Dvector& xref, my_Dvector& yref){xi = move(xref); yi = move(yref);check_consistency();}
    const my_Dvector& ref_to_xi() const {return xi;}
    const my_Dvector& ref_to_yi() const {return yi;}

    inline double get_y(const double& x) const
    {
        if      (x > xi.back())  {return 0.0;}
        else if (x < xi.front()) {return 0.0;}

        long i;
        if (equally_spaced){i = (long) ((x - xi[0])/(xi[1]-xi[0]));}
        else {for (i=0;i<xi.size();++i){if (x<xi[i]){break;}}}

        return  yi[i]   +   (yi[i+1]-yi[i]) / (xi[i+1]-xi[i]) * (x-xi[i]);
    }
};

my_Ivector bounds_for_threads(const int parts, const int ini, const int final);

my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2);

// this function follows matlab's convention of same, not numpy's:
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2);
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2);

#endif
