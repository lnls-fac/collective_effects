
#ifndef _ESSENTIALS_H
#define _ESSENTIALS_H

#include <vector>
#include <iostream>
using namespace std;

#define TWOPI  6.28318530717959
#define light_speed 299792458         // [m/s]   - definition

typedef vector<double> my_Dvector;

double interpola(const my_Dvector& si, const my_Dvector& Wi, const double& s);
// {
//     unsigned long i;
//
//     if      (s >= si.back())  {return 0.0;}
//     else if (s <= si.front()) {return 0.0;}
//
//     i = (unsigned long) ((s - si[0])/(si[1]-si[0]));
//     // fprintf(stdout," %05lu",i);
//     return (Wi[i+1] - Wi[i]) / (si[i+1] - si[i]) * s;
// }

#endif
