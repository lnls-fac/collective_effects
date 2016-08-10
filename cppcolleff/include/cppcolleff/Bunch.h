#ifndef _BUNCH_H
#define _BUNCH_H

#include <cppcolleff/essentials.h>

struct Bunch_t {
    long num_part;
    double Ib;  // bunch current;
    my_Dvector de;
    my_Dvector xx;
    my_Dvector xl;
    my_Dvector ss;
    Bunch_t (const long part):
        num_part(part),
        de(part,0.0),   xx(part,0.0),
        xl(part,0.0),   ss(part,0.0){};
    // ~Bunch_t() = default;

    void InsertionSort();
};


#endif
