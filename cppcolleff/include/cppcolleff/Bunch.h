#ifndef _BUNCH_H
#define _BUNCH_H

#include <list>
#include <random>
#include <algorithm>
#include <cppcolleff/essentials.h>

struct Bunch_t {
    long num_part;
    double Ib;  // bunch current;
    my_PartVector particles;
    Bunch_t (const long part): num_part(part), particles(part,0.0) {};
    ~Bunch_t() = default;

    void generate_bunch();
    void sort();
};

#endif
