#ifndef _BUNCH_H
#define _BUNCH_H

#include <random>  //std::generator and some distributions
#include <algorithm> //std::sort
#include <utility> //std::swap
#include <cppcolleff/essentials.h>

class Bunch_t {
  private:
    bool is_sorted;
  public:
    long num_part;
    double Ib;  // bunch current;
    my_PartVector particles;
    Bunch_t (const long part): num_part(part), particles(part,0.0) {};
    ~Bunch_t() = default;

    void generate_bunch();
    void general_sort();
    void insertion_sort();
    void selection_sort();
    void sort();
};

#endif
