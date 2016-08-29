#ifndef _BUNCH_H
#define _BUNCH_H

#include <random>  //std::generator and some distributions
#include <parallel/algorithm> //std::sort
#include <utility> //std::swap
#include <cppcolleff/essentials.h>

class Bunch_t {
  public:
    bool is_sorted;
    long num_part;
    double Ib;  // bunch current;
    my_PartVector particles;
    Bunch_t (const long part): num_part(part), particles(part,0.0) {};
    Bunch_t (const long part,const double I): num_part(part), particles(part,0.0), Ib(I) {};
    ~Bunch_t() = default;

    my_Dvector get_xx() const;
    my_Dvector get_xl() const;
    my_Dvector get_de() const;
    my_Dvector get_ss() const;

    void generate_bunch();
    void general_sort();
    void insertion_sort();
    void selection_sort();
    void sort();
};

#endif
