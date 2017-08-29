#ifndef _BUNCH_H
#define _BUNCH_H

#include <cppcolleff/essentials.h>

class Bunch_t {
    public:
        static const int XX = 0;
        static const int XL = 1;
        static const int DE = 2;
        static const int SS = 3;
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

        void general_sort();
        void insertion_sort();
        void selection_sort();
        void sort();

        void add_offsets(
            const double xx,
            const double de = 0.0,
        	const double xl = 0.0,
            const double ss = 0.0);

        void scale_longitudinal(const double scale);
        void scale_transverse(const double scale);

        my_Dvector calc_particles_distribution(
            const my_Dvector& spos,
            const int plane = SS) const;

        void to_file(const char* filename) const;
        void from_file(const char* filename);
};

#endif
