#ifndef _WAKE_H
#define _WAKE_H

#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

struct WakePl
{
    bool wake_function, resonator;
    Interpola_t WF; // can be wake functions or wake potentials;
    Convolve_t WFC;
    my_Dvector wr, Rs, Q;
    WakePl(): wake_function(false), resonator(false) {};
    ~WakePl() = default;
    my_Dvector get_wake_at_points(const my_Dvector& spos, const double& stren) const;
    void to_stream(ostream& fp, const bool isFile = true) const;
    void to_file(const char* filename) const;
    void from_file(const char* filename);
    void show_properties() const;
};

class Wake_t
{
    public:
        static const int LL = 0;
        static const int XD = 1;
        static const int XQ = 2;
        WakePl Wd, Wq, Wl;
        Wake_t() {};
        ~Wake_t() = default;

        my_Dvector apply_kicks(
            Bunch_t& bun,
            const double stren,
            const double betax,
            ThreadPool& pool);

        void to_file(const char* filename) const;
        void from_file(const char* filename);
        void show_properties() const;

    private:
        double apply_wake_function_kick(
            Bunch_t& bun,
            double stren,
            int Ktype);

        double apply_wake_resonator_kick(
            my_PartVector& p,
            int Ktype, // 0 for longituinal, 1 dipolar, 2 for quadrupolar
            double stren,
            my_Ivector& lims,
            ThreadPool& pool) const;
};

#endif
