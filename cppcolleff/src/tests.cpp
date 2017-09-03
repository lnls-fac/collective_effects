#include <numeric>
#include <future>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

int main()
{
    set_num_threads(32);
    long N = 50000000;
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    // p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    for (long i=0;i<N;++i){
        in[i][0] = (i>N/2) ? (N-i):i;
        in[i][1] = 2;
    }

    typedef chrono::high_resolution_clock clock_;
    typedef chrono::duration<double, ratio<1> > s_;
    chrono::time_point<clock_> beg_ = clock_::now();

    fftw_execute(p); /* repeat as needed */

    cout << "ET: " << chrono::duration_cast<s_> (clock_::now()-beg_).count() << " s" << endl;

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    return 0;
}
