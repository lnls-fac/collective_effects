#include <numeric>
#include <future>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

int main()
{
    my_Dvector a, b, c, d;
    a.push_back(1);
    a.push_back(2);
    a.push_back(3);
    // a.push_back(2);
    a.push_back(1);

    b.push_back(1);
    b.push_back(1);
    // b.push_back(1);
    // b.push_back(2);
    b.push_back(1);

    Convolve_t A = Convolve_t();
    A.prepare(a, b);

    // c = convolution_full(a, b);
    // d = A.execute();
    // d = convolution_fft(a, b);
    
    c = convolution_same(a, b);
    d = A.execute_same();
    // d = convolution_fft_same(a, b);
    cout << setw(10) << "trad" << endl;
    for (int i=0; i<c.size(); ++i){
        cout << setw(10) << c[i] << endl;
    }
    cout << setw(10) << "fft" << endl;
    for (int i=0; i<d.size(); ++i){
        cout << setw(10) << d[i] << endl;
    }
    return 0;
}
