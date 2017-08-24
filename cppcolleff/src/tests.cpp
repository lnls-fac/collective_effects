#include <numeric>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>

uniform_int_distribution<int> dist(0,100000);


void run_thread()
{
    default_random_engine gen1(1);
    for (int i=0;i<10; ++i){
        cout << dist(gen1) << endl;
}

int main()
{
    default_random_engine gen1(1);
    vector<thread> ths;


    for (int i=0;i<3; ++i){
        cout << dist(gen1) << endl;
    }
}
