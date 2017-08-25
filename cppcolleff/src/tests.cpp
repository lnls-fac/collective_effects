#include <numeric>
#include <future>
#include <cppcolleff/essentials.h>
#include <cppcolleff/Bunch.h>
#include <cppcolleff/ThreadPool/ThreadPool.h>

// default_random_engine gen1(1);

int run_thread(int seed, int& maxx)
{
    uniform_int_distribution<int> dist(0,100000);
    default_random_engine gen1(seed);
    for (int i=0;i<maxx; ++i){
        double r = dist(gen1);
    }
    return 1;
}

int main()
{
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > s_;

    std::vector< std::future<int> > results;

    int first = 1000;
    int nr_th = 16;
    int maxx = 100000000/first;
    ThreadPool pool(nr_th);
    std::chrono::time_point<clock_> beg_ = clock_::now();
    for (int j=0;j<first;++j){
            results.emplace_back(pool.enqueue(run_thread, j, ref(maxx)));
    }
    int b = 0;
    for(auto && result: results)
        b +=  result.get();

    cout << "b: " << b << endl;
    cout << "ET: " << chrono::duration_cast<s_> (clock_::now()-beg_).count() << " s" << endl;
    cout << thread::hardware_concurrency() << endl;
}
