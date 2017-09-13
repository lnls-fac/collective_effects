#include <cppcolleff/essentials.h>

unsigned long seed = 1;
int num_threads = 1;

void set_num_threads(int nr)
{
    num_threads = nr;
    omp_set_num_threads(nr);
    int ret = fftw_init_threads();
    if (ret == 0) {cout << "error" << endl; exit(1);}
    fftw_plan_with_nthreads(nr);
}

int get_num_threads() { return num_threads;}

void set_seed_num(int nr) {seed = nr;}

unsigned long get_seed(){return seed;}

my_Ivector get_bounds(const int ini, const int fin)
{
    return get_bounds(ini, fin, num_threads);
}

my_Ivector get_bounds(const int ini, const int fin, const int nr)
{
    my_Ivector bounds (nr+1,0);
    int incr = (fin-ini) / nr;
    int rem  = (fin-ini) % nr;
    int ad   = 0;
    int N = ini;
    for (auto& bnd:bounds) {
        bnd = N;
        if (rem-- > 0) {ad = 1;} else {ad = 0;}
        N += incr + ad;
    }
    return bounds;
}


void Interpola_t::check_consistency()
{
    if (xi.size() != yi.size()){
        fprintf(stdout,"yi must be the same size as xi.\n");
        exit(1);
    }

    my_Dvector x, y;
    equally_spaced = true;
    double ds0 = xi[1]-xi[0];
    x.push_back(xi[0]);
    y.push_back(yi[0]);
    for (auto i=1;i<xi.size(); ++i){
        double&& ds = (xi[i] - xi[i-1]);
        if (ds <= 1.0e-15) continue;
        x.push_back(xi[i]);
        y.push_back(yi[i]);
        if (abs(ds - ds0) > abs(ds0)*1e-10) {equally_spaced = false;}
    }
    swap(xi, x);
    swap(yi, y);
}

my_Dvector trapz_cumul_integral(const my_Dvector& x, const my_Dvector& y)
{
	my_Dvector iy;
    iy.reserve(y.size());
	iy.push_back(0.0);
	for (int i=1;i<x.size();++i){
		double&& idistri = (y[i]+y[i-1]) * (x[i]-x[i-1]) / 2;
		iy.push_back(iy.back() + idistri);
	}
	return iy;
}

// this function follows matlab's convention of same, not numpy's.
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2)
{
    ThreadPool pool (get_num_threads());
    return convolution_same_orig(vec1, vec2, pool);
}
my_Dvector convolution_same_orig(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool)
{
    my_Dvector conv (vec1.size(),0.0);
    int init = vec2.size()/2;

    unsigned int nr_th = get_num_threads();
    my_Ivector lims (get_bounds(0,conv.size()));
    std::vector< std::future<int> > results;

    //lambda function
    auto func = [&](my_Dvector& v, int com, int ter)
    {
        for (int ii=com;ii<ter;++ii){
            for (int j=0,k=(ii+init);j<vec2.size() && k>=0;++j,--k){
                if (k<vec1.size()){ v[ii] += vec1[k]*vec2[j];}
            }
        }
        return 1;
    }; //
    for (unsigned int i=0;i<nr_th;++i) results.emplace_back(pool.enqueue(
                        func, ref(conv), lims[i], lims[i+1]));
    for(auto && result: results) result.get();

    return conv;
}
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2)
{
    ThreadPool pool (get_num_threads());
    return convolution_full_orig(vec1, vec2, pool);
}
my_Dvector convolution_full_orig(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool)
{
    my_Dvector conv (vec1.size()+vec2.size()-1,0.0);

    unsigned int nr_th = get_num_threads();
	my_Ivector lims (get_bounds(0,conv.size()));
	std::vector< std::future<int> > results;

    //lambda function
    auto func = [&](my_Dvector& v, int com, int ter)
    {
        for (int ii=com;ii<ter;++ii){
            for (int j=0,k=ii;j<vec2.size() && k>=0;++j,--k){
                if (k<vec1.size()){ v[ii] += vec1[k]*vec2[j];}
            }
        }
        return 1;
    }; //
	for (unsigned int i=0;i<nr_th;++i) results.emplace_back(pool.enqueue(
                        func, ref(conv), lims[i], lims[i+1]));
    for(auto && result: results) result.get();

    return conv;
}


/* The functions below implement the convolutions just like above, but with an
optimization for the case when the vectors to be convoluted have zeros at the
beginning and end. The speedup is of the order o 16 times for the convolutions
between the wake function and the bunch distribution.*/
static void _trim_vec(const my_Dvector& vec, int& mini, int& maxi)
{
    bool minb (true), maxb (true);
    for (int i=0,j=vec.size(); i<vec.size(); ++i,--j){
        if (minb && (vec[i]  >1e-30 || vec[i]  <-1e-30)) {minb = false; mini = i;}
        if (maxb && (vec[j-1]>1e-30 || vec[j-1]<-1e-30)) {maxb = false; maxi = j;}
        if ((!minb) && (!maxb)) {break;}
    }
}

static void _convolution(const my_Dvector& vec1, const my_Dvector& vec2, my_Dvector& conv, int init, ThreadPool& pool)
{
    // To do the convolution I trim the beginning and ending points if they are zero (1e-30);
    int min1 (0), max1 (vec1.size());
    _trim_vec(vec1, min1, max1);

    int min2 (0), max2 (vec2.size());
    _trim_vec(vec2, min2, max2);

    int ini = max(min1+min2-init,0);
    int fin = min(int(conv.size()), max1 + max2 - 1);

    unsigned int nr_th = get_num_threads();
    unsigned int nr_job = max((fin-ini)/100, min(int(nr_th), fin-ini));
	my_Ivector lims (get_bounds(ini, fin, nr_job));
	std::vector< std::future<int> > results;

    //lambda function
    auto func = [&](my_Dvector& v, int com, int ter)
    {
        for (int ii=com;ii<ter;++ii){
            const int delta = max(ii+init-min2-max1,0);
            for (int j=min2+delta,k=(ii+init-min2-delta);j<max2 && k>=min1;++j,--k){
                v[ii] += vec1[k]*vec2[j];
            }
        }
        return 1;
    }; //

	for (unsigned int i=0;i<nr_job;++i) results.emplace_back(pool.enqueue(
        func, ref(conv), lims[i], lims[i+1]));
    for(auto && result: results) result.get();


}
// this function follows matlab's convention of same, not numpy's.
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2)
{
    ThreadPool pool (get_num_threads());
    return convolution_same(vec1, vec2, pool);
}
my_Dvector convolution_same(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool)
{
    my_Dvector conv (vec1.size(),0.0);
    int init = vec2.size()/2;

    _convolution(vec1, vec2, conv, init, pool);
    return conv;
}
my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2)
{
    ThreadPool pool (get_num_threads());
    return convolution_full(vec1, vec2, pool);
}
my_Dvector convolution_full(const my_Dvector& vec1, const my_Dvector& vec2, ThreadPool& pool)
{
    my_Dvector conv (vec1.size()+vec2.size()-1,0.0);
    int init (0);

    _convolution(vec1, vec2, conv, init, pool);
    return conv;
}

void Convolve_t::free_memory()
{
    if (in1 == nullptr) return;
    fftw_free(in1); in1 = nullptr;
    fftw_free(in2); in2 = nullptr;
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(pr);
}
void Convolve_t::create_plans(const long n1, const long n2, const bool meas)
{
    if (n1 != N1 || n2 != N2) free_memory();
    if (in1 != nullptr) return;
    N1 = n1;
    N2 = n2;
    N = N1 + N2;
    measure = meas;
    auto flag = measure ? FFTW_MEASURE: FFTW_ESTIMATE;
    in1 = fftw_alloc_real(N);
    in2 = fftw_alloc_real(N);
    p1 = fftw_plan_r2r_1d(N, in1, in1, FFTW_R2HC, flag);
    p2 = fftw_plan_r2r_1d(N, in2, in2, FFTW_R2HC, flag);
    pr = fftw_plan_r2r_1d(N, in1, in1, FFTW_HC2R, flag);
}

void Convolve_t::prepare(
    const my_Dvector& vec1,
    const my_Dvector& vec2,
    const bool meas)
{
    create_plans(vec1.size(), vec2.size(), meas);
    copy(vec1.begin(), vec1.end(), in1);
    fill(in1+N1, in1+N, 0.0);
    copy(vec2.begin(), vec2.end(), in2);
    fill(in2+N2, in2+N, 0.0);
}

my_Dvector Convolve_t::execute()
{
    fftw_execute(p1);
    fftw_execute(p2);

    in1[0] *= in2[0];
    long stop = N/2 + 1;
    if (N%2 == 0){
        in1[N/2] *= in2[N/2];
        stop = N/2;
    }
    for (long i=1, j=N-1; i<stop; ++i, --j){
        double&& tmp = in1[i] * in2[i];
        tmp -= in1[j]*in2[j];
        in1[j] *= in2[i];
        in1[j] += in1[i]*in2[j];
        in1[i] = tmp;
    }

    fftw_execute(pr);

    my_Dvector res;
    res.reserve(N-1);
    for (long i=0; i<N-1; ++i) res.push_back(in1[i]/N);
    return res;
    }

my_Dvector Convolve_t::execute_same()
{
    my_Dvector&& res = execute();
    my_Dvector res2;
    res2.resize(N1);

    long&& ini = N2/2;
    long&& fin = ini + N1;
    copy(res.begin()+ini, res.begin()+fin, res2.begin());
    return res2;
}

my_Dvector convolution_fft(const my_Dvector& vec1, const my_Dvector& vec2)
{
    Convolve_t conv(vec1.size(), vec2.size());
    conv.prepare(vec1, vec2);
    return conv.execute();
}

my_Dvector convolution_fft_same(const my_Dvector& vec1, const my_Dvector& vec2)
{
    Convolve_t conv(vec1.size(), vec2.size());
    conv.prepare(vec1, vec2);
    return conv.execute_same();
}

void save_distribution_to_file(
    const char* filename,
    const my_Dvector& distr,
    const double ini,
    const double fin,
    const int nbin,
    const char* unit,
    const char* pl)
{
    ofstream fp(filename);
    if (fp.fail()) exit(1);
    fp.setf(fp.left | fp.scientific);
    fp.precision(15);
    fp << setw(30) << "% initial" << ini << " " << unit << endl;
    fp << setw(30) << "% final" << fin  << " " << unit << endl;
    fp << setw(30) << "% nbins" << nbin << endl;
    fp << "# " << pl << " " << unit << endl;
    fp.setf(fp.left | fp.showpos | fp.scientific);
    for (auto& p:distr) fp << p << endl;
    fp.close();
}

Interpola_t load_distribution_from_file(const char* filename)
{
    ifstream fp(filename);
	if (fp.fail()) return Interpola_t();
    my_Dvector distr;

    int nbin;
  	string line;
    double ini, fim;
	unsigned long line_count = 0;
	while (getline(fp, line)) {
  		line_count++;
  		istringstream ss(line);
		char c = ss.get();
		while (c == ' ') c = ss.get();
  		if (c == '#' || c == '\n') continue;
  		if (c == '%') {
			string cmd;
	  		ss >> cmd;
	  		if (cmd.compare("initial") == 0) {ss >> ini;}
	  		else if (cmd.compare("final") == 0){ss >> fim;}
            else if (cmd.compare("nbins") == 0){ss >> nbin; distr.reserve(nbin);}
            continue;
  		}
		ss.unget();
        double x;
  		ss >> x;
        distr.push_back(x);
	}
	fp.close();
    my_Dvector s(nbin, 0.0);
    double ds((fim-ini)/(nbin-1));
    s[0] = ini;
    for (int i=1;i<distr.size(); ++i) {s[i] = s[i-1]; s[i] += ds;}
    return Interpola_t(s, distr);
}
