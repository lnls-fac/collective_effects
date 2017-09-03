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

    equally_spaced = true;
    double ds0 = xi[1]-xi[0];
    for (auto i=1;i<xi.size(); ++i){
        double&& ds = (xi[i] - xi[i-1]);
        if (ds <= 0.0) {
            fprintf(stdout,"xi must be strictly increasing.\n");
            exit(1);
        }
        if (abs(ds - ds0) > abs(ds0)*1e-10) {equally_spaced = false;}
    }
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

my_Dvector convolution_fft(const my_Dvector& vec1 const my_Dvector& vec2)
    set_num_threads(32);
    fftw_complex *in1, in2;
    fftw_plan p;

    long N = vec1.size();
    in1 = (double*) fftw_malloc(sizeof(double) * 2*(N/2+1));
    in2 = (double*) fftw_malloc(sizeof(double) * 2*(N/2+1));
    inr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
    // p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p1 = fftw_plan_dft_r2c_1d(N, in1, in1, FFTW_MEASURE);
    p2 = fftw_plan_dft_r2c_1d(N, in2, in2, FFTW_MEASURE);
    pr = fftw_plan_dft_c2r_1d(N, in1, in1, FFTW_MEASURE);

    p1 = fftw_plan_r2r_1d(N, in1, in1, FFTW_R2HC, FFTW_MEASURE);
    p2 = fftw_plan_r2r_1d(N, in2, in2, FFTW_R2HC, FFTW_MEASURE);
    p1 = fftw_plan_r2r_1d(N, in1, in1, FFTW_HC2R, FFTW_MEASURE);


    for (long i=0;i<N;++i){
        in1[i][0] = vec1[i];
        in2[i][0] = vec2[i];
    }

    fftw_execute(p1); /* repeat as needed */
    fftw_execute(p2); /* repeat as needed */

    for (long k=0, i=0, j=1; k<(N/2+1); ++k, i+=2, j+=2){
        inr[k][0] = in1[i]*in2[i] - in2[j]*in2[j];
        inr[k][1] = in1[i]*in2[j] + in1[j]*in2[i];
    }

    fftw_execute(pr);

    my_Dvector res;
    res.reserve(N);
    for (long i=0;i<N; i++) res.push_back(in1[i][0]/N)

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return res;
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

my_Dvector load_distribution_from_file(const char* filename, my_Dvector& lims)
{
    my_Dvector distr;
    ifstream fp(filename);
	if (fp.fail()) return distr;

    int nbin;
  	string line;
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
	  		if (cmd.compare("initial") == 0) {ss >> lims[0];}
	  		else if (cmd.compare("final") == 0){ss >> lims[1];}
            else if (cmd.compare("nbins") == 0){ss >> nbin; distr.reserve(nbin);}
            continue;
  		}
		ss.unget();
        double x;
  		ss >> x;
        distr.push_back(x);
	}
	fp.close();
    return distr;
}
