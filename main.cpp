#include <iostream>
#include <random>


#include "counter.hpp"
#include "rans.hpp"
#include "timer.hpp"

#include <vector>
#include <cassert>


using u8 = cnt::u8;


template <typename T>
class GeometricDistribution {
private:
    std::random_device rd;
    std::mt19937 gen;
    std::geometric_distribution<T> dist;
public:
    GeometricDistribution(double p): gen(rd()), dist(p) {}
    T operator()() { return dist(gen); }
};


template<class Generator>
auto fill_array_randomly(Generator& g, u8 *v, int n, bool inversion=false) {
    for (int i=0; i<n; ++i) {
        v[i] = (! inversion) ? (g() & 255) : (255 - (g() & 255));
    }
}


static timer_n::Timer timer;

int main() {
    using namespace std;

    std::srand(std::time(0));

    double dt = 0;
    auto calc_perf = [&dt](int n) {
        return (1.e3 * n / dt); // MB/s
    };
    rans::Rans rans;
    int iters = 0;
    constexpr int N_max = 8*1024*1024;
    const int N_performance_is_mearured = N_max / 4;
    std::vector<u8> v; v.reserve(N_max);
    std::vector<u8> output; output.reserve(rans.required_bytes(N_max));
    std::vector<u8> v_decoded; v_decoded.reserve(N_max);
    while (true) {
        size_t N = (((size_t)std::rand()) % N_max) + 1;
        const bool measure_perf = (N >= N_performance_is_mearured);
        std::cout << std::endl;
        std::cout << "New iteration: required output size: " << rans.required_bytes(N) << ", N: " << N << std::endl;
        constexpr int Q = 128;
        const bool inv_f = std::rand() % 2;
        cout << " Test is started, inversion flag: " << inv_f << ", please, wait..." << endl;
        std::vector<double> CRs {};
        std::vector<double> comp_perfs {};
        std::vector<double> decomp_perfs {};
        v.resize(N); v_decoded.resize(N);
        output.resize( rans.required_bytes(N) );
        for (int q=0; q<=Q; ++q) {
            const double par = double(q) / double(Q);
            GeometricDistribution<int> r(par);
            fill_array_randomly(r, v.data(), N, inv_f);
            int out_size;
            timer.reset();
            rans.encode(v.data(), N, output.data(), out_size);
            dt = timer.elapsed_ns();
            const double CR = double(N) / double(out_size);
            CRs.push_back(CR);
            assert( out_size <= (int)output.size() );
            if (measure_perf) {
                // cout << "Compression: ";
                comp_perfs.push_back(calc_perf(N));
            }
            assert(CR >= 1.00000);
            //
            timer.reset();
            rans.decode(output.data(), out_size, v_decoded.data(), N);
            dt = timer.elapsed_ns();
            const bool is_equal = std::equal(v.begin(), v.end(), v_decoded.begin());
            assert( is_equal );
            if (measure_perf) {
                // cout << " Decompression: ";
                decomp_perfs.push_back(calc_perf(N));
            }
        }
        // Get some statistics
        std::sort(CRs.begin(), CRs.end());
        std::sort(comp_perfs.begin(), comp_perfs.end());
        std::sort(decomp_perfs.begin(), decomp_perfs.end());
        const int NN = CRs.size();
        const int NNc = comp_perfs.size();
        assert(comp_perfs.size() == decomp_perfs.size());
        cout << endl;
        cout << "Input length: " << N << " bytes" << endl;
        cout << "Compression Ratio, CR:" << endl;
        cout << " Min: " << CRs[0] << ", Max: " << CRs[NN-1] << ", Median: " << CRs[NN / 2] << endl;
        if (measure_perf) {
            cout << "Decompression performance, MB/s:" << endl;
            cout << " Min: " << decomp_perfs[0] << ", Max: " << decomp_perfs[NNc-1] << ", Median: " << decomp_perfs[NNc / 2] << endl;
            cout << "Compression performance, MB/s:" << endl;
            cout << " Min: " << comp_perfs[0] << ", Max: " << comp_perfs[NNc-1] << ", Median: " << comp_perfs[NNc / 2] << endl;
        } else {
            cout << "Performance measurement is skipped due to small input size: " << N << " that is less than: " << N_performance_is_mearured << endl;
        }
        iters++;
        cout << "Test is Ok! Total iterations: " << iters << endl;
    }

    return 0;
}