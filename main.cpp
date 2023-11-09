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
    double current_perf;
    auto disp_perf = [&dt, &current_perf](int n, bool display=true) {
        current_perf = (1.e3 * n / dt);
        if (display) {
            std::cout << " performance: " << current_perf << " MB/s" << std::endl;
        }
    };

    rans::Rans rans;
    // cnt::Counter<65536> c;
    int iters = 0;
    constexpr int N_max = 8*1024*1024;
    const int N_performance_is_mearured = N_max / 4;
    while (true) {
        size_t N = (((size_t)std::rand()) % N_max) + 1;
        // size_t N = 32;
        std::vector<u8> v(N);
        std::vector<u8> output(rans.required_bytes(N));

        const bool measure_perf = (N >= N_performance_is_mearured);

        std::cout << std::endl;
        std::cout << "New iteration: required output size: " << output.size() << ", N: " << N << std::endl;

        std::vector<u8> v_decoded(N);
        constexpr int Q = 128;
        std::vector<double> perfomances(Q);
        bool inv_f = std::rand() % 2;

        cout << " Test is started, inversion flag: " << inv_f << ", please, wait..." << endl;
        std::vector<double> CRs {};
        double min_decomp_perf = 1.e18;
        double max_decomp_perf = 0;
        double min_comp_perf = 1.e18;
        double max_comp_perf = 0;
        double min_CR = 1.e18;
        double max_CR = 0;
        for (int q=0; q<=Q; ++q) {
            const double par = double(q) / double(Q);
            GeometricDistribution<int> r(par);
            fill_array_randomly(r, v.data(), N, inv_f);
            int out_size;
            timer.reset();
            rans.encode(v.data(), N, output.data(), out_size);
            // auto f { c.count(v.data(), N) };
            dt = timer.elapsed_ns();
            // cnt::u64 sum_f = 0;
            // for (auto el : f) {
                // sum_f += el;
            // }
            // assert(sum_f == c.get_L());
            // cout << " geometric distr. p: " << par << ", sum of frequencies: " << sum_f << ", ";
            const double CR = double(N) / double(out_size);
            CRs.push_back(CR);
            min_CR = std::min(min_CR, CR);
            max_CR = std::max(max_CR, CR);
            // cout << " geom. distr. p: " << par << ", compressed size: " << out_size << ", bytes. CR: " << CR << ", ";
            if (measure_perf) {
                // cout << "Geometric distribution p: " << par << ", CR: " << CR << endl;
                // cout << "Compression: ";
                disp_perf(N, false);
                min_comp_perf = std::min(min_comp_perf, current_perf);
                max_comp_perf = std::max(max_comp_perf, current_perf);
            }
            assert(CR >= 1.00000);
            //
            timer.reset();
            rans.decode(output.data(), out_size, v_decoded.data(), N);
            dt = timer.elapsed_ns();
            // cout << "  decompressed size: " << v_decoded.size() << ", bytes" << ", ";
            const bool is_equal = std::equal(v.begin(), v.end(), v_decoded.begin());
            // if (! is_equal) {
            //     std::cout << "Original vector: ";
            //     for (const auto& el : v) {
            //         std::cout << int(el) << ", ";
            //     }
            //     std::cout << std::endl;
            //     std::cout << "Decompressed vector: ";
            //     for (const auto& el : v_decoded) {
            //         std::cout << int(el) << ", ";
            //     }
            //     std::cout << std::endl;
            // }
            assert( is_equal );
            if (measure_perf) {
                // cout << " Decompression: ";
                disp_perf(N, false);
                min_decomp_perf = std::min(min_decomp_perf, current_perf);
                max_decomp_perf = std::max(max_decomp_perf, current_perf);
            }
        }
        std::sort(CRs.begin(), CRs.end());
        const double CR_median = CRs[CRs.size()/2];
        cout << endl;
        cout << "Input length: " << N << " bytes" << endl;
        cout << "Median CR: " << CR_median << endl;
        cout << "Min CR: " << min_CR << endl;
        cout << "Max CR: " << max_CR << endl;
        if (measure_perf) {
            cout << "Min decompression performance: " << min_decomp_perf << " MB/s" << endl;
            cout << "Max decompression performance: " << max_decomp_perf << " MB/s" << endl;
            cout << "Min compression performance: " << min_comp_perf << " MB/s" << endl;
            cout << "Max compression performance: " << max_comp_perf << " MB/s" << endl;
        } else {
            cout << "Performance measurement is skipped due to small input size: " << N << " that is less than: " << N_performance_is_mearured << endl;
        }
        iters++;
        cout << "Test is Ok! Total iterations: " << iters << endl;
    }

    return 0;
}