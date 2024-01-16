#include <iostream>
#include <random>
#include <vector>
#include <cassert>

#include "counter.hpp"
#include "rans.hpp"
#include "timer.hpp"

using u8 = cnt::u8;

static std::random_device rd;

template <typename T>
class GeometricDistribution {
private:
    std::mt19937 gen;
    std::geometric_distribution<T> dist;
public:
    GeometricDistribution(double p): gen(rd()), dist(p) {}
    GeometricDistribution(GeometricDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
    T operator()() { return dist(gen); }
};


template<class Generator>
auto fill_array_randomly(Generator& g, u8 *v, size_t n, bool inversion=false) {
    const auto div = n / 2;
    const auto rem = n % 2;
    for (size_t i=0; i<div; ++i) {
        v[2*i] = (! inversion) ? (g() & 255) : (255 - (g() & 255));
        v[2*i+1] = (! inversion) ? (g() & 255) : (255 - (g() & 255));
    }
    for (size_t i=0; i<rem; ++i) {
        v[div*2 + i] = (! inversion) ? (g() & 255) : (255 - (g() & 255));
    }
}


static timer_n::Timer timer;

int main() {
    using namespace std;
    io_u::io_utils io;

    const auto seed = std::time(0);
    std::srand(seed);
    cout << "Welcome!\n";
    cout << "Seed: " << seed << "\n";
    cout << "Endianess: " << (io.is_little_endian() ? "LE\n" : "Not LE\n");
    cout << "Endianess: " << (io.is_big_endian() ? "BE\n" : "Not BE\n");
    assert( io.is_little_endian() ^ io.is_big_endian());
    //
    constexpr int Q = 128;
    std::vector<GeometricDistribution<int>> gd{};
    for (int q=0; q<=Q; ++q) {
        double par = double(q) / double(Q);
        gd.emplace_back( GeometricDistribution<int>(par) );
    }
    //

    double dt = 0;
    auto calc_perf = [&dt](int n) {
        return (1.e3 * n / dt); // MB/s
    };
    rans::Rans rans;
    long long iters = 0;
    double perf_iters = 0;
    constexpr size_t N_max = 8*1024*1024; // slow test
    // constexpr size_t N_max = 4*1024; // fast test
    const size_t N_performance_is_mearured = N_max / 4;
    std::vector<u8> v; v.reserve(N_max);
    std::vector<u8> output; output.reserve(rans.required_bytes(N_max));
    std::vector<u8> v_decoded; v_decoded.reserve(N_max);
    double global_min_CR = 1.e18;
    double global_median_compression_perf = 0;
    double global_median_decompression_perf = 0;
    while (true) {
        size_t N = (((size_t)std::rand()) % N_max) + 1;
        const bool measure_perf = (N >= N_performance_is_mearured);
        std::cout << std::endl;
        std::cout << "New iteration: required output size: " << rans.required_bytes(N) << ", N: " << N << std::endl;
        const bool inv_f = std::rand() % 2;
        cout << " Test is started, inversion flag: " << inv_f << ", please, wait..." << endl;
        std::vector<double> CRs {};
        std::vector<double> comp_perfs {};
        std::vector<double> decomp_perfs {};
        v.resize(N); v_decoded.resize(N);
        output.resize( rans.required_bytes(N) );
        for (int q=0; q<=Q; ++q) {
            // const double par = double(q) / double(Q);
            // GeometricDistribution<int> r(par);
            fill_array_randomly(gd[q], v.data(), N, inv_f);
            rans::u32 out_size;
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
        const auto NN = CRs.size();
        const auto NNc = comp_perfs.size();
        assert(comp_perfs.size() == decomp_perfs.size());
        global_min_CR = std::min(global_min_CR, CRs[0]);
        cout << endl;
        cout << "Input length: " << N << " bytes" << endl;
        cout << "Compression Ratio, CR:" << endl;
        cout << " Min: " << CRs[0] << ", Max: " << CRs[NN-1] << ", Median: " << CRs[NN / 2] << endl;
        if (measure_perf) {
            cout << "Decompression performance, MB/s:" << endl;
            cout << " Min: " << decomp_perfs[0] << ", Max: " << decomp_perfs[NNc-1] << ", Median: " << decomp_perfs[NNc / 2] << endl;
            cout << "Compression performance, MB/s:" << endl;
            cout << " Min: " << comp_perfs[0] << ", Max: " << comp_perfs[NNc-1] << ", Median: " << comp_perfs[NNc / 2] << endl;

            global_median_compression_perf += (comp_perfs[NNc / 2] - global_median_compression_perf) / (perf_iters + 1);
            global_median_decompression_perf += (decomp_perfs[NNc / 2] - global_median_decompression_perf) / (perf_iters + 1);
            perf_iters += 1;
        } else {
            cout << "Performance measurement is skipped due to small input size: " << N << " that is less than: " << N_performance_is_mearured << endl;
        }
        iters++;
        cout << "Test is Ok! Total iterations: " << iters << ", global min CR: " << global_min_CR << endl;
        if (perf_iters > 0) {
            cout << "Average median performance: compression: " << global_median_compression_perf << 
                    ", decompression: " << global_median_decompression_perf << 
                    " MB/s, perf. iters: " << perf_iters << endl;
        }
    }

    return 0;
}