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
auto fill_array_randomly(Generator& g, u8 *v, int n) {
    for (int i=0; i<n; ++i) {
        v[i] = g() & 255;
    }
}


static timer_n::Timer timer;

int main() {
    using namespace std;

    double dt = 0;
    auto disp_perf = [&dt](int n) {
        cout << " performance: " << (1.e3 * n / dt) << " MB/s" << endl;
    };

    rans::Rans rans;
    // cnt::Counter<65536> c;
    constexpr int N_max = 8*1024*1024;

    while (true) {

        size_t N = (((size_t)std::rand()) % N_max) + 1;
        std::vector<u8> v(N);

        std::vector<u8> output(rans.required_bytes(N));

        std::cout << std::endl;
        std::cout << "New iteration: required size: " << output.size() << ", N: " << N << std::endl;

        std::vector<u8> v_decoded(N);
        constexpr int Q = 32;
        std::vector<double> perfomances(Q);

        cout << " Test is started. please, wait..." << endl;
        for (int q=0; q<=Q; ++q) {
            const double par = double(q) / double(Q);
            GeometricDistribution<int> r(par);

            fill_array_randomly(r, v.data(), N);
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
            cout << " geom. distr. p: " << par << ", compressed size: " << out_size << ", bytes. CR: " << double(N) / (out_size + sizeof(cnt::u32) + cnt::M*sizeof(cnt::u16)) << ", ";
            disp_perf(N);
            //
            timer.reset();
            rans.decode(output.data(), out_size, v_decoded.data(), N);
            dt = timer.elapsed_ns();
            assert( std::equal(v.begin(), v.end(), v_decoded.begin()) );
            cout << "  decompressed size: " << v_decoded.size() << ", bytes" << ", ";
            disp_perf(N);
        }

        cout << "Test is Ok!" << endl;
    }

    return 0;
}