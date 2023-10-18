#include <iostream>
#include <random>

#include "counter.hpp"
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
        cout << "dt: " << dt << " ns. Perf: " << (1.e3 * n / dt) << " MB/s" << endl;
    };

    constexpr int N = 1024*1024*32;
    std::vector<u8> v(N);

    constexpr int Q = 16;
    std::vector<double> perfomances(Q);

    cnt::Counter c;

    cout << "Test is started. please, wait..." << endl;
    for (int q=0; q<=Q; ++q) {
        const double par = double(q) / double(Q);
        GeometricDistribution<int> r(par);

        // timer.reset();
        fill_array_randomly(r, v.data(), N);
        // dt = timer.elapsed_ns();
        // disp_perf(N);
        //
        timer.reset();
        auto f { c.count(v.data(), N) };
        dt = timer.elapsed_ns();
        cnt::u64 sum_f = 0;
        for (auto el : f) {
            sum_f += el;
        }
        assert(sum_f == N);
        cout << "par: " << par << ", sum f: " << sum_f << ", ";
        disp_perf(N);
    }

    // for (int i=0; i<cnt::M; ++i) {
    //     cout << c.get_f(i) << ", ";
    // } 
    // cout << endl;

    cout << "Test is Ok!" << endl;

    return 0;
}