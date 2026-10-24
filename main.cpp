#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <array>
#include <algorithm>

#include "counter.hpp"
#include "rans.hpp"
#include "timer.hpp"

using u8 = cnt::u8;

double calculate_entropy(const std::array<cnt::u32, cnt::M>& f, size_t L) {
    double entropy = 0;
    for (auto count : f) {
        if (count > 0) {
            double p = static_cast<double>(count) / L;
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

static std::random_device rd{};
// Отдельный ГПСЧ для размера блока и флага инверсии: std::rand() имеет
// маленький период и общий глобальный state.
static std::mt19937_64 rng64{rd()};

template <typename T>
class GeometricDistribution {
private:
    std::geometric_distribution<T> dist;
    std::mt19937 gen;
public:
    GeometricDistribution(double p): dist(p), gen(rd()) {}
    GeometricDistribution(GeometricDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
    T operator()() { return dist(gen); }
};

template<class Generator>
void fill_array_randomly(Generator& g, u8 *v, size_t n, bool inversion = false) {
    const u8 xor_mask = inversion ? 0xFF : 0x00;
    // Младшие 8 бит k ~ Geom(p) дают в точности геометрическое распределение,
    // усечённое до [0, 255] и перенормированное:
    //   P(b) = p (1-p)^b / (1 - (1-p)^256).
    // Свёрнутый хвост k >= 256 — геометрическая прогрессия, поэтому это не
    // потеря точности, а точная выборка из усечённого распределения.
    for (size_t i = 0; i < n; ++i) {
        v[i] = static_cast<u8>(g() & 0xFF) ^ xor_mask;
    }
}

static timer_n::Timer timer;

constexpr std::array<qnt::Mode, 4> kModes = {
    qnt::Mode::Precise, qnt::Mode::Pow2, qnt::Mode::Pow23, qnt::Mode::Pow235
};

struct ModeStats {
    std::vector<double> CRs;
    std::vector<double> redundancies;
    std::vector<double> comp_perfs;
    std::vector<double> decomp_perfs;
    double global_min_CR = 1e18;
    double global_median_comp = 0;
    double global_median_decomp = 0;
};

int main() {
    io_u::io_utils io;
    std::cout << "Welcome!\n";
    std::cout << "Endianess: " << (io.is_little_endian() ? "LE\n" : "Not LE\n");
    std::cout << "Endianess: " << (io.is_big_endian() ? "BE\n" : "Not BE\n");
    assert( io.is_little_endian() ^ io.is_big_endian() );
    constexpr int Q = 128;
    std::vector<GeometricDistribution<int>> gd;
    gd.reserve(Q + 1);
    for (int q = 0; q <= Q; ++q) {
        double par = std::max(0.0001, double(q) / Q);
        gd.emplace_back( par );
    }

    // Отдельные прогоны для фиксированных p: от почти равномерного
    // (p -> 0, геометрическое распределение заворачивается по модулю 256)
    // до сильно скошенного (p -> 1, почти все байты равны 0).
    constexpr std::array<double, 9> kFixedP = {
        0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.6, 0.9
    };
    std::vector<GeometricDistribution<int>> fixed_gd;
    fixed_gd.reserve(kFixedP.size());
    for (double p : kFixedP) fixed_gd.emplace_back(p);
    double dt = 0;
    auto calc_perf = [&dt](size_t n) {
        return (1.e3 * static_cast<double>(n) / dt); // MB/s
    };
    rans::Rans rans;
    long long iters = 0;
    constexpr size_t N_max = 8*1024*1024; // slow test
    // constexpr size_t N_max = 256*1024; // medium test
    // constexpr size_t N_max = 512; // ultra fast test
    const size_t N_performance_is_measured = N_max / 4;
    std::vector<u8> v; v.reserve(N_max);
    std::vector<u8> output; output.reserve(rans.required_bytes(N_max));
    std::vector<u8> v_decoded; v_decoded.reserve(N_max);

    std::array<ModeStats, kModes.size()> global;
    double perf_iters = 0;

    for (;;) {
        size_t N = (((size_t)rng64()) % N_max) + 1;
        const bool measure_perf = (N >= N_performance_is_measured);
        std::cout << std::endl;
        std::cout << "New iteration: required output size: " << rans.required_bytes(static_cast<rans::u32>(N)) << ", N: " << N << std::endl;
        const bool inv_f = (rng64() & 1) != 0;
        std::cout << " Test is started, inversion flag: " << inv_f << ", please, wait..." << std::endl;

        std::array<ModeStats, kModes.size()> stats;
        v.resize(N); v_decoded.resize(N);
        output.resize(rans.required_bytes(static_cast<rans::u32>(N)));

        for (int q = 0; q <= Q; ++q) {
            fill_array_randomly(gd[q], v.data(), N, inv_f);

            // Истинная энтропия считается один раз для всех режимов.
            cnt::Counter<65536> stat_cnt;
            auto [f_norm, L] = stat_cnt.count(v.data(), N);
            const double entropy = calculate_entropy(f_norm, L);
            const bool entropy_ok = entropy > 0.001;
            const double ideal_cr = entropy_ok ? 8.0 / entropy : 0.0;

            for (size_t mi = 0; mi < kModes.size(); ++mi) {
                const qnt::Mode mode = kModes[mi];
                rans::u32 out_size = 0;
                timer.reset();
                rans.encode(v.data(), (rans::u32)N, output.data(), out_size, mode);
                dt = timer.elapsed_ns();

                double cr = double(N) / double(out_size);
                stats[mi].CRs.push_back(cr);
                if (entropy_ok) {
                    stats[mi].redundancies.push_back((ideal_cr / cr - 1.0) * 100.0);
                } else {
                    stats[mi].redundancies.push_back(0.0);
                }

                assert( out_size <= (int)output.size() );
                if (measure_perf) stats[mi].comp_perfs.push_back(calc_perf(N));

                timer.reset();
                rans.decode(output.data(), out_size, v_decoded.data(), static_cast<rans::u32>(N));
                dt = timer.elapsed_ns();
                const bool is_equal = std::equal(v.begin(), v.end(), v_decoded.begin());
                if (!is_equal) {
                    std::cout << "Mode " << qnt::mode_name(mode) << " FAILED round-trip!\n";
                }
                assert( is_equal );
                if (measure_perf) stats[mi].decomp_perfs.push_back(calc_perf(N));
            }
        }

        for (size_t mi = 0; mi < kModes.size(); ++mi) {
            auto& s = stats[mi];
            auto& gm = global[mi];
            std::sort(s.CRs.begin(), s.CRs.end());
            std::sort(s.redundancies.begin(), s.redundancies.end());
            std::sort(s.comp_perfs.begin(), s.comp_perfs.end());
            std::sort(s.decomp_perfs.begin(), s.decomp_perfs.end());

            if (s.CRs.front() < gm.global_min_CR) gm.global_min_CR = s.CRs.front();
            if (measure_perf && !s.comp_perfs.empty()) {
                gm.global_median_comp   += (s.comp_perfs[s.comp_perfs.size()/2] - gm.global_median_comp) / (perf_iters + 1.0);
                gm.global_median_decomp += (s.decomp_perfs[s.decomp_perfs.size()/2] - gm.global_median_decomp) / (perf_iters + 1.0);
            }
        }
        if (measure_perf) perf_iters += 1.0;

        printf("\n============================ ITERATION %-5lld ============================\n", iters + 1);
        printf(" N: %-8zu bytes | Inv: %s | Tests per mode: %zu\n", N, inv_f ? "Yes" : "No ", kModes.size());
        printf("--------------------------------------------------------------------------------------------\n");
        printf(" %-11s | %8s | %8s | %9s | %10s | %9s\n", "MODE", "CR min", "CR med", "Red% min", "Red% med", "Enc MB/s");
        printf("-------------|----------|----------|-----------|------------|----------\n");
        for (size_t mi = 0; mi < kModes.size(); ++mi) {
            auto& s = stats[mi];
            const double enc = (measure_perf && !s.comp_perfs.empty()) ? s.comp_perfs[s.comp_perfs.size()/2] : 0.0;
            printf(" %-11s | %8.3f | %8.3f | %9.2f | %10.2f | %9.0f\n",
                   qnt::mode_name(kModes[mi]),
                   s.CRs.front(), s.CRs[s.CRs.size()/2],
                   s.redundancies.front(), s.redundancies[s.redundancies.size()/2],
                   enc);
        }
        printf("--------------------------------------------------------------------------------------------\n");
        printf(" %-11s", "Dec MB/s:");
        for (size_t mi = 0; mi < kModes.size(); ++mi) {
            auto& s = stats[mi];
            const double dec = (measure_perf && !s.decomp_perfs.empty()) ? s.decomp_perfs[s.decomp_perfs.size()/2] : 0.0;
            printf(" %s=%.0f", qnt::mode_name(kModes[mi]), dec);
        }
        printf("\n");

        // Отдельные прогоны для фиксированных p: один и тот же N и inversion
        // flag для всех режимов, поэтому строки можно сравнивать напрямую.
        if (measure_perf) {
            printf("\n Fixed-p benchmarks (N=%zu, inv=%d):\n", N, inv_f ? 1 : 0);
            printf(" %-8s %6s | %-11s | %8s %9s %9s %9s\n",
                   "p", "H", "MODE", "CR", "Red%", "Enc", "Dec");
            printf("----------|--------|-------------|----------|----------|----------|----------\n");
            for (size_t pi = 0; pi < kFixedP.size(); ++pi) {
                fill_array_randomly(fixed_gd[pi], v.data(), N, inv_f);

                cnt::Counter<65536> stat_cnt;
                auto [f_norm, L] = stat_cnt.count(v.data(), N);
                const double H = calculate_entropy(f_norm, L);
                const bool h_ok = H > 0.001;
                const double ideal_cr = h_ok ? 8.0 / H : 0.0;

                for (size_t mi = 0; mi < kModes.size(); ++mi) {
                    rans::u32 out_size = 0;
                    timer.reset();
                    rans.encode(v.data(), (rans::u32)N, output.data(), out_size, kModes[mi]);
                    dt = timer.elapsed_ns();
                    const double enc = calc_perf(N);

                    timer.reset();
                    rans.decode(output.data(), out_size, v_decoded.data(), (rans::u32)N);
                    dt = timer.elapsed_ns();
                    const double dec = calc_perf(N);

                    const bool is_equal = std::equal(v.begin(), v.end(), v_decoded.begin());
                    if (!is_equal) {
                        std::cout << "Mode " << qnt::mode_name(kModes[mi]) << " FAILED round-trip!\n";
                    }
                    assert( is_equal );

                    const double cr = double(N) / double(out_size);
                    const double red = h_ok ? (ideal_cr / cr - 1.0) * 100.0 : 0.0;
                    printf(" %-8.4f %6.2f | %-11s | %8.3f %9.2f %9.0f %9.0f\n",
                           kFixedP[pi], H, qnt::mode_name(kModes[mi]), cr, red, enc, dec);
                }
            }
        }

        iters++;
        printf(" Global min CR:");
        for (size_t mi = 0; mi < kModes.size(); ++mi) {
            printf(" %s=%.3f", qnt::mode_name(kModes[mi]), global[mi].global_min_CR);
        }
        printf("\n");
        printf(" Global median:");
        for (size_t mi = 0; mi < kModes.size(); ++mi) {
            printf(" %s Enc=%.0f/Dec=%.0f", qnt::mode_name(kModes[mi]),
                   global[mi].global_median_comp, global[mi].global_median_decomp);
        }
        printf("  (over %.0f measurements)\n", perf_iters);
        printf("============================================================\n");
    }

    return 0;
}
