#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cassert>

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

template <typename T>
class GeometricDistribution {
private:
    std::geometric_distribution<T> dist;
    std::mt19937 gen;
public:
    GeometricDistribution(double p): gen(rd()), dist(p) {}
    GeometricDistribution(GeometricDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
    T operator()() { return dist(gen); }
};

template<class Generator>
void fill_array_randomly(Generator& g, u8 *v, size_t n, bool inversion = false) {
    // Выносим инверсию наружу или делаем её маской для скорости
    const u8 xor_mask = inversion ? 0xFF : 0x00;

    for (size_t i = 0; i < n; ++i) {
        // g() возвращает T, приводим к u8 и применяем XOR для инверсии
        // (255 - x) для байта это то же самое, что x ^ 0xFF
        v[i] = static_cast<u8>(g() & 0xFF) ^ xor_mask;
    }
}

static timer_n::Timer timer;


int main() {
    io_u::io_utils io;
    std::cout << "Welcome!\n";
    std::cout << "Endianess: " << (io.is_little_endian() ? "LE\n" : "Not LE\n");
    std::cout << "Endianess: " << (io.is_big_endian() ? "BE\n" : "Not BE\n");
    assert( io.is_little_endian() ^ io.is_big_endian());
    constexpr int Q = 128;
    std::vector<GeometricDistribution<int>> gd;
    gd.reserve(Q + 1);
    for (int q = 0; q <= Q; ++q) {
        // Нужен вариант с очень низкой вероятностью вместо 0:
        double par = std::max(0.0001, double(q) / Q);
        gd.emplace_back( par );
    }
    double dt = 0;
    auto calc_perf = [&dt](int n) {
        return (1.e3 * n / dt); // MB/s
    };
    rans::Rans rans;
    long long iters = 0;
    double perf_iters = 0;
    constexpr size_t N_max = 8*1024*1024; // slow test
    // constexpr size_t N_max = 256*1024; // meduim test
    // constexpr size_t N_max = 512; // ultra fast test
    const size_t N_performance_is_mearured = N_max / 4;
    std::vector<u8> v; v.reserve(N_max);
    std::vector<u8> output; output.reserve(rans.required_bytes(N_max));
    std::vector<u8> v_decoded; v_decoded.reserve(N_max);
    double global_min_CR = 1.e18;
    double global_median_compression_perf = 0;
    double global_median_decompression_perf = 0;
    for (;;) {
        size_t N = (((size_t)std::rand()) % N_max) + 1;
        const bool measure_perf = (N >= N_performance_is_mearured);
        std::cout << std::endl;
        std::cout << "New iteration: required output size: " << rans.required_bytes(N) << ", N: " << N << std::endl;
        const bool inv_f = std::rand() % 2;
        std::cout << " Test is started, inversion flag: " << inv_f << ", please, wait..." << std::endl;
        std::vector<double> CRs {};
        std::vector<double> comp_perfs {};
        std::vector<double> decomp_perfs {};
        v.resize(N); v_decoded.resize(N);
        output.resize( rans.required_bytes(N) );
        std::vector<double> redundancies {};
        redundancies.reserve(Q + 1);
        for (int q = 0; q <= Q; ++q) {
            fill_array_randomly(gd[q], v.data(), N, inv_f);

            rans::u32 out_size;
            timer.reset();
            rans.encode(v.data(), N, output.data(), out_size);
            dt = timer.elapsed_ns();

            double cr = double(N) / double(out_size);
            CRs.push_back(cr);

            // Расчет Redundancy
            cnt::Counter<65536> stat_cnt;
            auto [f_norm, L] = stat_cnt.count(v.data(), N);
            double entropy = calculate_entropy(f_norm, L);

            if (entropy > 0.001) {
            double ideal_cr = 8.0 / entropy;
                redundancies.push_back((ideal_cr / cr - 1.0) * 100.0);
            } else {
                redundancies.push_back(0.0); // Почти нулевая энтропия
            }

            assert( out_size <= (int)output.size() );
            if (measure_perf)
                comp_perfs.push_back(calc_perf(N));                        
            timer.reset();
            rans.decode(output.data(), out_size, v_decoded.data(), N);
            dt = timer.elapsed_ns();
            const bool is_equal = std::equal(v.begin(), v.end(), v_decoded.begin());
            if (!is_equal) {
                std::cout << "Input: ";
                for (auto el : v) {
                    std::cout << int(el) << ", ";
                }
                std::cout << std::endl;
                std::cout << "Output size: " << out_size << std::endl;
                std::cout << "Output: ";
                for (auto el : output) {
                    std::cout << int(el) << ", ";
                }
                std::cout << std::endl;
            }
            assert( is_equal );
            if (measure_perf)
                decomp_perfs.push_back(calc_perf(N));            
        }
        // Get some statistics
        std::sort(CRs.begin(), CRs.end());
        std::sort(comp_perfs.begin(), comp_perfs.end());
        std::sort(decomp_perfs.begin(), decomp_perfs.end());
        std::sort(redundancies.begin(), redundancies.end());

        // ... (после сортировки векторов) ...

        const auto NN  = CRs.size();          // Сколько всего тестов (обычно 129)
        const auto NNc = comp_perfs.size();   // В скольких из них замеряли скорость

        // ОБНОВЛЕНИЕ ГЛОБАЛЬНОГО МИНИМУМА
        // CRs[0] — это самый маленький CR в текущей итерации (после сортировки)
        if (CRs[0] < global_min_CR) {
            global_min_CR = CRs[0];
        }

        if (measure_perf && NNc > 0) {
            // Накопительное среднее медиан
            global_median_compression_perf += (comp_perfs[NNc / 2] - global_median_compression_perf) / (perf_iters + 1.0);
            global_median_decompression_perf += (decomp_perfs[NNc / 2] - global_median_decompression_perf) / (perf_iters + 1.0);
            perf_iters += 1.0; // Убедитесь, что perf_iters — это double
        }

        // Увеличиваем общую ширину разделителей (примерно до 70-75 символов)
        printf("\n============================ ITERATION %-5lld ============================\n", iters + 1);
        printf(" N: %-8zu bytes | Global Min CR: %-6.3f | Inv: %s\n", 
                N, global_min_CR, inv_f ? "Yes" : "No ");
        printf("----------------------------------------------------------------------------\n");
        // Заголовки: METRIC (15), остальное по 14
        printf(" %-15s | %-14s | %-14s | %-14s\n", "METRIC", "MIN", "MAX", "MEDIAN");
        printf(" ----------------|----------------|----------------|----------------\n");

        // Для чисел используем %14.3f (всего 14 символов, 3 после запятой)
        printf(" Comp. Ratio     | %14.3f | %14.3f | %14.3f\n", CRs[0], CRs[NN-1], CRs[NN/2]);
        printf(" Redundancy (%%)  | %14.2f | %14.2f | %14.2f\n", redundancies[0], redundancies[NN-1], redundancies[NN/2]);

        if (measure_perf && !comp_perfs.empty()) {
            // Для скорости используем %14.0f, так как дробная часть в MB/s обычно не критична
            printf(" Enc (MB/s)      | %14.0f | %14.0f | %14.0f\n", comp_perfs[0], comp_perfs[NNc-1], comp_perfs[NNc/2]);
            printf(" Dec (MB/s)      | %14.0f | %14.0f | %14.0f\n", decomp_perfs[0], decomp_perfs[NNc-1], decomp_perfs[NNc/2]);
        }

        iters++;
        printf("----------------------------------------------------------------------------\n");

        // Вывод финальной строки итерации
        printf("PASSED. Iterations: %-6lld | Global CR > %-6.3f\n", iters, global_min_CR);
        if (perf_iters > 0) {        
            printf("Global Median: Enc %10.0f MB/s | Dec %10.0f MB/s (over %5.0f measurements)\n", 
                    global_median_compression_perf, global_median_decompression_perf, perf_iters);
        }
        printf("============================================================\n");

    }

    return 0;
}