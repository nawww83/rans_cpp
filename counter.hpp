#pragma once

#include <cstdint>
#include <type_traits>
#include <array>
#include <cassert>
#include <algorithm>


namespace cnt {

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

constexpr int M = 256; // byte-wise counter

template <typename T>
static bool all_bytes_equal(T x) {
    return ((((x << 8) ^ x) >> 8) == 0);
}

template <typename T>
static int int_log(T x) {
    int log_x = 0;
    while (true) {
        if ((1ull << log_x) >= (u64)x) {
            break;
        }
        log_x += 1;
    }
    return log_x;
}


template <u32 Lmax, u32 Lmin=256>
class Counter {
public:
    explicit Counter() = default;
    auto count(const u8* input, int n) {
        std::array<u32, M> f {};
        const std::pair<int, int> q { n / 4, n % 4};
        std::array<u32, 4> idx;
        for (int i=0; i<q.first; ++i) {
            const auto A = *(u32*)(input + 4*i);
            const bool eq = all_bytes_equal(A);
            idx = {u8(A >> 0), u8(A >> 8), u8(A >> 16), u8(A >> 24)};
            f[idx[0]] += eq ? 4 : 1;
            f[idx[1]] += eq ? 0 : 1;
            f[idx[2]] += eq ? 0 : 1;
            f[idx[3]] += eq ? 0 : 1;
        }
        for (int i=0; i<q.second; ++i) {
            f[input[q.first*4 + i]] += 1;
        }
        // renormalize
        std::array<u64, M+1> cs {0};
        // std::cout << "Unnormalized frequencies: ";
        u32 min_f = u32(-1); // minimal non-zero frequency
        for (int i=0; i<M; ++i) {
            cs[i+1] = cs[i] + f[i];
            if (f[i] == 0) {
                continue;
            }
            min_f = (f[i] < min_f) ? f[i] : min_f;
            // std::cout << f[i] << ", ";
        }
        // std::cout << std::endl;
        assert( min_f != u32(-1) );
        const u64 N = cs[M];
        u32 L = (N / min_f) + 1; // minimal L when dynamic range correction is not needed: L > floor( N / minimal f)
        // std::cout << "Counter L: " << L << ", min f: " << min_f << ", N: " << N << std::endl;
        const int logL = int_log(L);
        L = (1u << logL);
        // std::cout << " CC L: " << L << std::endl;
        L = std::min(L , Lmax);
        L = std::max(L, Lmin);
        // std::cout << " CCC L: " << L << std::endl;
        int exceed = 0;
        for (int i=0; i<M; ++i) {
            const bool is_non_zero = (f[i] > 0);
            cs[i + 1] *= L;
            cs[i + 1] /= N;
            f[i] = cs[i + 1] - cs[i];
            const bool is_zero = (f[i] == 0);
            if (is_non_zero && is_zero) { // dynamic range failure
                f[i] = 1;
                exceed++;
            }
        }
        // correct the dynamic range failures
        for (int i=0; i<M; ++i) {
            if (exceed > 0) {
                if (f[i] > 1) {
                    f[i] -= 1;
                    exceed -= 1;
                }
            }
        }
        //
        return f;
    }
};

}