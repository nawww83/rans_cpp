#pragma once

#include <cstdint>
#include <type_traits>
#include <array>
#include <cassert>
#include <algorithm> // std::clamp
#include <cstring>   // std::memcpy
#include <cstdint>
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace cnt {

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

constexpr int M = 256; // byte-wise counter

static int int_log(uint64_t x) {
    if (x <= 1) return 0;
#ifdef _MSC_VER
    unsigned long leading_zero = 0;
    // MSVC ищет индекс самого значимого бита (справа налево)
    if (_BitScanReverse64(&leading_zero, x - 1))
        return static_cast<int>(leading_zero + 1);
    return 0;
#else
    // Clang/GCC считают количество нулей слева
    return 64 - __builtin_clzll(x - 1);
#endif
}

template <u32 Lmax, u32 Lmin=256>
class Counter {
public:
    explicit Counter() = default;
    struct Result {
        std::array<u32, M> frequencies;
        u32 Lscale;
    };
    Result count(const u8* input, size_t n) {
        std::array<u32, M> f{};
        const std::pair<size_t, size_t> q { n / 4, n % 4};
        for (size_t i = 0; i < q.first; ++i) {
            u32 A;
            std::memcpy(&A, input + i * 4, 4);
            f[u8(A)]++;
            f[u8(A >> 8)]++;
            f[u8(A >> 16)]++;
            f[u8(A >> 24)]++;
        }
        for (size_t i = 0; i < q.second; ++i) {
            f[input[q.first * 4 + i]]++;
        }
        // Renormalization.
        std::array<u64, M + 1> cs {0};
        u32 min_f = u32(-1); // minimal non-zero frequency
        for (int i = 0; i < M; ++i) {
            cs[i + 1] = cs[i] + f[i];
            if (f[i] > 0 && f[i] < min_f) min_f = f[i];
        }
        assert( min_f != u32(-1) );
        const u64 N = cs[M];
        u32 L = (u32)(N / min_f) + 1; // Min. L when dynamic range correction is not needed: L > floor( N / minimal f)
        L = 1u << int_log(L);
        L = std::clamp(L, Lmin, Lmax);
        int exceed = 0;
        for (int i = 0; i < M; ++i) {
            const bool is_non_zero = f[i] > 0;
            cs[i + 1] *= L;
            cs[i + 1] /= N;
            f[i] = cs[i + 1] - cs[i];
            const bool is_zero = f[i] == 0;
            if (is_non_zero && is_zero) { // dynamic range failure
                f[i] = 1;
                exceed++;
            }
        }
        // correct the dynamic range failures
        for (int i = 0; i < M; ++i) {
            if (exceed > 0) {
                if (f[i] > 1) {
                    f[i] -= 1;
                    exceed -= 1;
                }
            }
        }
        return {f, L};
    }
};

}