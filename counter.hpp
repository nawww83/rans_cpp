#pragma once

#include <cstdint>
#include <type_traits>
#include <array>


namespace cnt {

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

constexpr int M = 256;

template <typename T>
static bool all_bytes_equal(T x) {
    return ((((x << 8) ^ x) >> 8) == 0);
}


template <u32 L>
class Counter {
public:
    explicit Counter() = default;
    auto get_L() const { return L; }
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
        for (int i=0; i<M; ++i) {
            cs[i+1] = cs.at(i) + f.at(i);
        }
        const u64 N = cs[M];
        int exceed = 0;
        for (int i=0; i<M; ++i) {
            const bool is_non_zero = (f.at(i) > 0);
            cs[i+1] *= L;
            cs[i+1] /= N;
            f[i] = cs.at(i+1) - cs.at(i);
            const bool is_zero = (f.at(i) == 0);
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