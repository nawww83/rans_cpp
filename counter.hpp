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
        return f;
    }
};

}