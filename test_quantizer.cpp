#include <array>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <random>
#include <vector>

#include "counter.hpp"
#include "quantizer.hpp"
#include "rans.hpp"

using namespace qnt;

static std::mt19937 rng(12345);

static void test_reciprocal() {
    // Проверяем ровно ту формулу, что в энкодере: floor(y / odd) ==
    // (y * ceil(2^48 / odd)) >> 48 для всех y < odd * Bflush.
    // После ренормализации x < f * Bflush, поэтому y = x >> tz < odd * Bflush.
    // Перебираем все нечётные odd < 2^16 — это покрывает любую частоту f.
    for (u32 odd = 1; odd < 65536; odd += 2) {
        const u64 recip = ((u64(1) << 48) + odd - 1) / odd;
        const u64 limit = (u64)odd * 65536; // y < limit
        auto check = [&](u32 y) {
            const u32 q = (u32)(((u64)y * recip) >> 48);
            if (q != y / odd) {
                std::printf("RECIP FAIL odd=%u y=%u got=%u want=%u\n", odd, y, q, y / odd);
                assert(false);
            }
        };
        // Границы диапазона и окрестности кратних odd.
        check(0);
        check((u32)(limit - 1));
        for (u32 k = 0; k < 8; ++k) {
            const u32 y = (u32)(limit - 1 - k);
            check(y);
            if (y > 0) check(y - 1);
        }
        for (u32 k = 0; k < 8; ++k) {
            const long long base = (long long)odd * k;
            for (long long d = -2; d <= 2; ++d) {
                const long long v = base + d;
                if (v >= 0 && (u64)v < limit) check((u32)v);
            }
        }
        for (int t = 0; t < 8; ++t) check((u32)(rng() % limit));
    }
    std::printf("reciprocal: OK\n");
}

static void test_quantizer() {
    for (int mode_i = 0; mode_i < 4; ++mode_i) {
        const Mode mode = (Mode)mode_i;
        for (int t = 0; t < 20; ++t) {
            std::array<u32, M> raw{};
            const int n = 1 + (int)(rng() % M);
            u32 N = 0;
            // Геометрически убывающие счётчики по случайному подмножеству.
            for (int i = 0; i < n; ++i) {
                raw[i] = 1 + (rng() % 1000000);
                N += raw[i];
            }
            // Сделаем распределение более скошенным.
            for (int i = 0; i < n; ++i) raw[i] = std::max<u32>(1, raw[i] / (u32)(i + 1));
            N = 0;
            for (int i = 0; i < M; ++i) N += raw[i];

            const u32 L = cnt::Counter<rans::Lmax>::choose_L(raw);
            std::array<u32, M> f{};
            quantize(raw, L, mode, f);

            u64 sum = 0;
            for (int i = 0; i < M; ++i) sum += f[i];
            assert(sum == L);

            for (int i = 0; i < M; ++i) {
                if (raw[i] > 0) assert(f[i] > 0);
                else { assert(f[i] == 0); continue; }
                if (mode == Mode::Pow2) {
                    assert((f[i] & (f[i] - 1)) == 0); // степень двойки
                } else if (mode == Mode::Pow23) {
                    u32 v = f[i];
                    while (v % 2 == 0) v /= 2;
                    while (v % 3 == 0) v /= 3;
                    assert(v == 1);
                } else if (mode == Mode::Pow235) {
                    u32 v = f[i];
                    while (v % 2 == 0) v /= 2;
                    while (v % 3 == 0) v /= 3;
                    while (v % 5 == 0) v /= 5;
                    assert(v == 1);
                }
            }
        }
    }
    std::printf("quantizer: OK\n");
}

static void test_roundtrip() {
    rans::Rans rans;
    for (int mode_i = 0; mode_i < 4; ++mode_i) {
        const Mode mode = (Mode)mode_i;
        for (int t = 0; t < 60; ++t) {
            const size_t N = 1 + (rng() % 200000);
            std::vector<rans::u8> v(N), out(rans::Rans().required_bytes((rans::u32)N));
            std::vector<rans::u8> dec(N);
            // Разные распределения, включая крайне скошенные.
            const unsigned kind = rng() % 4;
            for (size_t i = 0; i < N; ++i) {
                if (kind == 0) v[i] = (rans::u8)(rng() & 0xFF);
                else if (kind == 1) v[i] = (rans::u8)(rng() % 8);
                else if (kind == 2) v[i] = (rng() % 100 == 0) ? (rans::u8)(rng() & 0xFF) : 0;
                else v[i] = (rans::u8)((i * 7) & 0xFF);
            }
            rans::u32 out_size = 0;
            rans.encode(v.data(), (rans::u32)N, out.data(), out_size, mode);
            assert(out_size <= out.size());
            rans.decode(out.data(), out_size, dec.data(), (rans::u32)N);
            if (v != dec) {
                std::printf("ROUNDTRIP FAIL mode=%s N=%zu kind=%u\n", mode_name(mode), N, kind);
                assert(false);
            }
        }
        std::printf("roundtrip %-10s: OK\n", mode_name(mode));
    }
}

static void test_edge_cases() {
    rans::Rans rans;
    std::vector<std::vector<rans::u8>> cases;
    cases.push_back({0});
    cases.push_back(std::vector<rans::u8>(1000, 0));
    cases.push_back({0, 1});
    cases.push_back({7, 7, 7, 8});
    std::vector<rans::u8> three(1000);
    for (size_t i = 0; i < three.size(); ++i) three[i] = (rans::u8)(i % 3);
    cases.push_back(three);
    for (int m = 0; m < 4; ++m) {
        for (auto& v : cases) {
            std::vector<rans::u8> out(rans.required_bytes((rans::u32)v.size()));
            std::vector<rans::u8> dec(v.size());
            rans::u32 os = 0;
            rans.encode(v.data(), (rans::u32)v.size(), out.data(), os, (qnt::Mode)m);
            assert(os <= out.size());
            rans.decode(out.data(), os, dec.data(), (rans::u32)v.size());
            if (v != dec) {
                std::printf("EDGE FAIL mode=%s N=%zu\n", mode_name((qnt::Mode)m), v.size());
                assert(false);
            }
        }
    }
    std::printf("edge cases: OK\n");
}

int main() {
    setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("start\n");
    test_reciprocal();
    std::printf("phase 2\n");
    test_quantizer();
    std::printf("phase 3\n");
    test_roundtrip();
    std::printf("phase 4\n");
    test_edge_cases();
    std::printf("ALL TESTS PASSED\n");
    return 0;
}
