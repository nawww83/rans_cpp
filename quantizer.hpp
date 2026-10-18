#pragma once

#include "counter.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <set>
#include <vector>

namespace qnt {

using u8  = cnt::u8;
using u16 = cnt::u16;
using u32 = cnt::u32;
using u64 = cnt::u64;

constexpr int M = cnt::M;

// Режимы нормализации/квантования частот.
//   Precise  - точная пропорциональная нормализация (аналог арифметического кодирования)
//   Pow2     - частоты только степени двойки (эквивалент кодов Хаффмана, сдвиг вместо деления)
//   Pow23    - частоты вида 2^a * 3^b
//   Pow235   - частоты вида 2^a * 3^b * 5^c
enum class Mode : u8 {
    Precise = 0,
    Pow2    = 1,
    Pow23   = 2,
    Pow235  = 3
};

inline const char* mode_name(Mode m) {
    switch (m) {
        case Mode::Precise: return "Precise";
        case Mode::Pow2:    return "Pow2";
        case Mode::Pow23:   return "Pow{2,3}";
        case Mode::Pow235:  return "Pow{2,3,5}";
    }
    return "?";
}

// Множество "гладких" чисел (все простые делители из заданного набора),
// не превосходящих limit. Хранится в возрастающем порядке.
class SmoothSet {
public:
    std::vector<u32> v;
    // Прямая таблица value -> индекс в v (или -1). Ускоряет поиск в горячем цикле.
    std::vector<int16_t> lut;

    void build(std::initializer_list<u32> primes, u32 limit) {
        std::set<u32> s{1};
        std::vector<u32> frontier{1};
        while (!frontier.empty()) {
            const u32 x = frontier.back();
            frontier.pop_back();
            for (u32 p : primes) {
                const u64 w = static_cast<u64>(x) * p;
                if (w <= limit && s.insert(static_cast<u32>(w)).second) {
                    frontier.push_back(static_cast<u32>(w));
                }
            }
        }
        v.assign(s.begin(), s.end());
        lut.assign(limit + 1, int16_t(-1));
        for (size_t i = 0; i < v.size(); ++i) lut[v[i]] = static_cast<int16_t>(i);
    }

    // Индекс наибольшего значения <= x (x >= 1). -1, если такого нет.
    int idx_le(u32 x) const {
        auto it = std::upper_bound(v.begin(), v.end(), x);
        return static_cast<int>(it - v.begin()) - 1;
    }

    // Индекс значения, равного x, либо -1.
    int index_of(u32 x) const {
        return (x < lut.size()) ? static_cast<int>(lut[x]) : -1;
    }

    bool contains(u32 x) const { return index_of(x) >= 0; }

    // Ближайшее (в логарифмическом масштабе) значение к вещественному r.
    u32 nearest(double r) const {
        if (r <= static_cast<double>(v.front())) return v.front();
        if (r >= static_cast<double>(v.back())) return v.back();
        const u32 x = static_cast<u32>(r);
        int j = idx_le(x);
        if (j < 0) j = 0;
        u32 best = v[j];
        double best_d = std::abs(std::log(static_cast<double>(best) / r));
        for (int k = j - 1; k <= j + 2; ++k) {
            if (k < 0 || k >= static_cast<int>(v.size())) continue;
            const double d = std::abs(std::log(static_cast<double>(v[k]) / r));
            if (d < best_d) { best_d = d; best = v[k]; }
        }
        return best;
    }
};

inline const SmoothSet& smooth23() {
    static const SmoothSet s = [] {
        SmoothSet t;
        t.build({2, 3}, 65536);
        return t;
    }();
    return s;
}

inline const SmoothSet& smooth235() {
    static const SmoothSet s = [] {
        SmoothSet t;
        t.build({2, 3, 5}, 65536);
        return t;
    }();
    return s;
}

// ---------------------------------------------------------------------------
// Length-limited Huffman (package-merge). Возвращает длины кодов.
// Гарантирует: depth[i] >= 1, depth[i] <= logL, sum 2^(logL-depth[i]) == 2^logL.
// ---------------------------------------------------------------------------
struct PMNode {
    u64 w;
    int sym;   // >= 0 для листа (индекс символа), -1 для пакета
    int l, r;  // дети пакета
};

inline void package_merge_depths(const std::vector<u32>& counts,
                                 int logL,
                                 std::vector<int>& depth) {
    const int n = static_cast<int>(counts.size());
    depth.assign(n, 0);
    if (n == 0) return;
    if (n == 1) { depth[0] = 1; return; }

    // Символы, отсортированные по возрастанию веса.
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return counts[a] < counts[b];
    });

    std::vector<PMNode> nodes;
    nodes.reserve(static_cast<size_t>(n) * logL * 2);

    std::vector<int> leaves(n);
    for (int i = 0; i < n; ++i) {
        nodes.push_back({static_cast<u64>(counts[order[i]]), order[i], -1, -1});
        leaves[i] = static_cast<int>(nodes.size()) - 1;
    }

    std::vector<int> lst = leaves;
    std::vector<int> packs, merged;
    for (int level = 2; level <= logL; ++level) {
        packs.clear();
        for (int k = 0; k + 1 < static_cast<int>(lst.size()); k += 2) {
            const u64 w = nodes[lst[k]].w + nodes[lst[k + 1]].w;
            nodes.push_back({w, -1, lst[k], lst[k + 1]});
            packs.push_back(static_cast<int>(nodes.size()) - 1);
        }
        merged.clear();
        merged.reserve(leaves.size() + packs.size());
        size_t a = 0, b = 0;
        while (a < leaves.size() && b < packs.size()) {
            if (nodes[leaves[a]].w <= nodes[packs[b]].w) merged.push_back(leaves[a++]);
            else                                          merged.push_back(packs[b++]);
        }
        while (a < leaves.size()) merged.push_back(leaves[a++]);
        while (b < packs.size())  merged.push_back(packs[b++]);
        lst.swap(merged);
    }

    const int need = 2 * n - 2;
    assert(static_cast<int>(lst.size()) >= need && "L too small for alphabet");

    std::vector<int> stack;
    stack.reserve(64);
    for (int i = 0; i < need; ++i) stack.push_back(lst[i]);
    while (!stack.empty()) {
        const int id = stack.back();
        stack.pop_back();
        if (nodes[id].sym >= 0) {
            depth[nodes[id].sym]++;
        } else {
            stack.push_back(nodes[id].l);
            stack.push_back(nodes[id].r);
        }
    }
    for (int i = 0; i < n; ++i) {
        assert(depth[i] >= 1 && depth[i] <= logL);
    }
}

// ---------------------------------------------------------------------------
// Точная пропорциональная нормализация (режим Precise).
// Повторяет поведение cnt::Counter<Lmax>::count().
// ---------------------------------------------------------------------------
inline void precise_freqs(const std::array<u32, M>& raw, u32 L,
                          std::array<u32, M>& f) {
    std::array<u64, M + 1> cs{0};
    for (int i = 0; i < M; ++i) cs[i + 1] = cs[i] + raw[i];
    const u64 N = cs[M];
    f = raw;
    int exceed = 0;
    for (int i = 0; i < M; ++i) {
        const bool is_non_zero = f[i] > 0;
        cs[i + 1] *= L;
        cs[i + 1] /= N;
        f[i] = static_cast<u32>(cs[i + 1] - cs[i]);
        if (is_non_zero && f[i] == 0) {
            f[i] = 1;
            ++exceed;
        }
    }
    for (int i = 0; i < M && exceed > 0; ++i) {
        if (f[i] > 1) { --f[i]; --exceed; }
    }
    assert(exceed == 0);
}

// ---------------------------------------------------------------------------
// Квантование частотами-степенями двойки (Huffman-эквивалент).
// ---------------------------------------------------------------------------
inline void pow2_freqs(const std::array<u32, M>& raw, u32 L, int logL,
                       std::array<u32, M>& f) {
    f.fill(0);
    std::vector<u32> counts;
    std::vector<int> syms;
    counts.reserve(M);
    syms.reserve(M);
    for (int i = 0; i < M; ++i) {
        if (raw[i] > 0) {
            counts.push_back(raw[i]);
            syms.push_back(i);
        }
    }
    if (counts.size() == 1) {
        // Единственный символ должен занимать весь диапазон L.
        f[syms[0]] = L;
        return;
    }
    std::vector<int> depth;
    package_merge_depths(counts, logL, depth);
    for (size_t i = 0; i < syms.size(); ++i) {
        f[syms[i]] = 1u << (logL - depth[i]);
    }
    // Проверка: суммарная частота обязана быть ровно L.
    u64 sum = 0;
    for (int i = 0; i < M; ++i) sum += f[i];
    assert(sum == L);
}

// ---------------------------------------------------------------------------
// Улучшение набора частот заменами пар с сохранением суммы и количества.
// На вход подаются частоты, уже являющиеся элементами SmoothSet.
// Замена (a,b) -> (a',b') при a'+b' = a+b сохраняет сумму и допустимость.
// ---------------------------------------------------------------------------
inline void smooth_refine(std::array<u32, M>& f,
                          const std::array<u32, M>& raw,
                          const SmoothSet& ss,
                          const std::vector<int>& syms,
                          int sweeps = 8) {
    const int n = static_cast<int>(syms.size());
    if (n < 2) return;

    // Предвычисленные логарифмы значений множества.
    std::vector<double> slog(ss.v.size());
    for (size_t i = 0; i < ss.v.size(); ++i) {
        slog[i] = std::log2(static_cast<double>(ss.v[i]));
    }

    std::vector<int> sidx(n);
    std::vector<u32> cur(n), cnt(n);
    for (int i = 0; i < n; ++i) {
        const int s = syms[i];
        cur[i] = f[s];
        cnt[i] = raw[s];
        sidx[i] = ss.index_of(cur[i]);
        assert(sidx[i] >= 0);
    }

    const int S = static_cast<int>(ss.v.size());
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        bool changed = false;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const u32 T = cur[i] + cur[j];
                const double target =
                    static_cast<double>(T) * cnt[i] / (static_cast<double>(cnt[i]) + cnt[j]);
                int k = ss.idx_le(static_cast<u32>(target));
                if (k < 0) k = 0;
                const double cur_cost = cnt[i] * slog[sidx[i]] + cnt[j] * slog[sidx[j]];
                double best_cost = cur_cost;
                int best_i = sidx[i], best_j = sidx[j];
                for (int d = -1; d <= 2; ++d) {
                    const int kk = k + d;
                    if (kk < 0 || kk >= S) continue;
                    const u32 a = ss.v[kk];
                    if (a >= T) continue;
                    const u32 b = T - a;
                    const int jb = ss.index_of(b);
                    if (jb < 0) continue;
                    const double cost = cnt[i] * slog[kk] + cnt[j] * slog[jb];
                    if (cost > best_cost + 1e-12) {
                        best_cost = cost;
                        best_i = kk;
                        best_j = jb;
                    }
                }
                if (best_i != sidx[i] || best_j != sidx[j]) {
                    f[syms[i]] = ss.v[best_i];
                    f[syms[j]] = ss.v[best_j];
                    sidx[i] = best_i;
                    sidx[j] = best_j;
                    cur[i] = ss.v[best_i];
                    cur[j] = ss.v[best_j];
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
}

// ---------------------------------------------------------------------------
// Жадное квантование: независимо выбираем ближайшее гладкое значение для
// каждой частоты, затем правим сумму одной или несколькими парными заменами.
// Идея: ближайшее гладкое значение даёт меньшую ошибку, чем pair-swap из Pow2,
// а сумма правится минимальным числом парных замен.
// ---------------------------------------------------------------------------
inline bool pair_change_sum(std::array<u32, M>& f,
                            const std::array<u32, M>& raw,
                            const SmoothSet& ss,
                            const std::vector<double>& slog,
                            long long delta) {
    if (delta == 0) return true;
    int n = 0;
    std::array<int, M> idx{};
    std::array<double, M> logf{};
    for (int i = 0; i < M; ++i) {
        if (f[i] > 0) {
            idx[n++] = i;
            logf[i] = std::log2(double(f[i]));
        }
    }
    if (n < 2) return false;

    const int S = static_cast<int>(ss.v.size());
    double best_gain = -1e300;
    int bi = -1, bj = -1;
    u32 ba = 0, bb = 0;
    for (int ii = 0; ii < n; ++ii) {
        const int i = idx[ii];
        for (int jj = ii + 1; jj < n; ++jj) {
            const int j = idx[jj];
            const long long T = (long long)f[i] + (long long)f[j] + delta;
            if (T < 2 || T > 2LL * 65536) continue;
            const double target = double(T) * raw[i] / double(raw[i] + raw[j]);
            int k = ss.idx_le(static_cast<u32>(std::min(target, 65536.0)));
            if (k < 0) k = 0;
            for (int d = -16; d <= 16; ++d) {
                const int kk = k + d;
                if (kk < 0 || kk >= S) continue;
                const u32 a = ss.v[kk];
                const long long b = T - a;
                if (b < 1 || b > 65536) continue;
                const int jb = ss.index_of(static_cast<u32>(b));
                if (jb < 0) continue;
                const double gain = raw[i] * (slog[kk] - logf[i]) +
                                    raw[j] * (slog[jb] - logf[j]);
                if (gain > best_gain) {
                    best_gain = gain;
                    bi = i; bj = j;
                    ba = a; bb = static_cast<u32>(b);
                }
            }
        }
    }
    if (bi < 0) return false;
    f[bi] = ba;
    f[bj] = bb;
    return true;
}

// Двухпроходная схема:
//   1) для каждой частоты берём наибольшее гладкое значение <= t_i,
//      где t_i = raw_i * L / N (нижняя граница брекета);
//   2) распределяем остаток R = L - sum(f) шагами вверх/вниз, выбирая символ
//      с максимальной маргинальной выгодой на единицу суммы (water-filling);
//   3) если остаток не удалось выбрать точно, правим его одной парной заменой.
// Возвращает false, если сумму не удалось свести к L.
inline bool waterfill_freqs(const std::array<u32, M>& raw, u32 L,
                            const SmoothSet& ss, std::array<u32, M>& f) {
    f.fill(0);
    u64 N64 = 0;
    int n = 0;
    std::array<int, M> idx{};
    for (int i = 0; i < M; ++i) {
        N64 += raw[i];
        if (raw[i] > 0) idx[n++] = i;
    }
    if (n == 0) return false;
    if (n == 1) { f[idx[0]] = L; return true; }

    const double N = double(N64);

    // Проход 1: нижняя граница брекета.
    long long sum = 0;
    for (int k = 0; k < n; ++k) {
        const int i = idx[k];
        const double t = double(raw[i]) * double(L) / N;
        int j = ss.idx_le(static_cast<u32>(std::floor(t)));
        if (j < 0) j = 0;
        f[i] = ss.v[j];
        sum += f[i];
    }

    // Проход 2: добираем/снимаем остаток по маргинальной выгоде.
    long long need = (long long)L - sum;
    const int S = static_cast<int>(ss.v.size());
    while (need != 0) {
        const bool up = need > 0;
        int best = -1;
        double best_prio = up ? -1e300 : 1e300;
        long long best_d = 0;
        for (int k = 0; k < n; ++k) {
            const int i = idx[k];
            const int j = ss.index_of(f[i]);
            if (j < 0) continue;
            const int jj = up ? j + 1 : j - 1;
            if (jj < 0 || jj >= S) continue;
            const u32 nv = ss.v[jj];
            const long long d = up ? (long long)nv - (long long)f[i]
                                   : (long long)f[i] - (long long)nv;
            if (d <= 0 || d > (up ? need : -need)) continue;
            const double gain =
                raw[i] * std::abs(std::log2(double(nv) / double(f[i])));
            const double prio = gain / double(d);
            if (up ? (prio > best_prio) : (prio < best_prio)) {
                best_prio = prio;
                best = k;
                best_d = d;
            }
        }
        if (best < 0) break;
        const int i = idx[best];
        const int j = ss.index_of(f[i]);
        f[i] = ss.v[up ? j + 1 : j - 1];
        need += up ? -best_d : best_d;
    }

    if (need == 0) return true;

    // Остаток (на практике 0..1) правим одной парной заменой.
    std::vector<double> slog(ss.v.size());
    for (size_t i = 0; i < ss.v.size(); ++i) slog[i] = std::log2(double(ss.v[i]));
    return pair_change_sum(f, raw, ss, slog, need);
}

// ---------------------------------------------------------------------------
// Основная точка входа: нормализация/квантование частот.
// raw - сырые счётчики (сумма N), L - целевая сумма (степень двойки).
// Результат f: сумма == L, f[i] > 0 там, где raw[i] > 0.
// ---------------------------------------------------------------------------
inline void quantize(const std::array<u32, M>& raw, u32 L, Mode mode,
                     std::array<u32, M>& f) {
    const int logL = cnt::int_log(L);
    assert((1u << logL) == L);

    std::vector<int> syms;
    syms.reserve(M);
    for (int i = 0; i < M; ++i) {
        if (raw[i] > 0) syms.push_back(i);
    }

    switch (mode) {
        case Mode::Precise:
            precise_freqs(raw, L, f);
            break;

        case Mode::Pow2:
            pow2_freqs(raw, L, logL, f);
            break;

        case Mode::Pow23:
        case Mode::Pow235: {
            const SmoothSet& ss = (mode == Mode::Pow23) ? smooth23() : smooth235();
            if (!waterfill_freqs(raw, L, ss, f)) {
                // Страховка: гарантированно корректный вариант.
                pow2_freqs(raw, L, logL, f);
                smooth_refine(f, raw, ss, syms, /*sweeps=*/8);
            }
            break;
        }
    }

    u64 sum = 0;
    for (int i = 0; i < M; ++i) sum += f[i];
    assert(sum == L && "quantized frequencies must sum to L");
}

// ---------------------------------------------------------------------------
// Компактная запись частот в заголовок: экспоненты в битовом потоке (LSB-first).
// ---------------------------------------------------------------------------
class BitWriter {
public:
    explicit BitWriter(u8* dst) : p(dst) {}
    void put(u32 value, int bits) {
        for (int i = 0; i < bits; ++i) {
            if ((bitpos & 7u) == 0) p[bitpos >> 3] = 0;
            const u32 bit = (value >> i) & 1u;
            p[bitpos >> 3] |= static_cast<u8>(bit << (bitpos & 7u));
            ++bitpos;
        }
    }
    u32 bytes() const { return (bitpos + 7) >> 3; }
private:
    u8* p;
    u32 bitpos = 0;
};

class BitReader {
public:
    explicit BitReader(const u8* src) : p(src) {}
    u32 get(int bits) {
        u32 v = 0;
        for (int i = 0; i < bits; ++i) {
            v |= static_cast<u32>((p[bitpos >> 3] >> (bitpos & 7u)) & 1u) << i;
            ++bitpos;
        }
        return v;
    }
    u32 bytes() const { return (bitpos + 7) >> 3; }
private:
    const u8* p;
    u32 bitpos = 0;
};

// Разложение f на простые множители 2, 3, 5: f == 2^a * 3^b * 5^c.
inline void pow_exponents(u32 f, u32& a, u32& b, u32& c) {
    a = 0; while ((f & 1u) == 0) { f >>= 1; ++a; }
    b = 0; while (f % 3 == 0) { f /= 3; ++b; }
    c = 0; while (f % 5 == 0) { f /= 5; ++c; }
    assert(f == 1 && "frequency is not a product of 2, 3, 5");
}

} // namespace qnt
