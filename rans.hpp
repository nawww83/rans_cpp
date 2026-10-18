#pragma once

#include "counter.hpp"
#include "io_utils.hpp"
#include "quantizer.hpp"

#include <algorithm>
#include <vector>

#if defined(_MSC_VER)
    // MSVC
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    // GCC и Clang
    #define RESTRICT __restrict__
#else
    // Остальные (или если не поддерживается)
    #define RESTRICT
#endif


namespace rans {

using u8 = cnt::u8;
using u16 = cnt::u16;
using u32 = cnt::u32;
using u64 = cnt::u64;

constexpr u32 Lmax = 65536;   // Doesn't change it!
constexpr u32 Bflush = 65536; // Doesn't change it!

static void assert_size(u32 in_size) {
    (void)in_size; // used only by asserts
    // Проверка на 0
    assert(in_size > 0 && "Size must be greater than zero");
    
    // Проверка на ограничение 2^31 (31-й бит должен быть 0)
    // 0x80000000U — это и есть 1 << 31
    assert((in_size & 0x80000000U) == 0 && "Size must be less than 2^31");
}

struct DecoderInfo {
    u16 f;        // Частота
    u16 cum_f;    // Накопленная частота
    u8 symbol;    // Символ-байт
};

struct FreqPair {
    u8 symbol;
    u32 count;
};

class Rans{
    public:
        explicit Rans() {}
        void encode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32& out_size,
                    qnt::Mode mode = qnt::Mode::Precise);
        void decode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32 out_size);
        u32 required_bytes(u32 in_size) const {
            assert_size(in_size);
            // 2 байта (M-1, mode|logL) + 5*u32 (4 состояния + размер) +
            // максимум 3 байта на символ (режим Precise).
            return in_size*2 + sizeof(u32)*5 + sizeof(u8)*2 + cnt::M*3;
        }
    private:
        io_u::io_utils io;
        cnt::Counter<Lmax> mCounter;
        std::array<u8, Lmax> mLut;
        std::array<DecoderInfo, cnt::M> dec_table;
};

inline void Rans::encode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32& out_size,
                         qnt::Mode mode) {
    assert_size(in_size);

    // 1. Считаем частоты и приводим их к выбранному режиму квантования.
    std::array<u32, cnt::M> freq;
    u32 L;
    if (mode == qnt::Mode::Precise) {
        // Точная нормализация (аналог арифметического кодирования).
        auto res = mCounter.count(input, in_size);
        freq = res.frequencies;
        L = res.Lscale;
    } else {
        const auto raw = cnt::Counter<Lmax>::histogram(input, in_size);
        L = cnt::Counter<Lmax>::choose_L(raw);
        qnt::quantize(raw, L, mode, freq);
    }
    const int logL = cnt::int_log(L);

    // 2. Группируем и сортируем (символ + его частота).
    struct SymFreq { u8 sym; u32 f; };
    std::array<SymFreq, cnt::M> sorted;
    u32 Ms = 0;
    for (int i = 0; i < cnt::M; ++i) {
        if (freq[i] > 0) sorted[Ms++] = {static_cast<u8>(i), freq[i]};
    }

    // Сортировка по убыванию частоты для кэш-локальности декодера
    std::sort(sorted.begin(), sorted.begin() + Ms, [](const SymFreq& a, const SymFreq& b) {
        return a.f > b.f;
    });

    // 3. Строим таблицы накопленных частот, маппинг и данные для спец-путей.
    std::array<u32, cnt::M> symbol_to_idx{};
    std::array<u32, cnt::M> fr{}; // Отсортированные частоты
    std::array<u32, cnt::M> cs{}; // Накопленные частоты
    std::array<u8,  cnt::M> sh{}; // Степень двойки в частоте
    std::array<u64, cnt::M> magic{}; // ceil(2^48 / odd): точное деление без коррекции
    std::array<u64, cnt::M> thr{}; // fr[i] * Bflush (порог ренормализации)
    std::array<u32, cnt::M> lmf{}; // L - fr[i]
    u32 cumulative = 0;
    for (u32 i = 0; i < Ms; ++i) {
        fr[i] = sorted[i].f;
        cs[i] = cumulative;
        cumulative += fr[i];
        symbol_to_idx[sorted[i].sym] = i;

        const u32 f = fr[i];
        int tz = 0;
        while (((f >> tz) & 1u) == 0) ++tz;
        sh[i] = static_cast<u8>(tz);
        const u32 odd = f >> tz;
        // y < odd * Bflush <= 2^32, поэтому (y * ceil(2^48/odd)) >> 48
        // равно floor(y / odd) без корректирующего умножения.
        magic[i] = ((u64(1) << 48) + odd - 1) / odd;
        thr[i] = static_cast<u64>(f) * Bflush;
        lmf[i] = L - f;
    }
    assert(cumulative == L);

    // 4. Записываем заголовок.
    //    Байт режима/масштаба: mode (2 бита) | logL (5 бит).
    //    Частоты квантованных режимов пишем экспонентами в битовом потоке:
    //      Precise: sym:8, f-1:16
    //      Pow2:    sym:8, a:5
    //      Pow23:   sym:8, a:5, b:4
    //      Pow235:  sym:8, a:5, b:4, c:3
    u8* out_ptr = output;
    *(out_ptr++) = static_cast<u8>(Ms - 1);
    *(out_ptr++) = static_cast<u8>((static_cast<u32>(mode) << 5) | static_cast<u32>(logL));
    {
        qnt::BitWriter bw(out_ptr);
        for (u32 i = 0; i < Ms; ++i) {
            bw.put(sorted[i].sym, 8);
            if (mode == qnt::Mode::Precise) {
                bw.put(fr[i] - 1, 16);
            } else {
                u32 a, b, c;
                qnt::pow_exponents(fr[i], a, b, c);
                bw.put(a, 5);
                if (mode != qnt::Mode::Pow2) bw.put(b, 4);
                if (mode == qnt::Mode::Pow235) bw.put(c, 3);
            }
        }
        out_ptr += bw.bytes();
    }

    // Резервируем место под финальные состояния и размер
    u8* const x_ptr = out_ptr; out_ptr += sizeof(u32);
    u8* const y_ptr = out_ptr; out_ptr += sizeof(u32);
    u8* const z_ptr = out_ptr; out_ptr += sizeof(u32);
    u8* const w_ptr = out_ptr; out_ptr += sizeof(u32);
    io.copy_to_mem<u32>(in_size, out_ptr);
    out_ptr += sizeof(u32);

    const u32 mask = Bflush - 1;

    // 5. Кодирование (Quad-Interleaved rANS). PATH выбирает способ деления на частоту:
    //    0 - обычное деление (Precise), 1 - сдвиг (Pow2), 2 - обратное умножение (Pow23/Pow235).
    auto run = [&]<int PATH>() {
        auto divf = [&](u32 v, u32 idx) -> u32 {
            if constexpr (PATH == 0) {
                return v / fr[idx];
            } else if constexpr (PATH == 1) {
                return v >> sh[idx];
            } else {
                const u32 y = v >> sh[idx];
                // Точное floor(y / odd): ошибка ceil-реципрока меньше 1/odd
                // на всём диапазоне y < odd * Bflush, коррекция не нужна.
                return static_cast<u32>((static_cast<u64>(y) * magic[idx]) >> 48);
            }
        };

        u32 x = L, y = L, z = L, w = L;

        u32 i = 0;
        for (; i + 3 < in_size; i += 4) {
            const u32 idx_x = symbol_to_idx[input[i]];
            const u32 idx_y = symbol_to_idx[input[i+1]];
            const u32 idx_z = symbol_to_idx[input[i+2]];
            const u32 idx_w = symbol_to_idx[input[i+3]];

            // Нормализация для всех четырех (выталкиваем в out_ptr).
            // thr[] = fr[] * Bflush и lmf[] = L - fr[] посчитаны заранее.
            if (x >= thr[idx_x]) { io.copy_to_mem<u16>(static_cast<u16>(x & mask), out_ptr); out_ptr += 2; x >>= 16; }
            if (y >= thr[idx_y]) { io.copy_to_mem<u16>(static_cast<u16>(y & mask), out_ptr); out_ptr += 2; y >>= 16; }
            if (z >= thr[idx_z]) { io.copy_to_mem<u16>(static_cast<u16>(z & mask), out_ptr); out_ptr += 2; z >>= 16; }
            if (w >= thr[idx_w]) { io.copy_to_mem<u16>(static_cast<u16>(w & mask), out_ptr); out_ptr += 2; w >>= 16; }

            // Непосредственно rANS кодирование
            x += divf(x, idx_x) * lmf[idx_x] + cs[idx_x];
            y += divf(y, idx_y) * lmf[idx_y] + cs[idx_y];
            z += divf(z, idx_z) * lmf[idx_z] + cs[idx_z];
            w += divf(w, idx_w) * lmf[idx_w] + cs[idx_w];
        }

        if (i < in_size) {
            const u32 idx_x = symbol_to_idx[input[i]];
            if (x >= thr[idx_x]) { io.copy_to_mem<u16>(static_cast<u16>(x & mask), out_ptr); out_ptr += 2; x >>= 16; }
            x += divf(x, idx_x) * lmf[idx_x] + cs[idx_x];
            i++;
        }
        if (i < in_size) {
            const u32 idx_y = symbol_to_idx[input[i]];
            if (y >= thr[idx_y]) { io.copy_to_mem<u16>(static_cast<u16>(y & mask), out_ptr); out_ptr += 2; y >>= 16; }
            y += divf(y, idx_y) * lmf[idx_y] + cs[idx_y];
            i++;
        }
        if (i < in_size) {
            const u32 idx_z = symbol_to_idx[input[i]];
            if (z >= thr[idx_z]) { io.copy_to_mem<u16>(static_cast<u16>(z & mask), out_ptr); out_ptr += 2; z >>= 16; }
            z += divf(z, idx_z) * lmf[idx_z] + cs[idx_z];
            i++;
        }

        // Сохраняем финальные состояния в заголовок
        io.copy_to_mem<u32>(x, x_ptr);
        io.copy_to_mem<u32>(y, y_ptr);
        io.copy_to_mem<u32>(z, z_ptr);
        io.copy_to_mem<u32>(w, w_ptr);

        out_size = static_cast<u32>(out_ptr - output);
    };

    switch (mode) {
        case qnt::Mode::Precise: run.template operator()<0>(); break;
        case qnt::Mode::Pow2:    run.template operator()<1>(); break;
        case qnt::Mode::Pow23:
        case qnt::Mode::Pow235:  run.template operator()<2>(); break;
    }
}

inline void Rans::decode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32 out_size) {
    (void)out_size; // used only by asserts
    const u8* in_ptr = input;
    // header: mode | logL, затем упакованная таблица частот
    int M = *(u8*)in_ptr; M += 1;
    in_ptr += sizeof(u8);
    const u8 mode_logL = *(in_ptr++);
    const u32 mode = mode_logL >> 5;   // 0=Precise, 1=Pow2, 2=Pow23, 3=Pow235
    const int logL = mode_logL & 0x1F;
    const u32 L = (1u << logL);
    assert( M <= cnt::M );
    u32 current_cum = 0;
    qnt::BitReader br(in_ptr);
    for (int i = 0; i < M; ++i) {
        const u8 sym = static_cast<u8>(br.get(8));
        u32 f;
        if (mode == 0) {
            f = br.get(16) + 1;
        } else {
            const u32 a = br.get(5);
            u32 b = 0, c = 0;
            if (mode >= 2) b = br.get(4);
            if (mode == 3) c = br.get(3);
            f = 1u << a;
            for (u32 t = 0; t < b; ++t) f *= 3;
            for (u32 t = 0; t < c; ++t) f *= 5;
        }

        // Заполняем таблицу
        dec_table[i] = { (u16)f, (u16)current_cum, sym };

        // LUT теперь указывает на индекс i (0, 1, 2...)
        for (u32 j = current_cum; j < current_cum + f; ++j) {
            mLut[j] = static_cast<u8>(i);
        }
        current_cum += f;
    }
    in_ptr += br.bytes();
    u32 y; u32 x; u32 z; u32 w;
    io.read_mem<u32>(x, in_ptr);
    in_ptr += sizeof(u32);
    io.read_mem<u32>(y, in_ptr);
    in_ptr += sizeof(u32);
    io.read_mem<u32>(z, in_ptr);
    in_ptr += sizeof(u32);
    io.read_mem<u32>(w, in_ptr);
    in_ptr += sizeof(u32);

    u32 decompressed_size;
    io.read_mem<u32>(decompressed_size, in_ptr);
    in_ptr += sizeof(u32);
    assert( decompressed_size == out_size );
    in_ptr = input + in_size;
    u8* out = output + decompressed_size;
    
    const u32 mask = L - 1;
    // Обработка остатка от деления на 4 (если есть)
    if ((size_t)(out - output) % 4 == 3) { 
        if (z < L) [[unlikely]] { in_ptr -= 2; z = (z << 16) | io.read_val<u16>(in_ptr); }
        const u8 iz = mLut[z & mask];
        const auto& inf_z = dec_table[iz];
        out--;
        out[0] = inf_z.symbol;
        z = (z >> logL) * inf_z.f + (z & mask) - inf_z.cum_f;
    }
    if ((size_t)(out - output) % 4 == 2) { 
        if (y < L) [[unlikely]] { in_ptr -= 2; y = (y << 16) | io.read_val<u16>(in_ptr); }
        const u8 iy = mLut[y & mask];
        const auto& inf_y = dec_table[iy];
        out--;
        out[0] = inf_y.symbol;
        y = (y >> logL) * inf_y.f + (y & mask) - inf_y.cum_f;
    }
    if ((size_t)(out - output) % 4 == 1) { 
        if (x < L) [[unlikely]] { in_ptr -= 2; x = (x << 16) | io.read_val<u16>(in_ptr); }
        const u8 ix = mLut[x & mask];
        const auto& inf_x = dec_table[ix];
        out--;
        out[0] = inf_x.symbol;
        x = (x >> logL) * inf_x.f + (x & mask) - inf_x.cum_f;
    }
    while (out > output) {
        // 1. Нормализация всех четырех (Branchless-style)
        if (w < L) [[unlikely]] { in_ptr -= 2; w = (w << 16) | io.read_val<u16>(in_ptr); }
        if (z < L) [[unlikely]] { in_ptr -= 2; z = (z << 16) | io.read_val<u16>(in_ptr); }
        if (y < L) [[unlikely]] { in_ptr -= 2; y = (y << 16) | io.read_val<u16>(in_ptr); }
        if (x < L) [[unlikely]] { in_ptr -= 2; x = (x << 16) | io.read_val<u16>(in_ptr); }

        // 2. Получаем индексы (Процессор может делать это параллельно)
        const u8 ix = mLut[x & mask];
        const u8 iy = mLut[y & mask];
        const u8 iz = mLut[z & mask];
        const u8 iw = mLut[w & mask];

        const auto& inf_x = dec_table[ix];
        const auto& inf_y = dec_table[iy];
        const auto& inf_z = dec_table[iz];
        const auto& inf_w = dec_table[iw];

        // 3. Запись 4 байт (LIFO порядок!)
        out -= 4;
        out[3] = inf_w.symbol;
        out[2] = inf_z.symbol;
        out[1] = inf_y.symbol;
        out[0] = inf_x.symbol;

        // 4. Шаги rANS (выполняются в суперскалярном режиме)
        w = (w >> logL) * inf_w.f + (w & mask) - inf_w.cum_f;
        z = (z >> logL) * inf_z.f + (z & mask) - inf_z.cum_f;
        y = (y >> logL) * inf_y.f + (y & mask) - inf_y.cum_f;
        x = (x >> logL) * inf_x.f + (x & mask) - inf_x.cum_f;
    }
}

} // namespace rans
