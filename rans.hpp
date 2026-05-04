#pragma once

#include "counter.hpp"
#include "io_utils.hpp"

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
        void encode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32& out_size);
        void decode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32 out_size);
        u32 required_bytes(u32 in_size) const {
            assert_size(in_size);
            return in_size*2 + sizeof(u32)*3 + sizeof(u8)*2 + cnt::M*(sizeof(u16) + sizeof(u8));
        }
    private:
        io_u::io_utils io;
        cnt::Counter<Lmax> mCounter;
        std::array<u8, Lmax> mLut;
        std::array<DecoderInfo, cnt::M> dec_table;
};

inline void Rans::encode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32& out_size) {
    assert_size(in_size);
    
    // 1. Считаем частоты (используем ваш счетчик)
    auto [raw_fr, L] = mCounter.count(input, in_size);
    const int logL = cnt::int_log(L);

    // 2. Группируем и сортируем (символ + его частота)
    struct FreqPair { u8 sym; u32 f; };
    std::vector<FreqPair> sorted;
    sorted.reserve(256);
    for (int i = 0; i < 256; ++i) {
        if (raw_fr[i] > 0) sorted.push_back({(u8)i, raw_fr[i]});
    }

    // Сортировка по убыванию частоты для кэш-локальности декодера
    std::sort(sorted.begin(), sorted.end(), [](auto& a, auto& b) {
        return a.f > b.f;
    });

    const u32 M = sorted.size();
    
    // 3. Строим таблицы накопленных частот и маппинг для кодирования
    std::array<u32, 256> symbol_to_idx;
    std::array<u32, 256> fr; // Отсортированные частоты
    std::array<u32, 257> cs;
    cs[0] = 0;
    for (u32 i = 0; i < M; ++i) {
        fr[i] = sorted[i].f;
        cs[i+1] = cs[i] + fr[i];
        symbol_to_idx[sorted[i].sym] = i; // Какой индекс в таблице у байта
    }

    // 4. Записываем заголовок
    u8* out_ptr = output;
    *(out_ptr++) = static_cast<u8>(M - 1);
    *(out_ptr++) = static_cast<u8>(logL);

    for (u32 i = 0; i < M; ++i) {
        *(out_ptr++) = sorted[i].sym;
        io.copy_to_mem<u16>(static_cast<u16>(fr[i] - 1), out_ptr);
        out_ptr += sizeof(u16);
    }

    // Резервируем место под финальные состояния и размер
    u8* const x_ptr = out_ptr; out_ptr += sizeof(u32);
    u8* const y_ptr = out_ptr; out_ptr += sizeof(u32);
    u8* const z_ptr = out_ptr; out_ptr += sizeof(u32);
    u8* const w_ptr = out_ptr; out_ptr += sizeof(u32);
    io.copy_to_mem<u32>(in_size, out_ptr);
    out_ptr += sizeof(u32);

    // 5. Кодирование (Quad-Interleaved rANS)
    // Инициализируем 4 состояния
    u32 x = L, y = L, z = L, w = L;

    const u32 mask = Bflush - 1;

    size_t i = 0;
    for (; i + 3 < in_size; i += 4) {
        u32 idx_x = symbol_to_idx[input[i]];
        u32 idx_y = symbol_to_idx[input[i+1]];
        u32 idx_z = symbol_to_idx[input[i+2]];
        u32 idx_w = symbol_to_idx[input[i+3]];

        // Нормализация для всех четырех (выталкиваем в out_ptr)
        if (x >= fr[idx_x] * Bflush) { io.copy_to_mem<u16>(x & mask, out_ptr); out_ptr += 2; x >>= 16; }
        if (y >= fr[idx_y] * Bflush) { io.copy_to_mem<u16>(y & mask, out_ptr); out_ptr += 2; y >>= 16; }
        if (z >= fr[idx_z] * Bflush) { io.copy_to_mem<u16>(z & mask, out_ptr); out_ptr += 2; z >>= 16; }
        if (w >= fr[idx_w] * Bflush) { io.copy_to_mem<u16>(w & mask, out_ptr); out_ptr += 2; w >>= 16; }

        // Непосредственно rANS кодирование
        x += (x / fr[idx_x]) * (L - fr[idx_x]) + cs[idx_x];
        y += (y / fr[idx_y]) * (L - fr[idx_y]) + cs[idx_y];
        z += (z / fr[idx_z]) * (L - fr[idx_z]) + cs[idx_z];
        w += (w / fr[idx_w]) * (L - fr[idx_w]) + cs[idx_w];
    }

    if (i < in_size) {
        u32 idx_x = symbol_to_idx[input[i]];
        // Нормализация
        if (x >= fr[idx_x] * Bflush) { io.copy_to_mem<u16>(x & mask, out_ptr); out_ptr += 2; x >>= 16; }
        // Непосредственно rANS кодирование
        x += (x / fr[idx_x]) * (L - fr[idx_x]) + cs[idx_x];
        i++;
    }
    if (i < in_size) {
        u32 idx_y = symbol_to_idx[input[i]];
        // Нормализация
        if (y >= fr[idx_y] * Bflush) { io.copy_to_mem<u16>(y & mask, out_ptr); out_ptr += 2; y >>= 16; }
        // Непосредственно rANS кодирование
        y += (y / fr[idx_y]) * (L - fr[idx_y]) + cs[idx_y];
        i++;
    }
    if (i < in_size) {
        u32 idx_z = symbol_to_idx[input[i]];
        // Нормализация
        if (z >= fr[idx_z] * Bflush) { io.copy_to_mem<u16>(z & mask, out_ptr); out_ptr += 2; z >>= 16; }
        // Непосредственно rANS кодирование
        z += (z / fr[idx_z]) * (L - fr[idx_z]) + cs[idx_z];
        i++;
    }

    // Сохраняем финальные состояния в заголовок
    io.copy_to_mem<u32>(x, x_ptr);
    io.copy_to_mem<u32>(y, y_ptr);
    io.copy_to_mem<u32>(z, z_ptr);
    io.copy_to_mem<u32>(w, w_ptr);

    out_size = out_ptr - output;
}

inline void Rans::decode(const u8* RESTRICT input, u32 in_size, u8* RESTRICT output, u32 out_size) {
    const u8* in_ptr = input;
    // header: frequency table in 16-bit format
    int M = *(u8*)in_ptr; M += 1;
    in_ptr += sizeof(u8);
    const int logL = *(u8*)in_ptr;
    in_ptr += sizeof(u8);
    const u32 L = (1u << logL);
    assert( M <= cnt::M );
    u32 current_cum = 0;
    for (int i = 0; i < M; ++i) {
        u8 sym = *(in_ptr++);
        u16 f0;
        io.read_mem<u16>(f0, in_ptr);
        in_ptr += sizeof(u16);
        u32 f = static_cast<u32>(f0) + 1;

        // Заполняем таблицу
        dec_table[i] = { (u16)f, (u16)current_cum, sym };

        // LUT теперь указывает на индекс i (0, 1, 2...)
        for (u32 j = current_cum; j < current_cum + f; ++j) {
            mLut[j] = static_cast<u8>(i);
        }
        current_cum += f;
    }
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
    const u8* data_ptr = in_ptr;
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
