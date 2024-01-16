#pragma once

#include "counter.hpp"
#include "io_utils.hpp"

#include <algorithm>
#include <vector>

namespace rans {

using u8 = cnt::u8;
using u16 = cnt::u16;
using u32 = cnt::u32;
using u64 = cnt::u64;

constexpr u32 Lmax = 65536;   // Doesn't change it!
constexpr u32 Bflush = 65536; // Doesn't change it!


static void assert_size(u32 in_size) {
    assert( in_size > 0 );
    assert( (u64)in_size < (1ull << (8*sizeof(u32) - 1)) ); // N < 2^31;
}

class Rans{
    public:
        explicit Rans() {
            lut.resize(Lmax);
        }
        void encode(const u8* input, u32 in_size, u8* output, u32& out_size);
        void decode(const u8* input, u32 in_size, u8* output, u32 out_size);
        u32 required_bytes(u32 in_size) const {
            assert_size(in_size);
            return in_size*2 + sizeof(u32)*3 + sizeof(u8)*2 + cnt::M*(sizeof(u16) + sizeof(u8));
        }
    private:
        io_u::io_utils io;
        cnt::Counter<Lmax> c;
        std::array<u32, cnt::M+1> cs;
        std::vector<u8> lut;
};

void Rans::encode(const u8* input, u32 in_size, u8* output, u32& out_size) {
    assert_size(in_size);
    assert( (1ull << (8*sizeof(u16)) == (u64)Bflush) );
    auto fr { c.count(input, in_size) };
    //
    std::array<u32, cnt::M> idx_to_symbol;
    for (int i=0; i<cnt::M; ++i) {
        idx_to_symbol[i] = i;
    }
    // Frequency sorting: O(N^2) but the size N is limited by 256
    for (int i=0; i<cnt::M; ++i) {
        for (int j=i; j<cnt::M; ++j) {
            if (fr[i] < fr[j]) {
                std::swap(fr[i], fr[j]);
                std::swap(idx_to_symbol[i], idx_to_symbol[j]);
            }
        }
    }
    //
    cs[0] = 0;
    int M = 0; // Alphabet length calculation
    for (int i=0; i<cnt::M; ++i) {
        cs[i+1] = cs[i] + fr[i];
        M += (fr[i] > 0);
    }
    assert(M <= cnt::M);
    const u32 L = cs[M];
    const int logL = cnt::int_log(L);
    assert( (int)L >= M );
    // std::cout << "L encode: " << L << std::endl;
    u8* out_ptr = output;
    // header: frequency table in 16-bit format
    // std::cout << " write frequency tables: M: " << M << ": ";
    *(u8*)out_ptr = M-1;
    out_ptr += sizeof(u8);
    *(u8*)out_ptr = (u8)logL;
    out_ptr += sizeof(u8);
    std::array<u32, cnt::M> symbol_to_idx;
    for (size_t i=0; i<M; ++i) {
        symbol_to_idx[idx_to_symbol[i]] = i;
        *(u8*)out_ptr = (u8)idx_to_symbol[i];
        out_ptr += sizeof(u8);
        // *(u16*)out_ptr = (u16)(fr[i] - 1);
        io.copy_to_mem_16(static_cast<u16>(fr[i]-1), out_ptr, sizeof(u16));
        out_ptr += sizeof(u16);
        // std::cout << " write sym: " << idx_to_symbol[i] << ", f: " << fr[i] << "; ";
        // std::cout << (u16)fr[i] << ", ";
    }
    // std::cout << " M: " << M << std::endl;
    u8* const x_ptr = out_ptr;
    out_ptr += sizeof(u32);
    u8* const y_ptr = out_ptr;
    out_ptr += sizeof(u32);
    // *(u32*)out_ptr = in_size;
    io.copy_to_mem_32(in_size, out_ptr, sizeof(u32));
    out_ptr += sizeof(u32);
    // std::cout << " write in size: " << in_size << std::endl;
    // x = [f, B*f) <=> [L, L*B); B - the parameter which determines output flushing
    u32 x = (in_size > 1) ? fr[symbol_to_idx[input[0]]] : 1;
    u32 y = (in_size > 1) ? fr[symbol_to_idx[input[1]]] : 1;
    size_t i = 0;
    while (i < in_size/2) {
        // std::cout << " encode: step i: " << i << std::endl;
        const int idx1 = symbol_to_idx[input[2*i]];
        const int idx2 = symbol_to_idx[input[2*i+1]];
        const u32 f1 = fr[idx1];
        const u32 f2 = fr[idx2];
        // assert(f > 0);
        auto update_xy = [&x, &y, f1, f2, L, this](u32 cum_fr1, u32 cum_fr2) { // classic rANS transformation: but we eliminate (x mod f) term
            x += ( (x / f1) * (L - f1) + cum_fr1 );
            y += ( (y / f2) * (L - f2) + cum_fr2 );
        };
        if ((u64)x >= (u64)(f1)*Bflush) {   // flush output
            // *(u16*)out_ptr = x & (Bflush - 1);
            io.copy_to_mem_16(static_cast<u16>(x & (Bflush - 1)), out_ptr, sizeof(u16));
            out_ptr += sizeof(u16);
            x >>= 8*sizeof(u16);
            // std::cout << " after flush: x: " << x << std::endl;
        }
        if ((u64)y >= (u64)(f2)*Bflush) {   // flush output
            // *(u16*)out_ptr = y & (Bflush - 1);
            io.copy_to_mem_16(static_cast<u16>(y & (Bflush - 1)), out_ptr, sizeof(u16));
            out_ptr += sizeof(u16);
            y >>= 8*sizeof(u16);
            //  std::cout << " after flush: y: " << y << std::endl;
        }
        // x += ( (x / f) * (L - f) + cs[input[i]] );
        update_xy(cs[idx1], cs[idx2]);
        // std::cout << "x: " << x << ", y: " << y << std::endl;
        i += 1;
    }
    if ((in_size % 2) == 1) {
        *(u8*)out_ptr = input[in_size-1];
        out_ptr += sizeof(u8);
    }
    // *(u32*)x_ptr = x;
    // *(u32*)y_ptr = y;
    io.copy_to_mem_32(x, x_ptr, sizeof(u32));
    io.copy_to_mem_32(y, y_ptr, sizeof(u32));
    // std::cout << " write x: " << x << std::endl;
    // std::cout << " write y: " << y << std::endl;
    out_size = out_ptr - output;
}

void Rans::decode(const u8* input, u32 in_size, u8* output, u32 out_size) {
    const u8* in_ptr = input;
    cs[0] = 0;
    // header: frequency table in 16-bit format
    int M = *(u8*)in_ptr; M += 1;
    in_ptr += sizeof(u8);
    const int logL = *(u8*)in_ptr;
    in_ptr += sizeof(u8);
    const u32 L = (1u << logL);
    // std::cout << "L decode: " << L << std::endl;
    // std::cout << " M: " << M << std::endl;
    // std::cout << " read frequency tables: M: " << M << ": ";
    assert( M <= cnt::M );
    std::array<u32, cnt::M> symbols;
    {
        bool use_lut = (M >= 8);
        lut.resize(use_lut ? L : 0);
    }
    {
        int i = -1;
        while (++i < M) {   
            symbols[i] = *(u8*)in_ptr;
            in_ptr += sizeof(u8);
            // const u32 f = (u32)(*(u16*)in_ptr) + 1;
            u16 f0;
            io.read_mem_16(f0, in_ptr, sizeof(u16));
            in_ptr += sizeof(u16);
            u32 f = static_cast<u32>(f0) + 1;
            // std::cout << " read sym: " << int(symbols[i]) << ", f: " << f << "; ";
            cs[i+1] = cs[i] + f;
            if (! lut.empty()) {
                for (u32 j=cs[i]; j<cs[i+1]; ++j) {
                    lut[j] = i & 255;
                }
            }
        }
    }
    // std::cout << std::endl;
    // u32 x = *(u32*)in_ptr;
    u32 x;
    io.read_mem_32(x, in_ptr, sizeof(u32));
    in_ptr += sizeof(u32);
    // u32 y = *(u32*)in_ptr;
    u32 y;
    io.read_mem_32(y, in_ptr, sizeof(u32));
    in_ptr += sizeof(u32);
    // std::cout << " read x: " << x << std::endl;
    // std::cout << " read y: " << y << std::endl;
    // const int decompressed_size = *(u32*)in_ptr;
    u32 decompressed_size;
    io.read_mem_32(decompressed_size, in_ptr, sizeof(u32));
    in_ptr += sizeof(u32);
    // std::cout << " read decomp size: " << decompressed_size << std::endl;
    assert( decompressed_size == out_size );
    //
    const u8* data_ptr = in_ptr;
    in_ptr = input + in_size;
    auto read_xy = [&x, &y, &in_ptr, data_ptr, logL, this]() {
        bool y_is_good = ((y >> logL) != 0); // y >= L
        bool x_is_good = ((x >> logL) != 0); // x >= L
        bool is_two_enough = (in_ptr >= (data_ptr + 2*sizeof(u16)));
        bool is_one_enough = (in_ptr >= (data_ptr + 1*sizeof(u16)));
        // std::cout << "is two enough: " << is_two_enough << ", is one enough: " << is_one_enough << std::endl;
        if ((!y_is_good) && (!x_is_good) && is_two_enough) {
            in_ptr -= sizeof(u16);
            y <<= 8*sizeof(u16);
            // y |= *(u16*)in_ptr;
            u16 tmp;
            io.read_mem_16(tmp, in_ptr, sizeof(u16));
            y |= tmp;
            in_ptr -= sizeof(u16);
            // std::cout << " read stream: y: " << y << std::endl;
            x <<= 8*sizeof(u16);
            // x |= *(u16*)in_ptr;
            io.read_mem_16(tmp, in_ptr, sizeof(u16));
            x |= tmp;
            // std::cout << " read stream: x: " << x << std::endl;
            return;
        }
        if ((!y_is_good) && (x_is_good) && is_one_enough) {
            in_ptr -= sizeof(u16);
            y <<= 8*sizeof(u16);
            // y |= *(u16*)in_ptr;
            u16 tmp;
            io.read_mem_16(tmp, in_ptr, sizeof(u16));
            y |= tmp;
            // std::cout << " read stream: only y: " << y << std::endl;
            return;
        }
        if ((y_is_good) && (!x_is_good) && is_one_enough) {
            in_ptr -= sizeof(u16);
            x <<= 8*sizeof(u16);
            // x |= *(u16*)in_ptr;
            u16 tmp;
            io.read_mem_16(tmp, in_ptr, sizeof(u16));
            x |= tmp;
            // std::cout << " read stream: only x: " << x << std::endl;
            return;
        }
    };
    u8* out = output + decompressed_size;
    if ((decompressed_size % 2) == 1) {
        in_ptr -= sizeof(u8);
        *(--out) = *in_ptr;
    }
    auto get_idx = [this](u32 xmod1, u32 xmod2) {
        std::pair<int, int> idxs {0, 0};
        if (! lut.empty()) {
            idxs = {lut[xmod1], lut[xmod2]};
            return idxs;   
        }
        bool f1 = cs[idxs.first + 1] <= xmod1 ;
        bool f2 = cs[idxs.second + 1] <= xmod2 ;
        while (f1 || f2) {
            idxs.first = f1 ? idxs.first + 1 : idxs.first;
            idxs.second = f2 ? idxs.second + 1 : idxs.second;
            f1 = (f1) ? cs[idxs.first + 1] <= xmod1 : f1;
            f2 = (f2) ? cs[idxs.second + 1] <= xmod2 : f2;
        }
        return idxs;
    };
    while (out != output) {
        // Linear search: symbol frequencies are sorted when they were encoded
        auto idxs = get_idx(x & (L - 1), y & (L - 1));   
        assert( idxs.first != -1 ); assert( idxs.second != -1 );
        *(--out) = symbols[idxs.second];
        *(--out) = symbols[idxs.first];
        x = (x >> logL) * (cs[idxs.first + 1] - cs[idxs.first]) + ((x & (L - 1)) - cs[idxs.first]);
        y = (y >> logL) * (cs[idxs.second + 1] - cs[idxs.second]) + ((y & (L - 1)) - cs[idxs.second]);
        // std::cout << " x: " << x << ", y: " << y << std::endl;
        read_xy();
    }   
}

} // namespace rans
