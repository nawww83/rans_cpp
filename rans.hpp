#pragma once

#include "counter.hpp"
#include <cassert>
#include <algorithm>


namespace rans{

using u8 = cnt::u8;
using u16 = cnt::u16;
using u32 = cnt::u32;
using u64 = cnt::u64;

constexpr u32 L = 65536;

static void assert_size(int in_size) {
    assert( (u64)in_size < (1ull << (8*sizeof(u32) - 1)) ); // N < 2^31;
}

class Rans{
    public:
        explicit Rans() = default;
        void encode(const u8* input, int in_size, u8* output, int& out_size);
        void decode(const u8* input, int in_size, u8* output, int out_size);
        int required_bytes(int in_size) const {
            assert_size(in_size);
            return in_size*2 + sizeof(u32)*2 + sizeof(u8)*1 + cnt::M*(sizeof(u16) + sizeof(u8));
        }
    private:
        cnt::Counter<L> c;
        std::array<u32, cnt::M+1> cs;
};

void Rans::encode(const u8* input, int in_size, u8* output, int& out_size) {
    assert_size(in_size);
    assert( (1ull << (8*sizeof(u16)) == (u64)L));
    auto fr { c.count(input, in_size) };
    //
    std::array<u8, cnt::M> idx_to_symbol;
    for (int i=0; i<cnt::M; ++i) {
        idx_to_symbol[i] = i;
    }
    //
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
    int M = 0;
    for (int i=0; i<cnt::M; ++i) {
        cs[i+1] = cs[i] + fr[i];
        M += (fr[i] > 0);
    }
    assert(M <= cnt::M);
    u32 x = 1;
    const u8* out_ptr = output;
    // header: frequency table in 16-bit format
    // std::cout << " write frequency tables: M: " << M << ": ";
    *(u8*)out_ptr = M-1;
    out_ptr += sizeof(u8);
    std::array<u8, cnt::M> symbol_to_idx;
    for (int i=0; i<M; ++i) {
        symbol_to_idx[idx_to_symbol[i]] = i;
        *(u8*)out_ptr = (u8)idx_to_symbol[i];
        out_ptr += sizeof(u8);
        *(u16*)out_ptr = (u16)(fr[i] - 1);
        out_ptr += sizeof(u16);
        // std::cout << " write sym: " << idx_to_symbol[i] << ", f: " << fr[i] << "; ";
        // std::cout << (u16)fr.at(i) << ", ";
    }
    // std::cout << " M: " << M << std::endl;
    const u8* x_ptr = out_ptr;
    out_ptr += sizeof(u32);
    *(u32*)out_ptr = in_size;
    out_ptr += sizeof(u32);
    // std::cout << " write in size: " << in_size << std::endl;
    //
    for (auto i=0; i<in_size; ++i) {
        int idx = symbol_to_idx[input[i]];
        const u32 f = fr[idx];
        assert(f > 0);
        auto update_x = [&x, f](u32 cum_fr) { // classic rANS transformation: but we eliminate (x mod f) term
            x += ( (x / f) * (L - f) + cum_fr );
        };
        if ((x >> (8*sizeof(u16))) < f) { // most likely if there is a good compression
            update_x(cs[idx]);
            continue;
        }
        {   // flush output
            *(u16*)out_ptr = x & (L - 1);
            // std::cout << " flush: " << (x & (base - 1)) << std::endl;
            out_ptr += sizeof(u16);
            x >>= 8*sizeof(u16);
        }
        // x += ( (x / f) * (L - f) + cs[input[i]] );
        update_x(cs[idx]);
    }
    *(u32*)x_ptr = x;
    // std::cout << " write x: " << x << std::endl;
    out_size = out_ptr - output;
}

void Rans::decode(const u8* input, int in_size, u8* output, int out_size) {
    const u8* in_ptr = input;
    cs[0] = 0;
    // header: frequency table in 16-bit format
    int M = *(u8*)in_ptr; M += 1;
    in_ptr += sizeof(u8);
    // std::cout << " M: " << M << std::endl;
    // std::cout << " read frequency tables: M: " << M << ": ";
    assert( M <= cnt::M );
    std::array<u8, cnt::M> symbols;
    {
        int i = 0;
        while (i < M) {   
            symbols[i] = *(u8*)in_ptr;
            in_ptr += sizeof(u8);
            u32 f = *(u16*)in_ptr;
            f += 1;
            in_ptr += sizeof(u16);
            // std::cout << " read sym: " << int(symbols[i]) << ", f: " << f << "; ";
            cs[i+1] = cs[i] + f;
            i++;
        }
    }
    // std::cout << std::endl;
    u32 x = *(u32*)in_ptr;
    in_ptr += sizeof(u32);
    // std::cout << " read x: " << x << std::endl;
    int decompressed_size = *(u32*)in_ptr;
    in_ptr += sizeof(u32);
    // std::cout << " read decomp size: " << decompressed_size << std::endl;
    assert( decompressed_size == out_size );
    //
    const u8* data_ptr = in_ptr;
    in_ptr = input + (in_size - sizeof(u16));
    int i = 1;
    while (i <= decompressed_size) {
        u32 xmod = x & (L - 1);
        u32 xdiv = x >> 8*sizeof(u16);
        u32 idx = 0;     
        while (cs[idx+1] <= xmod) {
            idx++;
        }
        // idx = std::upper_bound(cs.begin(), cs.begin() + M, xmod) - cs.begin() - 1;      
        assert(idx < cnt::M);
        output[decompressed_size - i] = symbols[idx]; i++;
        x = xdiv*(cs[idx+1] - cs[idx]) + xmod - cs[idx];
        if ((in_ptr >= data_ptr) && (x < L)) {
            x <<= 8*sizeof(u16);
            // std::cout << " read stream: " << *(u16*)in_ptr << std::endl;
            x |= *(u16*)in_ptr;
            in_ptr -= sizeof(u16);
        }
    }
}

} // namespace rans
