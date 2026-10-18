#pragma once

#include <cstdint>
#include <bit>
#include <cstring>
#include <type_traits>

namespace io_u {

class io_utils {
public:
    // В C++20 используем std::endian, для более старых версий - простая проверка
    static constexpr bool is_little_endian() {
#if __cplusplus >= 202002L
        return std::endian::native == std::endian::little;
#else
        return uint16_t(0x0001) == 0x0001; 
#endif
    }

    static constexpr bool is_big_endian() {
#if __cplusplus >= 202002L
        return std::endian::native == std::endian::big;
#else
        return !is_little_endian(); 
#endif
    }

    // Универсальный bswap для целых чисел (портируемо, без интринсиков компилятора).
    template<typename T>
    static constexpr T byteswap(T val) {
        static_assert(std::is_integral_v<T>, "Only integers supported");
        if constexpr (sizeof(T) == 1) {
            return val;
        } else if constexpr (sizeof(T) == 2) {
            const uint16_t v = static_cast<uint16_t>(val);
            return static_cast<T>((v << 8) | (v >> 8));
        } else if constexpr (sizeof(T) == 4) {
            const uint32_t v = static_cast<uint32_t>(val);
            return static_cast<T>(((v & 0x000000FFu) << 24) |
                                  ((v & 0x0000FF00u) << 8)  |
                                  ((v & 0x00FF0000u) >> 8)  |
                                  ((v & 0xFF000000u) >> 24));
        } else if constexpr (sizeof(T) == 8) {
            const uint64_t v = static_cast<uint64_t>(val);
            return static_cast<T>(((v & 0x00000000000000FFull) << 56) |
                                  ((v & 0x000000000000FF00ull) << 40) |
                                  ((v & 0x0000000000FF0000ull) << 24) |
                                  ((v & 0x00000000FF000000ull) << 8)  |
                                  ((v & 0x000000FF00000000ull) >> 8)  |
                                  ((v & 0x0000FF0000000000ull) >> 24) |
                                  ((v & 0x00FF000000000000ull) >> 40) |
                                  ((v & 0xFF00000000000000ull) >> 56));
        } else {
            static_assert(sizeof(T) <= 8, "Unsupported integer size");
        }
    }

    template<typename T>
    static void copy_to_mem(T x, uint8_t* buffer) {
        if constexpr (!is_little_endian()) {
            x = byteswap(x);
        }
        std::memcpy(buffer, &x, sizeof(T));
    }

    template<typename T>
    static void read_mem(T& x, const uint8_t* buffer) {
        std::memcpy(&x, buffer, sizeof(T));
        if constexpr (!is_little_endian()) {
            x = byteswap(x);
        }
    }

    // "Быстрая" версия внутри io_utils специально для горячих циклов
    template<typename T>
    [[nodiscard]] static inline T read_val(const uint8_t* buffer) {
        T x;
        std::memcpy(&x, buffer, sizeof(T));
        if constexpr (!is_little_endian()) {
            x = byteswap(x);
        }
        return x;
    }

};

} // namespace io_u
