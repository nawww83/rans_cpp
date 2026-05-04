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

    // Универсальный bswap для целых чисел
    template<typename T>
    static constexpr T byteswap(T val) {
        static_assert(std::is_integral_v<T>, "Only integers supported");
        if constexpr (sizeof(T) == 1) return val;
        else if constexpr (sizeof(T) == 2) return __builtin_bswap16(static_cast<uint16_t>(val));
        else if constexpr (sizeof(T) == 4) return __builtin_bswap32(static_cast<uint32_t>(val));
        else if constexpr (sizeof(T) == 8) return __builtin_bswap64(static_cast<uint64_t>(val));
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
