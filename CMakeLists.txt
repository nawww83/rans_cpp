cmake_minimum_required(VERSION 3.16)

project(rans_cpp
    VERSION 1.0.0
    DESCRIPTION "Byte-wise rANS codec with selectable frequency quantization"
    LANGUAGES CXX)

# --- C++20 ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release for single-config generators.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

set(RANS_WARNING_FLAGS "")
if(MSVC)
    # /utf-8          - source files contain UTF-8 comments
    # /Zc:__cplusplus - report the real __cplusplus value
    # /permissive-    - conforming mode
    set(RANS_WARNING_FLAGS /W4 /utf-8 /Zc:__cplusplus /permissive-)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(RANS_WARNING_FLAGS -Wall -Wextra)
endif()

# --- Benchmark / demo executable ---
add_executable(rans main.cpp)
target_include_directories(rans PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(rans PRIVATE ${RANS_WARNING_FLAGS})

# --- Tests ---
add_executable(rans_tests test_quantizer.cpp)
target_include_directories(rans_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(rans_tests PRIVATE ${RANS_WARNING_FLAGS})

# The tests are assert-based, so keep asserts enabled even in Release builds.
if(MSVC)
    target_compile_options(rans_tests PRIVATE /UNDEBUG)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(rans_tests PRIVATE -UNDEBUG)
endif()

enable_testing()
add_test(NAME quantizer COMMAND rans_tests)

# --- Optional native tuning ---
option(RANS_NATIVE "Enable CPU-specific tuning (-march=native or /arch:AVX2)" OFF)
if(RANS_NATIVE)
    if(MSVC)
        target_compile_options(rans PRIVATE /arch:AVX2)
        target_compile_options(rans_tests PRIVATE /arch:AVX2)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(rans PRIVATE -march=native)
        target_compile_options(rans_tests PRIVATE -march=native)
    endif()
endif()
