# rans_cpp

Byte-wise range ANS codec implementation in C++.

A compact and fast entropy coder designed for high-throughput compression and decompression workloads. It emphasizes efficient integer arithmetic, low-overhead symbol handling, and strong performance across practical workloads.

## Why rans_cpp?

- High throughput on modern CPUs
- Byte-oriented range ANS design
- Optimized for 32-bit arithmetic and compact frequency tables
- Fast symbol counting and frequency sorting
- Efficient implementation with a focus on performance

## Features

- Quad interleave mode
- 16-bit flushing and 32-bit arithmetic
- Symbol frequency sorting
- 16-bit per-frequency entries in the frequency table
- Automatic scale selection via parameter L
- Free-remainder rANS transformation using only `div` (no `mod` required)
- Fast symbol counter

## Build

```bash
g++ main.cpp -O3 -std=c++20 -o rans
```

## Quick usage

```bash
./rans
```

The project is intended as a C++ implementation of byte-wise range ANS compression and decompression routines, suitable for benchmarking and experimentation.

## Performance

Benchmark results on representative workloads:

### Desktop (4.5 GHz, Linux native, default turbo boost)
- Intel Core i7 12th Gen
- Compression: 350 MB/s per core (77.7 MB/s/GHz)
- Decompression: 510 MB/s per core (113.3 MB/s/GHz)

### Laptop (1.99 GHz, WSL Ubuntu, turbo boost disabled)
- Intel Core i7 8th Gen
- Compression: 120 MB/s per core (60 MB/s/GHz)
- Decompression: 195 MB/s per core (98 MB/s/GHz)

### Conditions
- Input size: 2–8 MB
- Compiler: GCC 11.4 on desktop, GCC 12.5 on WSL Ubuntu

## Notes

This implementation is optimized for speed and is designed for experimentation and performance-oriented use cases. The code focuses on efficient range ANS execution while keeping the design compact and approachable.
