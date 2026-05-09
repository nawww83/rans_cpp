# rans_cpp
Byte-wise range ANS codec C++ implementation

# Features
- Quad interleave mode
- 16-bit flushing, 32-bit arithmetic
- Symbol frequency sorting
- 16-bit per frequency in the frequency table
- Auto scale selection (parameter L)
- Free-remainder rANS transformation: only "div" is used without "mod" function
- Fast symbol counter.

# Build
g++ main.cpp -O3 -std=c++20 -o rans

# Average performance (variate the Geometric distribution parameter)
4.5 GHz (desktop, Linux native, default turbo boost):
- Compression: 350 MB/s per core (77.7 MB/s/GHz)
- Decompression: 510 MB/s per core (113.3 MB/s/GHz).

1.99 Ghz (laptop, WSL Ubuntu, switch off turbo boost):
- Compression: 120 MB/s per core (60 MB/s/GHz)
- Decompression: 195 MB/s per core (98 MB/s/GHz).

Under the conditions:
- The input size is in [2, 8] MB
- Compilers:
  * GCC 11.4 (desktop)
  * GCC 12.5 (WSL Ubuntu).
